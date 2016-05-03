#include "tools.h"
#include "script.h"

#include <QCoreApplication>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <QDebug>

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include <queue>
#include <map>
#include <vector>
#include <thread>
#include <mutex>

#include <QtSql>

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#pragma GCC diagnostic pop
using namespace Eigen;

const double PI = M_PI;

using namespace std;

// Thread loop for multithreading
void threadLoop(SimulationTools * tools);

bool tryPop(SimulationTools * tools, DataPoint & point);

// Initial conditions to allow step-by-step solving
void initialize(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools);

// Performs a time step
void step(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools);

// Applies an operation (pi-pulse, ...)
void applyOperation(const Operation & operation, Solution & solution, SimulationTools * tools);

// Compute the phase of the spin
void computePhase(const SimulationTools & tools);

int main(int argc, char * argv[])
{
    auto startTime = chrono::high_resolution_clock::now();

    // Necessary to use QScriptEngine and QSqlDatabase
    QCoreApplication app(argc, argv); Q_UNUSED(app);

    // Wrapper for various simulation settings
    SimulationTools tools;

    // Reading the settings of the simulation, exit if the settings are invalid or if the file is not found
    if(argc < 2 || !processFile(QString(argv[1]), tools)) {
        qDebug() << "No/invalid script provided.";
        return -1;
    }

    qDebug() << tools.queue().size() << "points to evaluate. Preparing the threads...";

    // Preparing & starting the threads
    unsigned int threadCount = thread::hardware_concurrency()-1;

    thread * threadPool[threadCount == 0 ? 1 : threadCount];

    for(unsigned int i = 0; i < threadCount; ++i) {
        if(i < tools.queue().size()-1)
            threadPool[i] = new thread(threadLoop, &tools);
        else
            threadPool[i] = NULL;
    }

    qDebug() << "All threads started" << flush;

    threadLoop(&tools);

    // Waiting for all thread to complete
    for(unsigned int i = 0; i < threadCount; i++) {
        if(threadPool[i]) {
            threadPool[i]->join();
            delete threadPool[i];
        }
    }

    // Compute the phase of the spin
    computePhase(tools);

    auto endTime = chrono::high_resolution_clock::now();

    qDebug() << "Simulation succesfully completed. Elapsed time (msec):" << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();

    return 0;
}

void threadLoop(SimulationTools * tools)
{
    DataPoint dataPoint;

    // Memory pre-allocation for the solution
    Solution solution(*tools);

    // Database
    QSqlDatabase db = QSqlDatabase::cloneDatabase(tools->database(), "isre-result-" + QString::number(hash<thread::id>()(this_thread::get_id())));
    if(!db.open()) {
        qDebug() << "Failed to open the database";
    }

    QSqlQuery query(db);
    // Energy bin rowids
    int * energy_bin_id  = new int[tools->energyBinCount()];

    // Main loop
    while(true) {
        // Obtaining the data to process
        if(!tryPop(tools, dataPoint)) {
            // Cleaning
            delete [] energy_bin_id;

            // Exiting
            return;
        }

        // Checking numerical value
        double timeStep = dataPoint.totalTime / dataPoint.timeStepCount;

        if(timeStep < 1e-6)
            qDebug() << "Warning: the time step is close to or below numerical precision 1e-7 (" << timeStep << ")";
        if(abs(timeStep*dataPoint.larmor) > 1e-1)
            qDebug() << "Warning: time step * larmor frequency is close to or higher than 1 (" << timeStep*dataPoint.larmor << ")";
        if(abs(timeStep*dataPoint.larmorInhomogeneity*tools->energy(tools->energyCount()-1)) > 1e-1)
            qDebug() << "Warning: time step * larmo inh. * max. energy  is close to or higher than 1 (" << timeStep*dataPoint.larmorInhomogeneity*tools->energy(tools->energyCount()-1) << ")";
        if(dataPoint.larmorInhomogeneity != 0 && abs(timeStep*dataPoint.larmorInhomogeneity*tools->energy(1)) < 1e-14)
            qDebug() << "Warning: larmor inh. * min. energy is close to or below numerical precision 1e-15 (" << timeStep*dataPoint.larmorInhomogeneity*tools->energy(1) << ")";
        if(abs(timeStep*dataPoint.exchange) > 1e-1)
            qDebug() << "Warning: timestep * exchange frequency is to close to or higher than 1 (" << timeStep*dataPoint.exchange << ")";

        // Numerical solution

        initialize(dataPoint, solution, tools);

        const Operation * next = NULL;

        while((next = tools->nextOperation(double(solution.timeStep) / double(dataPoint.timeStepCount) * dataPoint.totalTime, dataPoint.totalTime))) {
            int nextTimestep = (next->isFixedTime ? next->time / dataPoint.totalTime : next->time ) * dataPoint.timeStepCount;
            while(solution.timeStep < nextTimestep)
                step(dataPoint, solution, tools);

            applyOperation(*next, solution, tools);
        }

        while(solution.timeStep < dataPoint.timeStepCount)
            step(dataPoint, solution, tools);

        // When using moving operations we are interested in only one point, hence we didn't save anything in solution.*Spin_*_save;
        if(tools->hasMovingOperation()) {
            solution.totalSpin_x_save[0] = solution.totalSpin_x;
            solution.totalSpin_y_save[0] = solution.totalSpin_y;
            solution.totalSpin_z_save[0] = solution.totalSpin_z;

            for(int i = 0; i < tools->energyBinCount(); ++i) {
                // TODO: change tools->energyBin(i).first so that is it corresponds to i_begin (if possible)
                int i_begin = 0;
                while(i_begin < tools->energyCount() && tools->energy(i_begin) < tools->energyBin(i).first)
                    ++i_begin;
                int i_end = i_begin;
                while(i_end < tools->energyCount() && tools->energy(i_end) <= tools->energyBin(i).second)
                    ++i_end;

                solution.partialSpin_x_save[i][0] = accumulate(solution.spin_x.begin() + i_begin, solution.spin_x.begin() + i_end, 0.0) / double(i_end-i_begin);
                solution.partialSpin_y_save[i][0] = accumulate(solution.spin_y.begin() + i_begin, solution.spin_y.begin() + i_end, 0.0) / double(i_end-i_begin);
                solution.partialSpin_z_save[i][0] = accumulate(solution.spin_z.begin() + i_begin, solution.spin_z.begin() + i_end, 0.0) / double(i_end-i_begin);
            }
        }

        // Saving data
        tools->databaseMutex()->lock();

        // Starting the SQL transaction
        query.clear();
        query.exec("BEGIN");

        // In case of spin echo we compute the solution point by point, while whithout it we compute all points at once
        int t_min = 0;
        int t_max = tools->hasMovingOperation() ? 0 : tools->timePointsCount();

        for(int t = t_min; t <= t_max; ++t) {
            // Time we are svaing the values for
            double time = tools->hasMovingOperation() ? dataPoint.totalTime : t * dataPoint.totalTime / tools->timePointsCount();

            double contrast = sqrt(pow(solution.totalSpin_x_save[t],2) + pow(solution.totalSpin_y_save[t],2) + pow(solution.totalSpin_z_save[t],2));

            // Total Spin
            query.prepare("INSERT INTO result (simulation_id, time, s_x, s_y, s_z, contrast, lost_1, lost_2) VALUES (:s, :t, :x, :y, :z, :c, :l1, :l2);");
            query.bindValue(":s", dataPoint.simulationId);
            query.bindValue(":t", time);
            query.bindValue(":x", solution.totalSpin_x_save[t]);
            query.bindValue(":y", solution.totalSpin_y_save[t]);
            query.bindValue(":z", solution.totalSpin_z_save[t]);
            query.bindValue(":c", contrast);
            query.bindValue(":l1", solution.spectatorPopulationF1_save[t]);
            query.bindValue(":l2", solution.spectatorPopulationF2_save[t]);

            query.exec();

            // Partial spin
            for(int i = 0; i < tools->energyBinCount(); ++i) {
                query.prepare("INSERT INTO energy_bin_result (energy_bin_id, time, s_x, s_y, s_z) VALUES (:i, :t, :x, :y, :z);");
                query.bindValue(":i", dataPoint.energyBinId.at(i));
                query.bindValue(":t", time);
                query.bindValue(":x", solution.partialSpin_x_save[i][t]);
                query.bindValue(":y", solution.partialSpin_y_save[i][t]);
                query.bindValue(":z", solution.partialSpin_z_save[i][t]);

                query.exec();
            }
        }

        // Commit all the SQL data
        query.exec("COMMIT");
        query.clear();

        tools->databaseMutex()->unlock();
    }
}

bool tryPop(SimulationTools * tools, DataPoint & point){
    lock_guard<mutex> lock(*tools->queueMutex()); Q_UNUSED(lock);

    if(tools->queue().empty()) {
        return false;
    }

    // Pops a point from the queue
    point = tools->queue().front();
    tools->queue().pop();

    return true;
}


// Initialize the first and second time step of the solution, as the computing the solution at time t requires to know the solution at time t-1 and t-2
void initialize(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools)
{
    solution.timeStep = 2;

    // Clearing all previous values
    memset(solution.totalSpin_x_save, 0, (tools->timePointsCount()+1)*sizeof(double));
    memset(solution.totalSpin_y_save, 0, (tools->timePointsCount()+1)*sizeof(double));
    memset(solution.totalSpin_z_save, 0, (tools->timePointsCount()+1)*sizeof(double));
    // no total_spin_n_save
    memset(solution.spectatorPopulationF1_save, 0, (tools->timePointsCount()+1)*sizeof(double));
    memset(solution.spectatorPopulationF2_save, 0, (tools->timePointsCount()+1)*sizeof(double));

    for(int i = 0; i < tools->energyBinCount(); ++i) {
        memset(solution.partialSpin_x_save[i], 0, (tools->timePointsCount()+1)*sizeof(double));
        memset(solution.partialSpin_y_save[i], 0, (tools->timePointsCount()+1)*sizeof(double));
        memset(solution.partialSpin_z_save[i], 0, (tools->timePointsCount()+1)*sizeof(double));
        // No partial_spin_n_save
    }

    // Values at the initial time
    solution.totalSpin_x_old = 1.0;
    solution.totalSpin_y_old = 0;
    solution.totalSpin_z_old = 0;
    solution.totalSpin_n_old = 1.0;

    solution.totalSpin_x_save[0] = 1.0;
    solution.totalSpin_y_save[0] = 0;
    solution.totalSpin_z_save[0] = 0;

    solution.spectatorPopulationF1_save[0] = 0;
    solution.spectatorPopulationF2_save[0] = 0;

    for(int i = 0; i < tools->energyBinCount(); ++i) {
        solution.partialSpin_x_save[i][0] = 1.0;
        solution.partialSpin_y_save[i][0] = 0;
        solution.partialSpin_z_save[i][0] = 0;
    }

    // Duration of the time steps in seconds
    double dt = double(dataPoint.totalTime) / double(dataPoint.timeStepCount);

    // Values at the first time step
    for(int i = 0; i < tools->energyCount(); ++i) {
        double b = dataPoint.larmor
                + dataPoint.larmorInhomogeneity * tools->energy(i)
                + tools->densityAtTimestep(0) * dataPoint.densityInhomogeneity * tools->densityAtEnergy(i);

        double g = tools->densityAtTimestep(0) * 0.25 * (dataPoint.lossFromF1 - dataPoint.lossFromF2);

        double norm = sqrt(1 + (dt*b)*(dt*b) + (dt*g)*(dt*g));

        solution.spin_x[i] = 1.0 / norm;
        solution.spin_y[i] = - dt * b / norm;
        solution.spin_z[i] = - dt * g / norm;
        solution.spin_n[i] = 1.0 - dt * tools->densityAtTimestep(0) * 0.25 * (dataPoint.lossFromF1 + dataPoint.lossFromF2);

        solution.totalSpin_x += solution.spin_x[i];
        solution.totalSpin_y += solution.spin_y[i];
        solution.totalSpin_z += solution.spin_z[i];
        solution.totalSpin_n += solution.spin_n[i];
    }
    solution.spectatorPopulationF1 = dt * dataPoint.lossFromF1;
    solution.spectatorPopulationF2 = dt * dataPoint.lossFromF2;

    solution.totalSpin_x /= tools->energyCount();
    solution.totalSpin_y /= tools->energyCount();
    solution.totalSpin_z /= tools->energyCount();
    solution.totalSpin_n /= tools->energyCount();
}

// Given the solution of the equation up to the time step solution.timestep - 1, computes the solution at the sime step solution.timestep and increases timestep by one
void step(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools)
{
    Matrix4d Re;
    Vector4d se;

    Matrix4d G;
    Vector4d s;

    Vector4d result;

    // The equation is result = (2-Re dt)^(-1) (2 + Re dt) se + (2-Re dt)^(-1) G s

    /* To avoid temporaries */

    Matrix4d R;  // R contains the energy-independent part of Re;
    Matrix4d RT; // = 2 + Re * dt
    Matrix4d RI; // = (2 - Re * dt).inverse();

    /* Total Spin */

    s <<    0.5 * (3.0*solution.totalSpin_x - solution.totalSpin_x_old),    // s(0) = Sx
            0.5 * (3.0*solution.totalSpin_y - solution.totalSpin_y_old),    // s(1) = Sy
            0.5 * (3.0*solution.totalSpin_z - solution.totalSpin_z_old),    // s(2) = Sz
            0.5 * (3.0*solution.totalSpin_n - solution.totalSpin_n_old);    // s(3) = n

    /* Aliases */

    // Timestep
    int t = solution.timeStep;
    // Duration of a time step in seconds
    double dt = dataPoint.totalTime / dataPoint.timeStepCount;
    double df = tools->densityAtTimestep(t);

    // Total atom losses to spectator states
    double dfls = df * 0.25 * (dataPoint.lossFromF1 + dataPoint.lossFromF2);
    // Difference of atom losses to spectator states
    double dfld = df * 0.25 * (dataPoint.lossFromF1 - dataPoint.lossFromF2);

    // TODO: Est-ce que la relaxation (colisions qui changent l'energie) est proportionelle au nombre d'atomes
    // ou seulement au nombre d'atomes dans F=1,2 ? (Réponse temporaire : nombre d'atomes total)
    double dfd = df * dataPoint.damping;
    // TODO: idem pour l'echange (Réponse temporaire : dans F=1,2 seulement)
    double dfe = df * dataPoint.exchange;

    /* Matrix & vector initialization: energy-independent part */

    // Terms propotional to the total density
    R <<    -dfd,       -dfe*s(2),   dfe*s(1),                      0,
             dfe*s(2),  -dfd,       -dfe*s(0),                      0,
            -dfe*s(1),   dfe*s(0),  -(dfls*s(3) + dfld*s(2) + dfd), -(dfld*s(3) + dfls*s(2)),
            0,          0,          -(dfld*s(3) + dfls*s(2)),       -(dfls*s(3) + dfld*s(2));

    // Larmor precession around +z
    R(0,1) -= dataPoint.larmor;
    R(1,0) += dataPoint.larmor;

    G << dfd,   0,      0,      0,
         0,     dfd,    0,      0,
         0,     0,      dfd,    0,
         0,     0,      0,      0;

    // Solving the equation of motion for the individual spins
    for(int i = 0; i < tools->energyCount(); ++i) {

        /* Matrix & vector initialization: energy dependent part */
        Re = R;

        // Change of larmor frequency dependant on the energy
        // For a harmonic trap, the larmor inhomogeneity is proprotional to ?, while the density proportional to ?
        // TODO: inhomogeneité de champ moyen : proportionel à n_total ou a n(F=1,2) ? (Réponse temporaire: n(F=1,2))
        double add = dataPoint.larmorInhomogeneity * tools->energy(i) + df * dataPoint.densityInhomogeneity * tools->densityAtEnergy(i);
        Re(0,1) -= add;
        Re(1,0) += add;

        se(0) = solution.spin_x[i];
        se(1) = solution.spin_y[i];
        se(2) = solution.spin_z[i];
        se(3) = solution.spin_n[i];

        /* Computation */

        // We can't inverse in place, so we use RT as a temporary to compute RI
        RT = 2.0 * Matrix4d::Identity() - dt * Re;
        RI = RT.inverse();
        RT = 2.0 * Matrix4d::Identity() + dt * Re;

        result = RI * (RT * se + dt * G * s);

        solution.spin_x[i] = result(0);
        solution.spin_y[i] = result(1);
        solution.spin_z[i] = result(2);
        solution.spin_n[i] = result(3);
    }

    // Spectator population update
    solution.spectatorPopulationF1 += dt * df * dataPoint.lossFromF1 * 0.5 * (s(3) + s(2));
    solution.spectatorPopulationF2 += dt * df * dataPoint.lossFromF2 * 0.5 * (s(3) - s(2));

    // Value of the total spin at the previous time step
    solution.totalSpin_x_old = solution.totalSpin_x;
    solution.totalSpin_y_old = solution.totalSpin_y;
    solution.totalSpin_z_old = solution.totalSpin_z;
    solution.totalSpin_n_old = solution.totalSpin_n;

    // Value of the total spin at the current time step
    solution.totalSpin_x = accumulate(solution.spin_x.begin(), solution.spin_x.end(), 0.0) / tools->energyCount();
    solution.totalSpin_y = accumulate(solution.spin_y.begin(), solution.spin_y.end(), 0.0) / tools->energyCount();
    solution.totalSpin_z = accumulate(solution.spin_z.begin(), solution.spin_z.end(), 0.0) / tools->energyCount();
    solution.totalSpin_n = accumulate(solution.spin_n.begin(), solution.spin_n.end(), 0.0) / tools->energyCount();

    // Saving the data - not needed for 'moving' spin echo as we are only interested in the final point

    if(!tools->hasMovingOperation() && t % ((int)std::floor(dataPoint.timeStepCount/tools->timePointsCount())) == 0) {
        t = t / std::floor(dataPoint.timeStepCount/tools->timePointsCount());

        solution.totalSpin_x_save[t] = solution.totalSpin_x;
        solution.totalSpin_y_save[t] = solution.totalSpin_y;
        solution.totalSpin_z_save[t] = solution.totalSpin_z;

        solution.spectatorPopulationF2_save[t] = solution.spectatorPopulationF1;
        solution.spectatorPopulationF2_save[t] = solution.spectatorPopulationF2;

        for(int i = 0; i < tools->energyBinCount(); ++i) {
            // TODO: change tools->energyBin(i).first so that is it corresponds to i_begin (if possible)
            int i_begin = 0;
            while(i_begin < tools->energyCount() && tools->energy(i_begin) < tools->energyBin(i).first)
                ++i_begin;
            int i_end = i_begin;
            while(i_end < tools->energyCount() && tools->energy(i_end) <= tools->energyBin(i).second)
                ++i_end;

            solution.partialSpin_x_save[i][t] = accumulate(solution.spin_x.begin() + i_begin, solution.spin_x.begin() + i_end, 0.0) / double(i_end-i_begin);
            solution.partialSpin_y_save[i][t] = accumulate(solution.spin_y.begin() + i_begin, solution.spin_y.begin() + i_end, 0.0) / double(i_end-i_begin);
            solution.partialSpin_z_save[i][t] = accumulate(solution.spin_z.begin() + i_begin, solution.spin_z.begin() + i_end, 0.0) / double(i_end-i_begin);
        }
    }

    solution.timeStep++;
}

// Applies an operation at the latest timestep (at solution.timestep-1)
// solution.timestep is not increased
void applyOperation(const Operation & operation, Solution & solution, SimulationTools * tools)
{
    for(int i = 0; i < tools->energyCount(); ++i) {
        Vector3d s(solution.spin_x[i], solution.spin_y[i], solution.spin_z[i]);

        s = operation.operationMatrix * s;

        solution.spin_x[i] = s(0);
        solution.spin_y[i] = s(1);
        solution.spin_z[i] = s(2);
    }

    Vector3d s_total(solution.totalSpin_x, solution.totalSpin_y, solution.totalSpin_z);

    s_total = operation.operationMatrix * s_total;

    solution.totalSpin_x = s_total(0);
    solution.totalSpin_y = s_total(1);
    solution.totalSpin_z = s_total(2);

    Vector3d s_total_old(solution.totalSpin_x_old, solution.totalSpin_y_old, solution.totalSpin_z_old);

    s_total_old = operation.operationMatrix * s_total;

    solution.totalSpin_x_old = s_total_old(0);
    solution.totalSpin_y_old = s_total_old(1);
    solution.totalSpin_z_old = s_total_old(2);
}

void computePhase(const SimulationTools & tools)
{
    QSqlQuery query(tools.database());
    query.exec("BEGIN");

    QSqlQuery queryId(tools.database());

    queryId.prepare("SELECT id FROM energy_bin WHERE simulation_id IN (SELECT id FROM simulation WHERE name = :n)");
    queryId.bindValue(":n", tools.simulationName());
    queryId.exec();

    while(queryId.next())
    {
        int energyBinId = queryId.record().value(0).toInt();

        QSqlQuery queryResult(tools.database());

        queryResult.prepare("SELECT * FROM energy_bin_result WHERE energy_bin_id = :i ORDER BY time ASC");
        queryResult.bindValue(":i", energyBinId);
        queryResult.exec();

        queryResult.next();

        auto norm_xy = [&]() -> double {
            return sqrt(pow(queryResult.record().value("s_x").toDouble(), 2) + pow(queryResult.record().value("s_y").toDouble(), 2));
        };

        double sx   = queryResult.record().value("s_x").toDouble();
        double sy   = queryResult.record().value("s_y").toDouble();
        double c    = norm_xy();
        double p    = 0;

        QSqlQuery queryUpdate(tools.database());

        // For each selected entry of the table
        while(true) {
            double sxo = sx;
            double syo = sy;
            double co = c;

            if(!queryResult.next()) {
                break;
            }

            sx  = queryResult.record().value("s_x").toDouble();
            sy  = queryResult.record().value("s_y").toDouble();
            c   = norm_xy();

            p += asin((sxo*sy - sx*syo)/(c*co));

            queryUpdate.prepare("UPDATE energy_bin_result SET phase = :p WHERE energy_bin_id = :i and time = :t");
            queryUpdate.bindValue(":p", p);
            queryUpdate.bindValue(":i", energyBinId);
            queryUpdate.bindValue(":t", queryResult.record().value("time").toDouble());
            queryUpdate.exec();
        }
    }

    queryId.prepare("SELECT id FROM simulation WHERE name = :n");
    queryId.bindValue(":n", tools.simulationName());
    queryId.exec();

    while(queryId.next()) {
        int simulationId = queryId.record().value(0).toInt();

        QSqlQuery queryResult(tools.database());

        queryResult.prepare("SELECT * FROM result WHERE simulation_id = :i ORDER BY time ASC");
        queryResult.bindValue(":i", simulationId);
        queryResult.exec();

        queryResult.next();

        auto norm_xy = [&]() -> double {
            return sqrt(pow(queryResult.record().value("s_x").toDouble(), 2) + pow(queryResult.record().value("s_y").toDouble(), 2));
        };

        QSqlQuery queryUpdate(tools.database());

        queryUpdate.prepare("UPDATE result SET phase = 0 WHERE simulation_id = :i and time = 0");
        queryUpdate.bindValue(":i", simulationId);
        queryUpdate.exec();

        double sx   = queryResult.record().value("s_x").toDouble();         // Spin s_x
        double sy   = queryResult.record().value("s_y").toDouble();         // Spin s_y
        double c    = norm_xy();                                            // Norm
        double p    = 0;                                                    // Phase

        // For each selected entry
        while(true) {
            double sxo = sx;
            double syo = sy;
            double co = c;

            if(!queryResult.next()) {
                break;
            }

            sx  = queryResult.record().value("s_x").toDouble();
            sy  = queryResult.record().value("s_y").toDouble();
            c   = norm_xy();

            p  += asin((sxo*sy - sx*syo)/(c*co));

            queryUpdate.prepare("UPDATE result SET phase = :p WHERE simulation_id = :i and time = :t");
            queryUpdate.bindValue(":p", p);
            queryUpdate.bindValue(":i", simulationId);
            queryUpdate.bindValue(":t", queryResult.record().value("time").toDouble());
            queryUpdate.exec();
        }
    }

    query.exec("COMMIT");
}
