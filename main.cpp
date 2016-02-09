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

const double PI = M_PI;

using namespace std;

// Thread loop for multithreading
void threadLoop(SimulationTools * tools);

// Initial conditions to allow step-by-step solving
void initialize(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools);

// Performs a time step
void step(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools);

// Realizes a spin echo
void spinEcho(Solution & solution, SimulationTools * tools);

// Compute the phase of the spin
void computePhase(const SimulationTools & tools);

int main(int argc, char * argv[])
{
    auto startTime = chrono::high_resolution_clock::now();

    // Necessary to use QScriptEngine and QSqlDatabase
    QCoreApplication app(argc, argv);

    // Wrapper for various simulation settings
    SimulationTools tools;

    // Reading the settings of the simulation, exit if the settings are invalid or if the file is not found
    if(argc < 2 || !processFile(QString(argv[1]), tools))
    {
        qDebug() << "No script provided.";
        return -1;
    }

    // Preparing & starting the threads
    unsigned int threadCount = thread::hardware_concurrency()-1;

    qDebug() << tools.queue().size() << "points to evaluate. Preparing the threads...";

    thread * threadPool[threadCount == 0 ? 1 : threadCount];

    for(unsigned int i = 0; i < threadCount; i++)
    {
        if(i < tools.queue().size()-1)
            threadPool[i] = new thread(threadLoop, & tools);
        else
            threadPool[i] = NULL;
    }

    qDebug() << "All threads started" << flush;

    threadLoop(&tools);

    // Waiting for all thread to complete
    for(unsigned int i = 0; i < threadCount; i++)
    {
        if(threadPool[i])
        {
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
    DataPoint   dataPoint;

    // Memory pre-allocation for the solution
    Solution    solution;

    solution.spin_x             = new real[tools->energyCount()];
    solution.spin_y             = new real[tools->energyCount()];
    solution.spin_z             = new real[tools->energyCount()];

    solution.totalSpin_x_save   = new real[tools->timePointsCount()+1];
    solution.totalSpin_y_save   = new real[tools->timePointsCount()+1];
    solution.totalSpin_z_save   = new real[tools->timePointsCount()+1];

    solution.partialSpin_x_save = new real*[tools->energyBinCount()+1];
    solution.partialSpin_y_save = new real*[tools->energyBinCount()+1];
    solution.partialSpin_z_save = new real*[tools->energyBinCount()+1];

    for(int i = 0; i < tools->energyBinCount(); i++)
    {
        solution.partialSpin_x_save[i] = new real[tools->timePointsCount()+1];
        solution.partialSpin_y_save[i] = new real[tools->timePointsCount()+1];
        solution.partialSpin_z_save[i] = new real[tools->timePointsCount()+1];
    }

    // Database
    QSqlDatabase db = QSqlDatabase::cloneDatabase(tools->database(), "isre-result-" + QString::number(hash<thread::id>()(this_thread::get_id())));
    if(!db.open())
    {
        qDebug() << "Failed to open the database";
    }

    QSqlQuery query(db);
    // Energy bin rowids
    int * energy_bin_id  = new int[tools->energyBinCount()];

    // Main loop
    while(true)
    {
        // Obtaining the data to process
        tools->queueMutex()->lock();

        if(tools->queue().empty())
        {
            // All data has been processed.
            tools->queueMutex()->unlock();

            // Cleaning
            delete [] solution.spin_x;
            delete [] solution.spin_y;
            delete [] solution.spin_z;

            delete [] solution.totalSpin_x_save;
            delete [] solution.totalSpin_y_save;
            delete [] solution.totalSpin_z_save;

            for(int i = 0; i < tools->energyBinCount(); i++)
            {
                delete [] solution.partialSpin_x_save[i];
                delete [] solution.partialSpin_y_save[i];
                delete [] solution.partialSpin_z_save[i];
            }

            delete [] solution.partialSpin_x_save;
            delete [] solution.partialSpin_y_save;
            delete [] solution.partialSpin_z_save;

            delete [] energy_bin_id;

            // Exiting
            return;
        }
        else
        {
            // Pops a point from the queue
            dataPoint = tools->queue().front();
            tools->queue().pop();

        }

        tools->queueMutex()->unlock();

        // Checking numerical value
        real timeStep = real(dataPoint.totalTime) / real(dataPoint.timeStepCount);

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

        for(int i = 0; i < tools->spinEchoCount(); i++)
        {
            // Gets the next spin echo time step
            int spinEchoTimeStep;

            // If tools->spinEcho(), spin echo times are reltaive to the simulation time
            if(tools->spinEcho())
                spinEchoTimeStep = tools->spinEchoTimeStep(i) * dataPoint.timeStepCount / tools->timeStepsCount();
            // if not, then the spin echo times are absolute
            else
                spinEchoTimeStep = tools->spinEchoTimeStep(i);

            // Free evolution until the next spin echo
            for(int t = solution.timeStep; t < spinEchoTimeStep && t <= dataPoint.timeStepCount; t++)
                step(dataPoint, solution, tools);

            // Performs the spin echo
            spinEcho(solution, tools);
        }

        // Free evolution until the end

        for(int t = solution.timeStep; t <= dataPoint.timeStepCount; t++)
            step(dataPoint, solution, tools);


        // When using spin echo we are interested in only one point, hence we didn't save anything in solution.*Spin_*_save;
        if(tools->spinEcho())
        {
            solution.totalSpin_x_save[0] = solution.totalSpin_x;
            solution.totalSpin_y_save[0] = solution.totalSpin_y;
            solution.totalSpin_z_save[0] = solution.totalSpin_z;

            for(int i = 0; i < tools->energyBinCount(); ++i)
            {
                solution.partialSpin_x_save[i][0] = 0;
                solution.partialSpin_y_save[i][0] = 0;
                solution.partialSpin_z_save[i][0] = 0;

                int j = 0;

                while(tools->energy(j) < tools->energyBin(i).first)
                    ++j;
                while(tools->energy(j) < tools->energyBin(i).second)
                {
                    solution.partialSpin_x_save[j][0] += solution.spin_x[j];
                    solution.partialSpin_y_save[j][0] += solution.spin_y[j];
                    solution.partialSpin_z_save[j][0] += solution.spin_z[j];
                    ++j;
                }
            }
        }

        // Saving data
        tools->databaseMutex()->lock();

        // Starting the SQL transaction
        query.clear();
        query.exec("BEGIN");

        // In case of spin echo we compute the solution point by point, while whithout it we compute all points at once
        int t_min = 0;
        int t_max = tools->spinEcho() ? 0 : tools->timePointsCount();

        for(int t = t_min; t <= t_max; ++t)
        {
            // Time we are svaing the values for
            double time = tools->spinEcho() ? dataPoint.totalTime : t * dataPoint.totalTime / tools->timePointsCount();

            // Total Spin
            double norm = sqrt(pow(solution.totalSpin_x_save[t], 2.0) + pow(solution.totalSpin_y_save[t], 2.0) + pow(solution.totalSpin_z_save[t], 2.0));

            query.prepare("INSERT INTO result (simulation_id, time, s_x, s_y, s_z, contrast) VALUES (:s, :t, :x, :y, :z, :c);");
            query.bindValue(":s", dataPoint.simulationId);
            query.bindValue(":t", time);
            query.bindValue(":x", solution.totalSpin_x_save[t]);
            query.bindValue(":y", solution.totalSpin_y_save[t]);
            query.bindValue(":z", solution.totalSpin_z_save[t]);
            query.bindValue(":c", norm);

            query.exec();

            // Partial spin
            for(int i = 0; i < tools->energyBinCount(); ++i)
            {
                norm = sqrt(pow(solution.partialSpin_x_save[i][t], 2.0) + pow(solution.partialSpin_y_save[i][t], 2.0) + pow(solution.partialSpin_z_save[i][t], 2.0));

                query.prepare("INSERT INTO energy_bin_result (energy_bin_id, time, s_x, s_y, s_z, contrast) VALUES (:i, :t, :x, :y, :z, :c);");
                query.bindValue(":i", dataPoint.energyBinId.at(i));
                query.bindValue(":t", time);
                query.bindValue(":x", solution.partialSpin_x_save[i][t]);
                query.bindValue(":y", solution.partialSpin_y_save[i][t]);
                query.bindValue(":z", solution.partialSpin_z_save[i][t]);
                query.bindValue(":c", norm);

                query.exec();
            }
        }

        // Commit all the SQL data
        query.exec("COMMIT");
        query.clear();

        tools->databaseMutex()->unlock();
    }
}

// Initialize the first and second time step of the solution, as the computing the solution at time t requires to know the solution at time t-1 and t-2
void initialize(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools)
{
    solution.timeStep = 2;

    // Duration of the time steps in seconds
    real timeStep = real(dataPoint.totalTime) / real(dataPoint.timeStepCount);

    // Clearing all previous values
    memset(solution.totalSpin_x_save, 0, (tools->timePointsCount()+1)*sizeof(real));
    memset(solution.totalSpin_y_save, 0, (tools->timePointsCount()+1)*sizeof(real));
    memset(solution.totalSpin_z_save, 0, (tools->timePointsCount()+1)*sizeof(real));

    for(int i = 0; i < tools->energyBinCount(); i++)
    {
        memset(solution.partialSpin_x_save[i], 0, (tools->timePointsCount()+1)*sizeof(real));
        memset(solution.partialSpin_x_save[i], 0, (tools->timePointsCount()+1)*sizeof(real));
        memset(solution.partialSpin_x_save[i], 0, (tools->timePointsCount()+1)*sizeof(real));
    }

    // Values at the initial time
    solution.totalSpin_x_old = real(1.0);
    solution.totalSpin_y_old = 0;
    solution.totalSpin_z_old = 0;

    solution.totalSpin_x_save[0] = real(1.0);
    solution.totalSpin_y_save[0] = 0;
    solution.totalSpin_z_save[0] = 0;

    for(int i = 0; i < tools->energyBinCount(); i++)
    {
        solution.partialSpin_x_save[i][0] = 1.0;
        solution.partialSpin_y_save[i][0] = 0;
        solution.partialSpin_z_save[i][0] = 0;
    }

    // Values at the first time step
    for(int i = 0; i < tools->energyCount(); i++)
    {
        real norm = sqrt(1.0 + pow(timeStep *(dataPoint.larmor + dataPoint.larmorInhomogeneity * tools->energy(i)), 2.0));

        solution.spin_x[i] = real(1.0)/norm;
        solution.spin_y[i] = -timeStep * (dataPoint.larmor + dataPoint.larmorInhomogeneity * tools->energy(i))/norm;
        solution.spin_z[i] = 0;

        solution.totalSpin_x += solution.spin_x[i];
        solution.totalSpin_y += solution.spin_y[i];
        solution.totalSpin_z += solution.spin_z[i];
    }

    solution.totalSpin_x /= tools->energyCount();
    solution.totalSpin_y /= tools->energyCount();
    solution.totalSpin_z /= tools->energyCount();
}

// Given the solution of the equation up to the time step solution.timestep - 1, computes the solution at the sime step solution.timestep and increases timestep by one
void step(const DataPoint & dataPoint, Solution & solution, SimulationTools * tools)
{
    int t = solution.timeStep;

    // Duration of the time steps in seconds
    real dt = real(dataPoint.totalTime) / real(dataPoint.timeStepCount);
    real dt2 = dt*dt;

    // Densiy factor (density at timestep t divided by density at t = 0)
    real df = tools->densityFactor(t);

    if(dataPoint.damping != 0)
    {
        // Damping appears only as dt * damping or as damping * damping
        real dtg =  df * dt * dataPoint.damping;
        real g2 = df * df * dataPoint.damping * dataPoint.damping;

        // x and y components of the total spin do not depend on the energy
        real tsx = df * dataPoint.exchange * 0.5 * (3.0*solution.totalSpin_x - solution.totalSpin_x_old);
        real tsx2 = tsx*tsx;
        real tsy = df * dataPoint.exchange * 0.5 * (3.0*solution.totalSpin_y - solution.totalSpin_y_old);
        real tsy2 = tsy*tsy;
        real tmp_tsz = df * dataPoint.exchange * 0.5 * (3.0*solution.totalSpin_z - solution.totalSpin_z_old) + dataPoint.larmor;

        real tsxdtg = dtg * (3.0*solution.totalSpin_x - solution.totalSpin_x_old);
        real tsydtg = dtg * (3.0*solution.totalSpin_y - solution.totalSpin_y_old);
        real tszdtg = dtg * (3.0*solution.totalSpin_z - solution.totalSpin_z_old);

        // Coefficient that do not depend on tsz
        real Nxx = 4.0 + 4.0*dtg + dt2*(tsx2 + g2);
        real Nyy = 4.0 + 4.0*dtg + dt2*(tsy2 + g2);

        // Equation of motion for the individual spins
        for(int i = 0; i < tools->energyCount(); i++)
        {
            // For a harmonic trap, the larmor inhomogeneity is proprotional to the energy, while the density proportional to exp(-energy)
            real tsz = tmp_tsz + dataPoint.larmorInhomogeneity*tools->energy(i) + df * dataPoint.densityInhomogeneity*exp(-tools->energy(i));
            real tsz2 = tsz*tsz;

            real dt2s2g2 = dt2*(tsx2 + tsy2 + tsz2 + g2);
            real det = (2.0 + dtg)*(4.0 + 4.0*dtg + dt2s2g2);

            real Mxx = 8.0 + (4.0-dt2s2g2)*dtg + dt2*(tsx2-tsy2-tsz2-g2);
            real Mxy = 4.0*(dt2*tsx*tsy - dt*tsz*(2.0 + dtg));
            real Mxz = 4.0*(dt2*tsx*tsz + dt*tsy*(2.0 + dtg));

            real Myx = 4.0*(dt2*tsx*tsy + dt*tsz*(2.0 + dtg));
            real Myy = 8.0 + (4.0-dt2s2g2)*dtg + dt2*(tsy2-tsx2-tsz2-g2);
            real Myz = 4.0*(dt2*tsy*tsz - dt*tsx*(2.0 + dtg));

            real Mzx = 4.0*(dt2*tsx*tsz - dt*tsy*(2.0 + dtg));
            real Mzy = 4.0*(dt2*tsy*tsz + dt*tsx*(2.0 + dtg));
            real Mzz = 8.0 + (4.0-dt2s2g2)*dtg + dt2*(tsz2-tsx2-tsy2-g2);

            real Nxy = dt2*tsx*tsy - dt*tsz*(2.0 + dtg);
            real Nxz = dt2*tsx*tsz + dt*tsy*(2.0 + dtg);

            real Nyx = dt2*tsx*tsy + dt*tsz*(2.0 + dtg);
            real Nyz = dt2*tsy*tsz - dt*tsx*(2.0 + dtg);

            real Nzx = dt2*tsx*tsz - dt*tsy*(2.0 + dtg);
            real Nzy = dt2*tsy*tsz + dt*tsx*(2.0 + dtg);
            real Nzz = 4.0 + 4.0*dtg + dt2*(tsz2 + g2);

            real sx = solution.spin_x[i];
            real sy = solution.spin_y[i];
            real sz = solution.spin_z[i];

            solution.spin_x[i] = ((Mxx*sx + Mxy*sy + Mxz*sz) + (Nxx*tsxdtg + Nxy*tsydtg + Nxz*tszdtg))/det;
            solution.spin_y[i] = ((Myx*sx + Myy*sy + Myz*sz) + (Nyx*tsxdtg + Nyy*tsydtg + Nyz*tszdtg))/det;
            solution.spin_z[i] = ((Mzx*sx + Mzy*sy + Mzz*sz) + (Nzx*tsxdtg + Nzy*tsydtg + Nzz*tszdtg))/det;
        }
    }
    else
    {
        // Alias for convenience

        // x and y components of the total spin do not depend on the energy
        real tsx = df * dataPoint.exchange * 0.5 * (3.0*solution.totalSpin_x - solution.totalSpin_x_old);
        real tsx2 = tsx*tsx;
        real tsy = df * dataPoint.exchange * 0.5 * (3.0*solution.totalSpin_y - solution.totalSpin_y_old);
        real tsy2 = tsy*tsy;
        real tmp_tsz = df * dataPoint.exchange * 0.5 * (3.0*solution.totalSpin_z - solution.totalSpin_z_old) + dataPoint.larmor;

        // Equation of motion for the individual spins
        for(int i = 0; i < tools->energyCount(); i++)
        {
            // Alias for convenience
            real tsz = tmp_tsz + dataPoint.larmorInhomogeneity*tools->energy(i) + df * dataPoint.densityInhomogeneity * exp(-tools->energy(i));
            real tsz2 = tsz*tsz;

            real det = 4.0 + dt2*(tsx2 + tsy2 + tsz2);

            real Mxx = 4.0 + dt2*(tsx2-tsy2-tsz2);
            real Mxy = 2.0*(dt2*tsx*tsy - dt*tsz*2.0);
            real Mxz = 2.0*(dt2*tsx*tsz + dt*tsy*2.0);

            real Myx = 2.0*(dt2*tsx*tsy + dt*tsz*2.0);
            real Myy = 4.0 + dt2*(tsy2-tsx2-tsz2);
            real Myz = 2.0*(dt2*tsy*tsz - dt*tsx*2.0);

            real Mzx = 2.0*(dt2*tsx*tsz - dt*tsy*2.0);
            real Mzy = 2.0*(dt2*tsy*tsz + dt*tsx*2.0);
            real Mzz = 4.0 + dt2*(tsz2-tsx2-tsy2);

            real sx = solution.spin_x[i];
            real sy = solution.spin_y[i];
            real sz = solution.spin_z[i];

            solution.spin_x[i] = (Mxx*sx + Mxy*sy + Mxz*sz)/det;
            solution.spin_y[i] = (Myx*sx + Myy*sy + Myz*sz)/det;
            solution.spin_z[i] = (Mzx*sx + Mzy*sy + Mzz*sz)/det;
        }
    }

    // Value of the total spin at the previous time step

    solution.totalSpin_x_old = solution.totalSpin_x;
    solution.totalSpin_y_old = solution.totalSpin_y;
    solution.totalSpin_z_old = solution.totalSpin_z;

    // Value of the total spin at the current time step

    solution.totalSpin_x = 0;
    solution.totalSpin_y = 0;
    solution.totalSpin_z = 0;

    for(int i = 0; i < tools->energyCount(); i++)
    {
        solution.totalSpin_x += solution.spin_x[i];
        solution.totalSpin_y += solution.spin_y[i];
        solution.totalSpin_z += solution.spin_z[i];
    }

    solution.totalSpin_x /= tools->energyCount();
    solution.totalSpin_y /= tools->energyCount();
    solution.totalSpin_z /= tools->energyCount();


    // Saving the data - not needed for spin echo as we are only interested in the final point

    if(!tools->spinEcho() && t % ((int)std::floor(dataPoint.timeStepCount/tools->timePointsCount())) == 0)
    {
        t = t / std::floor(dataPoint.timeStepCount/tools->timePointsCount());

        solution.totalSpin_x_save[t] = solution.totalSpin_x;
        solution.totalSpin_y_save[t] = solution.totalSpin_y;
        solution.totalSpin_z_save[t] = solution.totalSpin_z;

        for(int i = 0; i < tools->energyBinCount(); i++)
        {
            int j = 0;

            while(j < tools->energyCount() && tools->energy(j) < tools->energyBin(i).first)
                j++;
            while(j < tools->energyCount() && tools->energy(j) <= tools->energyBin(i).second)
            {
                solution.partialSpin_x_save[i][t] += solution.spin_x[j];
                solution.partialSpin_y_save[i][t] += solution.spin_y[j];
                solution.partialSpin_z_save[i][t] += solution.spin_z[j];
                j++;
            }
            solution.partialSpin_x_save[i][t] /= tools->energyBinWeight(i);
            solution.partialSpin_y_save[i][t] /= tools->energyBinWeight(i);
            solution.partialSpin_z_save[i][t] /= tools->energyBinWeight(i);
        }
    }

    solution.timeStep++;
}

// Applies a spin echo (pi pulse) around the y-axis at the latest timestep (at solution.timestep-1). solution.timestep is not increased.
void spinEcho(Solution & solution, SimulationTools * tools)
{
    // Pi pulse around the y-axis sends Sx to -Sx and Sz to -Sz

    for(int i = 0; i < tools->energyCount(); i++)
    {
        solution.spin_x[i] = -solution.spin_x[i];
        solution.spin_z[i] = -solution.spin_z[i];
    }

    solution.totalSpin_x = -solution.totalSpin_x;
    solution.totalSpin_z = -solution.totalSpin_z;

    solution.totalSpin_x_old = -solution.totalSpin_x_old;
    solution.totalSpin_z_old = -solution.totalSpin_z_old;
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

        double sx   = queryResult.record().value("s_x").toDouble();
        double sy   = queryResult.record().value("s_y").toDouble();
        double c    = queryResult.record().value("contrast").toDouble();
        double p    = 0;

        QSqlQuery queryUpdate(tools.database());

        while(true)
        {
            double sxo = sx;
            double syo = sy;
            double co = c;

            if(!queryResult.next()) break;

            sx  = queryResult.record().value("s_x").toDouble();
            sy  = queryResult.record().value("s_y").toDouble();
            c   = queryResult.record().value("contrast").toDouble();

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

    while(queryId.next())
    {
        int simulationId = queryId.record().value(0).toInt();

        QSqlQuery queryResult(tools.database());

        queryResult.prepare("SELECT * FROM result WHERE simulation_id = :i ORDER BY time ASC");
        queryResult.bindValue(":i", simulationId);
        queryResult.exec();

        queryResult.next();

        QSqlQuery queryUpdate(tools.database());

        queryUpdate.prepare("UPDATE result SET phase = 0 WHERE simulation_id = :i and time = 0");
        queryUpdate.bindValue(":i", simulationId);
        queryUpdate.exec();

        double sx   = queryResult.record().value("s_x").toDouble();         // Spin s_x
        double sy   = queryResult.record().value("s_y").toDouble();         // Spin s_y
        double c    = queryResult.record().value("contrast").toDouble();    // Norm
        double p    = 0;                                                    // Phase

        while(true)
        {
            double sxo = sx;
            double syo = sy;
            double co = c;

            if(!queryResult.next()) break;

            sx  = queryResult.record().value("s_x").toDouble();
            sy  = queryResult.record().value("s_y").toDouble();
            c   = queryResult.record().value("contrast").toDouble();

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
