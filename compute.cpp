#include "compute.h"

#include "strings.h"
#include "result.h"

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
using namespace std;
using namespace ScriptObject;
constexpr double PI = M_PI;

void threadLoop(SimulationContext & context)
{
    double threadID = double(hash<thread::id>()(this_thread::get_id())); Q_UNUSED(threadID);

    // Aliases
    const Settings & settings = context.settings();

    // Memory pre-allocation for the solution
    Solution solution(settings);

    unsigned int simulationIndex = 0;

    while(context.next(simulationIndex)) {
        // Obtaining the data to process
        double * data = context.data(simulationIndex);
        // Operations are returned in a time-ordered fashion
        SimulationContext::Operations operations = context.operations(simulationIndex);
        // Where to store the simulation results
        const vector<ResultPtr> & results = context.results(simulationIndex);

        context.evaluateFunctions(ComputeFor::Initialization,   Compute::AtStart, data);
        context.evaluateFunctions(ComputeFor::Step,             Compute::AtStart, data);
        context.evaluateFunctions(ComputeFor::Operation,        Compute::AtStart, data);

        // Numerical solution
        initialize(data, solution, context, results);

        for(const ScriptObject::Operation & operation : operations) {
            int timestep = operation.isFixedTime ?
                        round(operation.time * data[Index::TimestepCount] / data[Index::Duration])
                      : round(operation.time * data[Index::TimestepCount]);

            while(solution.timeStep <= timestep) {
                step(data, solution, context, results);
            }

            // It's guaranted that the last operation is the identity happening at time data[Index::Duration]
            applyOperation(operation, solution, settings);
        }

        // Getting database IDs
        //vector<int> simulationIds = context.ids(simulationIndex);

        // Saving data
        /*context.databaseLock().lock();

        // Starting the SQL transaction
        db.transaction();

        // Prepare the queries
        QSqlQuery sqlSimulationResult(db);
        sqlSimulationResult.prepare("INSERT INTO result (id, time, s_x, s_y, s_z, contrast, lost_1, lost_2, phase) VALUES (:s, :t, :x, :y, :z, :c, :l1, :l2, :p);");

        QSqlQuery sqlEnergyBinResult(db);
        sqlEnergyBinResult.prepare("INSERT INTO result (id, time, s_x, s_y, s_z, contrast, phase) VALUES (:i, :t, :x, :y, :z, :c, :p);");

        for(unsigned int t = 0; t < solution.save.size(); ++t) {
            // Time we are saving the values for
            double time = solution.save(t).timestep * data[Index::Duration] / data[Index::TimestepCount];

            // Total Spin
            sqlSimulationResult.bindValue(":s",  simulationIds.back());
            sqlSimulationResult.bindValue(":t",  time);
            sqlSimulationResult.bindValue(":x",  solution.save(t).total.spin(0));
            sqlSimulationResult.bindValue(":y",  solution.save(t).total.spin(1));
            sqlSimulationResult.bindValue(":z",  solution.save(t).total.spin(2));
            // NOTE: Experimentally measured contrast is sqrt(s_x² + s_y²)
            sqlSimulationResult.bindValue(":c",  solution.save(t).total.spin.head<3>().norm());
            sqlSimulationResult.bindValue(":l1", solution.save(t).total.atomLosses(1));
            sqlSimulationResult.bindValue(":l2", solution.save(t).total.atomLosses(2));
            sqlSimulationResult.bindValue(":p",  solution.save(t).total.phase);

            sqlSimulationResult.exec();

            // Partial spin
            for(unsigned int i = 0; i < settings.energyBins().size(); ++i) {
                sqlEnergyBinResult.bindValue(":i", simulationIds.at(i));
                sqlEnergyBinResult.bindValue(":t", time);
                sqlEnergyBinResult.bindValue(":x", solution.save(t).partial(i).spin(0));
                sqlEnergyBinResult.bindValue(":y", solution.save(t).partial(i).spin(1));
                sqlEnergyBinResult.bindValue(":z", solution.save(t).partial(i).spin(2));
                sqlEnergyBinResult.bindValue(":c", solution.save(t).partial(i).spin.head<3>().norm());
                sqlEnergyBinResult.bindValue(":p", solution.save(t).partial(i).phase);

                sqlEnergyBinResult.exec();
            }
        }

        // If the simulation has no moving operations, then at this point it is complete and we can remove the
        // marker.
        if(context.allOperationsAreFixedTime()) {
            QSqlQuery sqlUpdateDescription(db);
            sqlUpdateDescription.prepare("DELETE FROM description WHERE id = :id and name = :n");
            sqlUpdateDescription.bindValue(":id", simulationIds.back());
            sqlUpdateDescription.bindValue(":n", Strings::Database::Description::SimulationIncomplete);
            sqlUpdateDescription.exec();
        }

        // Commit all the SQL data
        db.commit();

        context.databaseLock().unlock();*/

        delete [] data;
    }
}

// Initialize the first and second time step of the solution, as the computing the solution at time t requires
// to know the solution at time t-1 and t-2
void initialize(double * data, Solution & solution, SimulationContext & context, const vector<ResultPtr> & results)
{
    // Alias
    const Settings & settings = context.settings();

    // Zero-ing
    solution.clear();

    // Values at time 0
    {
        solution.timeStep       = 0;
        data[Index::Timestep]   = 0;
        data[Index::Time]       = 0;

        context.evaluateFunctions(ComputeFor::Initialization, Compute::OnceForEachStep, data);

        auto spin_it = solution.current.spin.begin();
        for(unsigned int i = 0; i < settings.energySamples().size(); ++i, ++spin_it) {
            Vector4d & spin = *spin_it;

            data[Index::EnergyIndex]    = i;
            data[Index::Energy]         = settings.energySamples().at(i);

            context.evaluateFunctions(ComputeFor::Initialization, Compute::OnceForEachEnergy, data);

            spin << data[Index::InitialSpinX],
                    data[Index::InitialSpinY],
                    data[Index::InitialSpinZ],
                    1.0;
        }

        // Save the values at the previous timestep
        solution.previous.total.spin  = solution.current.total.spin;
        solution.previous.total.phase = solution.current.total.phase;

        // Update the current values
        solution.current.total.spin  = accumulate(solution.current.spin.begin(), solution.current.spin.end(), Vector4d::Zero().eval()) / solution.current.spin.size();
        solution.current.total.phase = 0;
        solution.current.total.atomLosses(1) = 0;
        solution.current.total.atomLosses(2) = 0;

        // Partial spins
        for(unsigned int i = 0; i < settings.energyBins().size(); ++i) {
            // Save the values at the previous timestep
            solution.previous.partial(i).spin  = solution.current.partial(i).spin;
            solution.previous.partial(i).phase = solution.current.partial(i).phase;
            // Update the current values
            const Settings::EnergyBin & bin = settings.energyBins().at(i);
            solution.current.partial(i).spin = accumulate(solution.current.spin.begin() + bin.first,
                                                          solution.current.spin.begin() + bin.second,
                                                          Vector4d::Zero().eval()
                                                          ) / (bin.second - bin.first);
            solution.current.partial(i).phase = 0;
        }
    }

    // Saves the values if necessary
    save(data, solution, settings, results);

    // Values at the first time step
    {
        // Duration of the time steps in seconds
        double dt = data[Index::Duration] / data[Index::TimestepCount];

        // Timestep
        solution.timeStep       = 1;
        data[Index::Timestep]   = 1;
        data[Index::Time]       = dt;

        context.evaluateFunctions(ComputeFor::Initialization, Compute::OnceForEachStep, data);

        for(unsigned int i = 0; i < settings.energySamples().size(); ++i) {
            Vector4d & spin = solution.current.spin.at(i);
            Vector4d & s    = solution.current.total.spin;

            data[Index::EnergyIndex]    = i;
            data[Index::Energy]         = settings.energySamples().at(i);

            context.evaluateFunctions(ComputeFor::Initialization, Compute::OnceForEachEnergy, data);

            Matrix4d RF, RE, RD, RL;

            // Pseudo magnetic field
            RF <<   0,                      -data[Index::Frequency],    0, 0,
                    data[Index::Frequency],  0,                         0, 0,
                    0,                       0,                         0, 0,
                    0,                       0,                         0, 0;

            // Exchange
            RE <<    0,                                 -data[Index::Exchange]*s(2)*s(3),    data[Index::Exchange]*s(1)*s(3),   0,
                     data[Index::Exchange]*s(2)*s(3),    0,                                 -data[Index::Exchange]*s(0)*s(3),   0,
                    -data[Index::Exchange]*s(1)*s(3),    data[Index::Exchange]*s(0)*s(3),    0,                                 0,
                     0,                                  0,                                  0,                                 0;
            // Damping
            RD <<   -data[Index::Damping],  0,                       0,                     0,
                     0,                    -data[Index::Damping],    0,                     0,
                     0,                     0,                      -data[Index::Damping],  0,
                     0,                     0,                       0,                     0;
            // Atom losses
            double dfls = 0.25 * (data[Index::AtomLossFromF1] + data[Index::AtomLossFromF2]);    // Total atom losses to spectator states
            double dfld = 0.25 * (data[Index::AtomLossFromF1] - data[Index::AtomLossFromF2]);    // Difference of atom losses to spectator states
            RL <<   0, 0,  0,                        0,
                    0, 0,  0,                        0,
                    0, 0, -(dfls*s(3) + dfld*s(2)), -(dfld*s(3) + dfls*s(2)),
                    0, 0, -(dfld*s(3) + dfls*s(2)), -(dfls*s(3) + dfld*s(2));

            spin = ((Matrix4d::Identity() + dt * (RF + RE + RD + RL)) * spin).eval();
        }

        // Save the values at the previous timestep
        solution.previous.total.spin  = solution.current.total.spin;
        solution.previous.total.phase = solution.current.total.phase;

        // Update current values
        solution.current.total.spin   = accumulate(solution.current.spin.begin(), solution.current.spin.end(), Vector4d::Zero().eval()) / solution.current.spin.size();
        solution.current.total.phase += asin((   solution.current.total.spin(0) * solution.previous.total.spin(1)
                                               - solution.previous.total.spin(0) * solution.current.total.spin(1))
                                             / (solution.current.total.spin.head<2>().norm() * solution.previous.total.spin.head<2>().norm()));
        solution.current.total.atomLosses(1) = dt * data[Index::AtomLossFromF1];
        solution.current.total.atomLosses(2) = dt * data[Index::AtomLossFromF2];

        // Partial spins
        for(unsigned int i = 0; i < settings.energyBins().size(); ++i) {
            // Save the values at the previous timestep
            solution.previous.partial(i).spin  = solution.current.partial(i).spin;
            solution.previous.partial(i).phase = solution.current.partial(i).phase;
            // Update the current values
            const Settings::EnergyBin & bin = settings.energyBins().at(i);
            solution.current.partial(i).spin = accumulate(solution.current.spin.begin() + bin.first,
                                                          solution.current.spin.begin() + bin.second,
                                                          Vector4d::Zero().eval()
                                                          ) / (bin.second - bin.first);
            solution.current.partial(i).phase += asin((  solution.current.partial(i).spin(0) * solution.previous.partial(i).spin(1)
                                                       - solution.previous.partial(i).spin(0) * solution.current.partial(i).spin(1)
                                                      ) / (solution.current.partial(i).spin.head<2>().norm() * solution.previous.partial(i).spin.head<2>().norm()));
        }

    }

    // Saves the values if necessary
    save(data, solution, settings, results);

    // Next timestep
    solution.timeStep = 2;
}

// Given the solution of the equation up to the time step solution.timestep - 1, computes the solution at the sime step solution.timestep and increases timestep by one
void step(double * data, Solution & solution, const SimulationContext & context, const std::vector<ResultPtr> & results)
{
    // The equation is result = (2-Re dt)^(-1) (2 + Re dt) se + (2-Re dt)^(-1) G s

    // Total Spin
    Vector4d s = 0.5*(3.0*solution.current.total.spin - solution.previous.total.spin);

    // Aliases
    const Settings & settings = context.settings();
    double dt = data[Index::Duration] / data[Index::TimestepCount];

    data[Index::Timestep]   = solution.timeStep;
    data[Index::Time]       = solution.timeStep * data[Index::Duration] / data[Index::TimestepCount];

    context.evaluateFunctions(ComputeFor::Step, Compute::OnceForEachStep, data);

    // test
    double h = accumulate(solution.current.spin.begin(), solution.current.spin.end(), 0.0, [](double old, const Vector4d & s) {return old + abs(s(2));}) / solution.current.spin.size();
    (void)h;

    // Solving the equation of motion for the individual spins
    for(unsigned int i = 0; i < settings.energySamples().size(); ++i) {
        Vector4d & spin = solution.current.spin.at(i);

        data[Index::EnergyIndex] = i;
        data[Index::Energy] = settings.energySamples().at(i);

        context.evaluateFunctions(ComputeFor::Step, Compute::OnceForEachEnergy, data);

        // Matrix initialization
        Matrix4d R, RD;
        // R = RF + RE + RD + RL;

        // Pseudo magnetic field
        /* 0,                      -data[Index::Frequency],    0, 0,
         * data[Index::Frequency],  0,                         0, 0,
         * 0,                       0,                         0, 0,
         * 0,                       0,                         0, 0;
         */

        R.setZero();
        R(0,1) += -data[Index::Frequency];
        R(1,0) += +data[Index::Frequency];

        // Exchange
        /*  0,                                 -data[Index::Exchange]*s(2)*s(3),    data[Index::Exchange]*s(1)*s(3),   0,
         *  data[Index::Exchange]*s(2)*s(3),    0,                                 -data[Index::Exchange]*s(0)*s(3),   0,
         * -data[Index::Exchange]*s(1)*s(3),    data[Index::Exchange]*s(0)*s(3),    0,                                 0,
         *  0,                                  0,                                  0,                                 0;
         */
        R(0,1) += -data[Index::Exchange]*s(2)*s(3);
        R(0,2) += +data[Index::Exchange]*s(1)*s(3);
        R(1,0) += +data[Index::Exchange]*s(2)*s(3);
        R(1,2) += -data[Index::Exchange]*s(0)*s(3);
        R(2,0) += -data[Index::Exchange]*s(1)*s(3);
        R(2,1) += +data[Index::Exchange]*s(0)*s(3);

        // Damping
        /* -data[Index::Damping],  0,                       0,                     0,
         *  0,                    -data[Index::Damping],    0,                     0,
         *  0,                     0,                      -data[Index::Damping],  0,
         *  0,                     0,                       0,                     0;
         */
        RD.setZero();
        RD(0,0) += -data[Index::Damping];
        RD(1,1) += -data[Index::Damping];
        RD(2,2) += -data[Index::Damping];

        R += RD;

        // Atom losses
        double dfls = 0.25 * (data[Index::AtomLossFromF1] + data[Index::AtomLossFromF2]);    // Total atom losses to spectator states
        double dfld = 0.25 * (data[Index::AtomLossFromF1] - data[Index::AtomLossFromF2]);    // Difference of atom losses to spectator states
        /* 0, 0,  0,                        0,
         * 0, 0,  0,                        0,
         * 0, 0, -(dfls*s(3) + dfld*s(2)), -(dfld*s(3) + dfls*s(2)),
         * 0, 0, -(dfld*s(3) + dfls*s(2)), -(dfls*s(3) + dfld*s(2));*/
        R(2,2) += -(dfls*s(3) + dfld*s(2));
        R(2,3) += -(dfld*s(3) + dfls*s(2));
        R(3,2) += -(dfld*s(3) + dfls*s(2));
        R(3,3) += -(dfls*s(3) + dfld*s(2));

        // Computation

        Matrix4d RT; // = 2 + R * dt
        Matrix4d RI; // = (2 - R * dt).inverse();

        // We can't inverse in place, so we use RT as a temporary to compute RI
        RT = 2.0 * Matrix4d::Identity() - dt * R;
        RI = RT.inverse();
        RT = 2.0 * Matrix4d::Identity() + dt * R;

        spin = RI * (RT * spin - dt * RD * s).eval();
    }

    // Update values at the previous timestep
    solution.previous.total.spin  = solution.current.total.spin;
    solution.previous.total.phase = solution.current.total.phase;

    // Update current values
    solution.current.total.spin = accumulate(solution.current.spin.begin(), solution.current.spin.end(), Vector4d::Zero().eval()) / solution.current.spin.size();

    // TODO : Move inside the loop
    solution.current.total.atomLosses(1) += dt * data[Index::AtomLossFromF1]* 0.5 * (s(3) + s(2));
    solution.current.total.atomLosses(2) += dt * data[Index::AtomLossFromF2]* 0.5 * (s(3) - s(2));

    // Partial spins
    for(unsigned int i = 0; i < settings.energyBins().size(); ++i) {
        // Save the values at the previous timestep
        solution.previous.partial(i).spin  = solution.current.partial(i).spin;
        solution.previous.partial(i).phase = solution.current.partial(i).phase;
        // Update the current values
        const Settings::EnergyBin & bin = settings.energyBins().at(i);
        solution.current.partial(i).spin = accumulate(solution.current.spin.begin() + bin.first,
                                                      solution.current.spin.begin() + bin.second,
                                                      Vector4d::Zero().eval()
                                                      ) / (bin.second - bin.first);
    }

    // Computes the phase of the spins
    phase(solution, settings);

    // Save the data if necessary
    save(data, solution, settings, results);

    ++solution.timeStep;
}

void save(double * data, Solution & solution, const Settings & settings, const vector<ResultPtr> & results)
{
    // Saving the data if needed
    if(solution.timeStep < static_cast<int>(data[Index::NextSave])) {
        // Do nothing
    } else {
        solution.save(data[Index::Time], results);

        data[Index::NextSave] = (round(data[Index::Time] / settings.saveInterval()) + 1) * settings.saveInterval() * data[Index::TimestepCount] / data[Index::Duration];
        data[Index::NextSave] = min(data[Index::NextSave], data[Index::TimestepCount]);
    }
}

// Computes the new phase of the spin (handles the 2-pi periodicity)
void phase(Solution & solution, const Settings & settings)
{
    double p = atan2(solution.current.total.spin(1), solution.current.total.spin(0));

    solution.current.total.phase = p - 2 * M_PI * round( (p - solution.previous.total.phase) / ( 2 * M_PI ) );

    if(abs(solution.current.total.phase - solution.previous.total.phase) / M_PI > 0.25) {
        qDebug() << "phase: large phase step at time" << solution.timeStep;
        qDebug() << "new =" << atan2(solution.current.total.spin(1), solution.current.total.spin(0)) / M_PI;
        qDebug() << "old =" << solution.previous.total.phase / M_PI;
        qDebug() << "sx" << solution.current.total.spin(0);
        qDebug() << "sy" << solution.current.total.spin(1);
    }

    for(unsigned int i = 0; i < settings.energyBins().size(); ++i) {
        double p = atan2(solution.current.partial(i).spin(1), solution.current.partial(i).spin(0));

        solution.previous.partial(i).phase = p - 2 * M_PI * round( (p - solution.previous.partial(i).phase) / ( 2 * M_PI ) );
    }
}

// Applies an operation at the latest timestep (at solution.timestep-1)
// solution.timestep is not increased
void applyOperation(const Operation & operation, Solution & solution, const Settings & settings)
{
    for(Vector4d & spin : solution.current.spin) {
       spin.head<3>() = operation.matrix * spin.head<3>().eval();
    }

    solution.current.total.spin.head<3>()  = operation.matrix * solution.current.total.spin.head<3>().eval();
    solution.previous.total.spin.head<3>() = operation.matrix * solution.previous.total.spin.head<3>().eval();

    // Update the phases (set them in the [-pi;+pi] range)
    solution.current.total.phase = atan2(solution.current.total.spin(1), solution.current.total.spin(0));

    for(unsigned int i = 0; i < settings.energyBins().size(); ++i) {
        solution.current.partial(i).spin.head<3>()  = operation.matrix * solution.current.partial(i).spin.head<3>().eval();

        solution.current.partial(i).phase = atan2(solution.current.partial(i).spin(1), solution.current.partial(i).spin(0));
    }
}

