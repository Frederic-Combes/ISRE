#include "tools.h"

#include <QDebug>

using namespace std;

Solution::Solution(const SimulationTools & simulation)
{
    spin_x.resize(simulation.energyCount());
    spin_y.resize(simulation.energyCount());
    spin_z.resize(simulation.energyCount());
    spin_n.resize(simulation.energyCount());

    totalSpin_x_save        = new double[simulation.timePointsCount()+1];
    totalSpin_y_save        = new double[simulation.timePointsCount()+1];
    totalSpin_z_save        = new double[simulation.timePointsCount()+1];
    spectatorPopulationF1_save  = new double[simulation.timePointsCount()+1];
    spectatorPopulationF2_save  = new double[simulation.timePointsCount()+1];

    partialSpin_x_save      = new double*[simulation.energyBinCount()+1];
    partialSpin_y_save      = new double*[simulation.energyBinCount()+1];
    partialSpin_z_save      = new double*[simulation.energyBinCount()+1];

    for(int i = 0; i < simulation.energyBinCount(); ++i)
    {
        partialSpin_x_save[i] = new double[simulation.timePointsCount()+1];
        partialSpin_y_save[i] = new double[simulation.timePointsCount()+1];
        partialSpin_z_save[i] = new double[simulation.timePointsCount()+1];
    }

    p_energyBinCount = simulation.energyBinCount();
}

Solution::~Solution()
{
    delete [] totalSpin_x_save;
    delete [] totalSpin_y_save;
    delete [] totalSpin_z_save;
    delete [] spectatorPopulationF1_save;
    delete [] spectatorPopulationF2_save;

    for(unsigned int i = 0; i < p_energyBinCount; ++i)
    {
        delete [] partialSpin_x_save[i];
        delete [] partialSpin_y_save[i];
        delete [] partialSpin_z_save[i];
    }

    delete [] partialSpin_x_save;
    delete [] partialSpin_y_save;
    delete [] partialSpin_z_save;


}

SimulationTools::SimulationTools()
    : p_queueMutex(new mutex()),
      p_dimension(0),
      p_energyCount(0),
      p_energy(NULL),
      p_hasMovingOperations(false),
      p_databaseMutex(new mutex())
{
    p_database = QSqlDatabase::addDatabase("QSQLITE", "isre-result");
    p_database.setDatabaseName("results.db");
    p_database.open();
}

SimulationTools::~SimulationTools()
{
    delete p_queueMutex;
    delete p_databaseMutex;

    delete [] p_energy;
}

void SimulationTools::init()
{
    // Compute the step at wich one must make a spin echo
    sort(p_tmp_spinEchoTime.begin(), p_tmp_spinEchoTime.end());

    if(p_spinEcho) {
        for(auto it = p_tmp_spinEchoTime.begin(); it != p_tmp_spinEchoTime.end(); ++it) {
            p_spinEchoTimeStep.push_back(*it * double(p_timeStepsCount));
        }
    } else {
        for(auto it = p_tmp_spinEchoTime.begin(); it != p_tmp_spinEchoTime.end(); ++it) {
            p_spinEchoTimeStep.push_back(*it * double(p_timeStepsCount) / p_totalTime);
        }
    }

    p_tmp_spinEchoTime.clear();
}

void SimulationTools::setSpaceDimension(int dimension)
{
    if(p_dimension != dimension) {
        p_dimension = dimension;

        sampleEnergy(dimension, p_energyCount);
    }
}

void SimulationTools::sampleEnergy(int spaceDimension, int sampleCount)
{
    p_dimension = spaceDimension;

    if(sampleCount <= 0) {
        return;
    }

    if(p_energy != NULL) {
        delete [] p_energy;
    }

    p_energyCount   = sampleCount;
    p_energy        = new double[p_energyCount];

    switch(p_dimension)
    {
    case 1:
        for(int i = 0; i < p_energyCount; i++)
            p_energy[i] = -log(1.0-double(i)/double(p_energyCount));
        break;
    case 2:
        p_energy[0] = 0;
        p_energy[1] = pow(2.0/p_energyCount, 0.5);
        for(int i = 2; i < p_energyCount; i++)
            p_energy[i] = p_energy[i-1] + exp(p_energy[i-1])/(p_energy[i-1]*double(p_energyCount));
        break;
    case 3:
        p_energy[0] = 0;
        p_energy[1] = pow(8.0/p_energyCount, 1.0/3.0);
        for(int i = 2; i < p_energyCount; i++)
            p_energy[i] = p_energy[i-1] + double(2.0)*exp(p_energy[i-1])/(p_energy[i-1]*p_energy[i-1]*double(p_energyCount));
        break;
    default:
        qDebug() << "Supported space dimensions are 1, 2 and 3 - this message should never appear.";
        return;
    }

    computeEnergyBinWeights();
}

void SimulationTools::computeEnergyBinWeights()
{
    for(auto it = p_energyBin.begin(); it != p_energyBin.end(); ++it) {
        int weight = 0;
        int i = 0;

        while(i < p_energyCount && p_energy[i] < it->first) {
            ++i;
        }
        while(i < p_energyCount && p_energy[i] <= it->second) {
            ++i;
            ++weight;
        }

        p_energyBinWeight.push_back(weight);
    }
}

void SimulationTools::addOperation(const Operation & operation)
{
    p_hasMovingOperations |= !operation.isFixedTime;

    p_operations.push_back(operation);
}

const Operation * SimulationTools::nextOperation(double currentTime, double totalTime) const
{
    double min = totalTime;
    const Operation * next = NULL;

    for(const Operation & operation : p_operations) {
        double deltaTime = operation.isFixedTime ? operation.time - currentTime : operation.time * totalTime - currentTime;

        if(deltaTime > 0 && deltaTime < min) {
            min = deltaTime;
            next = &operation;
        }
     }

    return next;
}
