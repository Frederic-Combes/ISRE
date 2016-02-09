#include "tools.h"

#include <QDebug>

using namespace std;

SimulationTools::SimulationTools()
    : p_queueMutex(new mutex()),
      p_energy(NULL),
      p_densityFactor(NULL),
      p_databaseMutex(new mutex())
{
    p_database = QSqlDatabase::addDatabase("QSQLITE", "isre-result");
    p_database.setDatabaseName("ISRE-RESULT");
    p_database.open();
}

SimulationTools::~SimulationTools()
{
    delete p_queueMutex;
    delete p_databaseMutex;

    delete [] p_energy;
    delete [] p_densityFactor;
}

void SimulationTools::init()
{
    // Compute the sampled energies

    p_energy = new double[p_energyCount];

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
            p_energy[i] = p_energy[i-1] + exp(p_energy[i-1])/(p_energy[i-1]*real(p_energyCount));
        break;
    case 3:
        p_energy[0] = 0;
        p_energy[1] = pow(8.0/p_energyCount, 1.0/3.0);
        for(int i = 2; i < p_energyCount; i++)
            p_energy[i] = p_energy[i-1] + real(2.0)*exp(p_energy[i-1])/(p_energy[i-1]*p_energy[i-1]*real(p_energyCount));
        break;
    default:
        qDebug() << "Supported space dimensions are 1, 2 and 3 - this message should never appear.";
        return;
    }

    // Compute the energy bin weight (the number of points inside each bin

    for(auto it = p_energyBin.begin(); it != p_energyBin.end(); it++)
    {
        int weight = 0,  i = 0;

        while(p_energy[i] < it->first)
            i++;
        while(i < p_energyCount && p_energy[i] <= it->second)
        {
            i++;
            weight++;
        }

        p_energyBinWeight.push_back(weight);
    }

    // Compute the step at wich one must make a spin echo
    std::sort(p_tmp_spinEchoTime.begin(), p_tmp_spinEchoTime.end());

    if(p_spinEcho)
    {
        for(auto it = p_tmp_spinEchoTime.begin(); it != p_tmp_spinEchoTime.end(); it++)
            p_spinEchoTimeStep.push_back((int((*it) * p_timeStepsCount)));
    }
    else
    {
        // TODO: Write the correct formula
        for(auto it = p_tmp_spinEchoTime.begin(); it != p_tmp_spinEchoTime.end(); it++)
            p_spinEchoTimeStep.push_back((int((*it) * p_timeStepsCount / p_totalTime)));
    }

    p_tmp_spinEchoTime.clear();
}
