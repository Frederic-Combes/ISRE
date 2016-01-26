#include "tools.h"

#include <QDebug>

using namespace std;

bool operator<(const DataPoint & a, const DataPoint & b)
{
    if(a.larmor < b.larmor) return true;
    else if( a.larmor > b.larmor) return false;

    if(a.larmorInhomogeneity < b.larmorInhomogeneity) return true;
    if(a.larmorInhomogeneity > b.larmorInhomogeneity) return false;

    if(a.exchange < b.exchange) return true;
    else if(a.exchange > b.exchange) return false;

    if(a.damping < b.damping) return true;
    else return false;
}

ThreadTools::ThreadTools()
    : p_queueMutex(new mutex()),
      p_fileMutex(new mutex()),
      p_resultMutex(new mutex()),
      p_energy(NULL),
      p_densityFactor(NULL)
{
    p_database = QSqlDatabase::addDatabase("QSQLITE", "isre-result");
    p_database.setDatabaseName("ISRE-RESULT");
    p_database.open();
}

ThreadTools::ThreadTools(const ThreadTools &)
{}

ThreadTools::~ThreadTools()
{
    p_file.close();

    delete p_queueMutex;
    delete p_fileMutex;
    delete p_resultMutex;

    delete [] p_energy;
    delete [] p_densityFactor;

}

void ThreadTools::init()
{
    /* Compute the sampled energies */

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

    /* Compute the energy bin weight (the number of points inside each bin */

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

    /* Compute the step at wich one must make a spin echo */
    std::sort(p_tmp_spinEchoTime.begin(), p_tmp_spinEchoTime.end());

    for(auto it = p_tmp_spinEchoTime.begin(); it != p_tmp_spinEchoTime.end(); it++)
        p_spinEchoTimeStep.push_back(((int) ((*it) * p_timeStepsCount)));

    p_tmp_spinEchoTime.clear();
}

bool ThreadTools::setOutput(QString filename)
{
    p_file.open(filename.toLocal8Bit().data(), ios::trunc | ios::out);

    if(!p_file.good())
    {
        qDebug() << "Could not open the file" << filename;
        return false;
    }

    return true;
}

void ThreadTools::insertResultPoint(const DataPoint & dataPoint, const Solution & solution)
{
    SpinEchoResult r;
    r.time = dataPoint.totalTime;
    r.totalSpin_x = solution.totalSpin_x;
    r.totalSpin_y = solution.totalSpin_y;
    r.totalSpin_z = solution.totalSpin_z;

    for(unsigned int i = 0; i < p_energyBin.size(); i++)
    {
        double sx = 0, sy = 0, sz = 0;

        int j = 0;

        while(p_energy[j] < p_energyBin.at(i).first)
            j++;
        while(p_energy[j] < p_energyBin.at(i).second)
        {
            sx += solution.spin_x[j];
            sy += solution.spin_y[j];
            sz += solution.spin_z[j];
            j++;
        }

        r.partialSpin_x.push_back(sx / p_energyBinWeight.at(i));
        r.partialSpin_y.push_back(sy / p_energyBinWeight.at(i));
        r.partialSpin_z.push_back(sz / p_energyBinWeight.at(i));
    }

    if(p_spinEchoResults.count(dataPoint) == 0)
    {
        p_spinEchoResults.insert(make_pair(dataPoint, vector<SpinEchoResult>()));

        SpinEchoResult intialPoint;

        intialPoint.time = 0;
        intialPoint.totalSpin_x = 1;
        intialPoint.totalSpin_y = 0;
        intialPoint.totalSpin_z = 0;

        for(unsigned int i = 0; i < p_energyBin.size(); i++)
        {
            intialPoint.partialSpin_x.push_back(1);
            intialPoint.partialSpin_y.push_back(0);
            intialPoint.partialSpin_z.push_back(0);
        }

        p_spinEchoResults.at(dataPoint).push_back(intialPoint);
    }

    p_spinEchoResults.at(dataPoint).push_back(r);
}

bool ThreadTools::allResultPointsComputed(const DataPoint & dataPoint) const
{
    if(p_spinEchoResults.find(dataPoint) == p_spinEchoResults.end()) return false;

    return (p_spinEchoResults.at(dataPoint).size() == (unsigned int)p_timePointsCount + 1);
}

vector<ThreadTools::SpinEchoResult> ThreadTools::result(const DataPoint & dataPoint)
{
    sort(p_spinEchoResults[dataPoint].begin(), p_spinEchoResults[dataPoint].end(), spinEchoResultTimeOrder);

    return p_spinEchoResults[dataPoint];
}

void ThreadTools::removeResult(const DataPoint & dataPoint)
{
    p_spinEchoResults.erase(dataPoint);
}
