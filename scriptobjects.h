#ifndef SCRIPTOBJECTS_H
#define SCRIPTOBJECTS_H

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/StdVector>
#pragma GCC diagnostic pop

// Script
#include <QScriptValue>

namespace ScriptObject {

class Operation
{
public:
    static Operation fromScriptValue(QScriptValue value);

public:
    alignas(16) Eigen::Matrix3d matrix;
    double time;
    bool   isFixedTime;
};

class Parameter
{
public:
    static Parameter fromScriptValue(QScriptValue value);

public:
    Parameter() : name(), values(), index(0) {}
    Parameter(QString n, std::vector<double> v) : name(n), values(v), index(0) {}

    QString             name;
    std::vector<double> values;
    unsigned int        index;  // TODO : move out of here, it is not part of the script
};

class Function
{
public:
    static Function fromScriptValue(QScriptValue value);

public:
    Function() {}
    Function(QString n, QString c) : name(n), code(c) {}

    QString name;
    QString code;
};

class Settings
{
public:
    static double defaultDensityEnergyDependenceDimension1(double energy);
    static double defaultDensityEnergyDependenceDimension2(double energy);
    static double defaultDensityEnergyDependenceDimension3(double energy);

public:
    static Settings fromScriptValue(QScriptValue value);

public:
    typedef std::pair<unsigned int, unsigned int> EnergyBin;

    inline int dimension() const {return p_dimension;}

    inline double duration() const      {return p_duration;}
    inline double timestep() const      {return p_duration / p_timestepCount;}
    inline int timestepCount() const    {return p_timestepCount;}
    inline double saveInterval() const        {return p_saveInterval;}
    inline int saveCount() const        {return std::ceil(p_duration/p_saveInterval);}

    const std::vector<double> & energySamples() const {return p_energySamples;}
    const std::vector<EnergyBin> & energyBins() const {return p_energyBins;}

private:
    void sampleEnergy(int sampleCount);
    void addEnergyBin(double energyMin, double energyMax);

private:
    int                         p_dimension;

    double                      p_duration;
    int                         p_timestepCount;
    double                      p_saveInterval;

    std::vector<double>         p_energySamples;
    std::vector<EnergyBin>      p_energyBins;
};

} // namespace ScriptObject

#endif // SCRIPTOBJECTS_H
