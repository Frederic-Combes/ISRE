#ifndef SCRIPTOBJECTS_H
#define SCRIPTOBJECTS_H

#include "result.h"

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
#include <QMetaType>
class QScriptContext;

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
    inline double saveInterval() const  {return p_saveInterval;}
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

class Simulation
{
public:
    Simulation();

    void setName(QString name) {p_name = name;}
    QString name() const {return p_name;}

    void setCode(int index, QString code);
    QString code(int index) const;

private:
    QString                     p_name;

    std::vector<QString>        p_code;
};

// NOTE: not really a "script object" as compared to ScriptObject::Parameter or ScriptObject::Function ...
class Result
{
    // The engine globalObject().data() stores a db connection
public:
    static QScriptValue scriptFunctionConstructor(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionMatches(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionExport(QScriptContext * context, QScriptEngine * engine);
    //static QScriptValue scriptFunctionSave(QScriptContext * context, QScriptEngine * engine);
    //static QScriptValue scriptFunctionErase(QScriptContext * context, QScriptEngine * engine);

public:
    Result();
    Result(unsigned int id, std::shared_ptr<SimulationResult> simulationResult, std::shared_ptr<SimulationDescription> simulationDescription);
    Result(const Result & result) = default;

    void setMetadata(QString name, QString value);
    std::pair<bool, QString> metadata(QString name) const;

private:
    /**
     * @brief Database ID of the simulation
     */
    unsigned int p_id;

    /**
     * @brief Simulation results wrapped by the script object
     */
    std::shared_ptr<SimulationResult> p_simulationResult;

    /**
     * @brief Metadata of the result
     */
    std::shared_ptr<SimulationDescription> p_simulationDescription;
};

} // namespace ScriptObject

Q_DECLARE_METATYPE(ScriptObject::Result)

#endif // SCRIPTOBJECTS_H
