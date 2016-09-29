#ifndef SIMULATION_H
#define SIMULATION_H

#include <utility>
#include <vector>
#include <map>
#include <thread>
#include <mutex>

// Database support
#include <QString>
#include <QSqlDatabase>

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/StdVector>
#pragma GCC diagnostic pop

#include "script.h"
#include "result.h"

#include <QScriptValue>
#include "scriptobjects.h"

class SimulationContext
{
public:
    typedef std::vector<std::reference_wrapper<const ScriptObject::Operation>> Operations;
public:
    SimulationContext();
    ~SimulationContext();

    void setOverwrite(bool overwrite = false);

    void setSettings(ScriptObject::Settings settings);
    void setSimulation(const ScriptObject::Simulation & simulation);
    void addParameter(ScriptObject::Parameter parameter);
    void addFunction(ScriptObject::Function function);
    void addOperation(ScriptObject::Operation operation);

    inline const ScriptObject::Settings & settings() const;
    inline const ScriptObject::Simulation & simulation() const;
    inline bool allOperationsAreFixedTime() const;

    const CompiledFunctions & functions(unsigned int index = 0) const;

    void evaluateFunctions(unsigned int index, int flag, double * data) const;

    int getIndex(double * & data);

    void prepare();
    std::vector<ScriptObject::Result> end();

    // TODO : auto release / regain the data ptr (as well as the Operations vector)!
    // HINT: use unique_ptr<double[], allocator> ?
    bool next(unsigned int & simulationIndex);          // Set simulationIndex to the next simulation index, and return true. Return false if there are no more simulations to run
    double * data(int simulationIndex) const;           // Returns the data needed for the simulation indexed by simulationIndex
    Operations operations(int simulationIndex) const;   // Returns the time-ordered list of operations needed for the simulation indexed by simulationIndex
    const std::vector<std::shared_ptr<SimulationResult>> & results(unsigned int simulationIndex) const;

private:
    void prepareFunctions(mathpresso::Context & context);
    void prepareResults();

    CompiledFunctions functionsRequiredToCompute(std::vector<QString> functionNames) const;
    CompiledFunctions functionsRequiredToCompute(const QString & functionName) const;

private:
    bool                                    p_overwrite;

    ScriptObject::Settings                  p_settings;
    ScriptObject::Simulation                p_simulation;
    std::vector<ScriptObject::Parameter>    p_parameters;
    std::vector<ScriptObject::Function>     p_functions;
    std::vector<ScriptObject::Operation>    p_operations;
    bool p_allOperationsAreFixedTime;

    CompiledFunctions                       p_compiledFunctions;
    std::vector<CompiledFunctions>          p_preparedFunctions;

    std::mutex          p_mutex;
    int                 p_dataSize;
    unsigned int        p_simulationCount;
    unsigned int        p_currentSimulation;

    QSqlDatabase        p_database;

    // TODO: turn into pairs
    std::vector<std::vector<std::shared_ptr<SimulationResult>>> p_simulationResults;
    std::vector<std::vector<std::shared_ptr<SimulationDescription>>> p_simulationDescriptions;

    std::vector<unsigned int> p_simulationIndexes;
};

namespace ComputeFor {
    enum {Initialization = 1, Step = 2, Operation = 3,
          Max};
}

/* Inline functions *****************************************************************************************/

const ScriptObject::Settings & SimulationContext::settings() const
{
    return p_settings;
}

const ScriptObject::Simulation & SimulationContext::simulation() const
{
    return p_simulation;
}

bool SimulationContext::allOperationsAreFixedTime() const
{
    return p_allOperationsAreFixedTime;
}

#endif // SIMULATION_H
