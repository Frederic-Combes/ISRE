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

// TODO : turn to ScriptObject
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

class SimulationContext
{
public:
    typedef std::vector<std::reference_wrapper<const ScriptObject::Operation>> Operations;
public:
    SimulationContext();
    ~SimulationContext();

    void setOverwrite(bool overwrite = false);

    void setSettings(ScriptObject::Settings settings);
    void setSimulation(const Simulation & simulation);
    void addParameter(ScriptObject::Parameter parameter);
    void addFunction(ScriptObject::Function function);
    void addOperation(ScriptObject::Operation operation);

    inline const ScriptObject::Settings & settings() const {return p_settings;}
    inline const Simulation & simulation() const {return p_simulation;}
    inline bool allOperationsAreFixedTime() const {return p_allOperationsAreFixedTime;}

    const CompiledFunctions & functions(unsigned int index = 0) const;

    void evaluateFunctions(unsigned int index, int flag, double * data) const;

    int getIndex(double * & data);

    // TODO: these 4 functions aren't needed?
    std::mutex &                databaseLock() {return p_dbLock;}
    QSqlDatabase &              database() {return p_database;}
    const std::vector<int> &    ids(int simulationIndex) {return p_indexToIds[simulationIndex];}
    const std::vector<int> &    ids() const {return p_ids;}

    void prepare();
    void end();

    // TODO : auto release / regain the data ptr (as well as the Operations vector)!
    // HINT: use unique_ptr<double[], allocator> ?
    bool next(unsigned int & simulationIndex);          // Set simulationIndex to the next simulation index, and return true. Return false if there are no more simulations to run
    double * data(int simulationIndex) const;           // Returns the data needed for the simulation indexed by simulationIndex
    Operations operations(int simulationIndex) const;   // Returns the time-ordered list of operations needed for the simulation indexed by simulationIndex
    const std::vector<ResultPtr> & results(unsigned int simulationIndex) const;

private:
    void prepareFunctions(mathpresso::Context & context);
    void prepareResults();

    CompiledFunctions functionsRequiredToCompute(std::vector<QString> functionNames) const;
    CompiledFunctions functionsRequiredToCompute(const QString & functionName) const;

private:
    bool                                    p_overwrite;

    ScriptObject::Settings                  p_settings;
    Simulation                              p_simulation;
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

    std::mutex          p_dbLock;
    QSqlDatabase        p_database;

    std::vector<std::vector<ResultPtr>> p_simulationResults;
    std::vector<unsigned int> p_simulationIndexes;
    std::map<unsigned int, std::vector<int> > p_indexToIds;
    std::vector<int> p_ids;
};

namespace ComputeFor {
    enum {Initialization = 1, Step = 2, Operation = 3,
          Max};
}

#endif // SIMULATION_H
