#include "simulation.h"
#include "strings.h"
#include "result.h"

#include <QString>
#include <QDate>

#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlError>

#include <QDebug>

#include <iostream>
#include <cmath>

#include "scriptobjects.h"

using namespace std;
using namespace ScriptObject;

namespace ScriptObject {

/*************************************************************************************************************
 * Class Simulation
 ************************************************************************************************************/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

Simulation::Simulation() {}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

void Simulation::setCode(int index, QString code)
{
    if((unsigned int)(index) >= p_code.size()) {
        p_code.resize(index+1);
    }

    p_code[index] = code;
}

QString Simulation::code(int index) const
{
    if((unsigned int)(index) >= p_code.size()) {
        return QString();
    }

    return p_code[index];
}

} // namespace ScriptObject

/*************************************************************************************************************
 * Class SimulationContext
 ************************************************************************************************************/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

SimulationContext::SimulationContext()
    : p_overwrite(false), p_allOperationsAreFixedTime(true), p_currentSimulation(0)
{
    // We insert the idendity operation at the end of the simulation
    Operation operation;

    operation.isFixedTime = false;
    operation.time = 1;
    operation.matrix = Eigen::Matrix3d::Identity();

    p_operations.push_back(operation);
}

SimulationContext::~SimulationContext()
{
    for(CompiledFunction * function : p_compiledFunctions) {
        delete function;
    }
}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

void SimulationContext::setOverwrite(bool overwrite)
{
    p_overwrite = overwrite;
}

void SimulationContext::setSettings(Settings settings)
{
    p_settings = settings;
}

void SimulationContext::setSimulation(const Simulation & simulation)
{
    p_simulation = simulation;
}

void SimulationContext::addParameter(Parameter parameter)
{
    p_parameters.push_back(parameter);
}

void SimulationContext::addFunction(Function function)
{
    p_functions.push_back(function);
}

void SimulationContext::addOperation(Operation operation)
{
    p_operations.push_back(operation);

    if(operation.isFixedTime == false) {
        p_allOperationsAreFixedTime = false;
    }
}

const CompiledFunctions & SimulationContext::functions(unsigned int index) const
{
    if(index == 0 || index > ComputeFor::Max-1) {
        return p_compiledFunctions;
    }
    else {
        return p_preparedFunctions[index];
    }
}

void SimulationContext::evaluateFunctions(unsigned int index, int flag, double * data) const
{
    for(const CompiledFunction * function : functions(index)) {
        if(function->flag() == flag) {
            data[function->index()] = function->eval(data);
        }
    }
}

int SimulationContext::getIndex(double * & data)
{
    lock_guard<mutex> lock(p_mutex); (void) lock;

    // Allocate data if the pointer is not valid
    if(data == NULL) {
        data = new double[p_dataSize];
    }

    unsigned int index = p_currentSimulation;

    for(Parameter parameter : p_parameters) {
        unsigned int i = index % parameter.values.size();
        index = (index - i) / parameter.values.size();

        data[parameter.index] = parameter.values[i];
    }

    if(index != 0 || p_parameters.size() == 0) {
        return -1;
    }

    data[Index::TimestepCount] = data[Index::Duration] / p_settings.duration() * p_settings.timestepCount();

    ++p_currentSimulation;

    return (p_currentSimulation-1) / p_parameters.front().values.size();
}

void SimulationContext::prepare()
{
    // FIXME: Issue when there are energy bins

    using namespace Strings::JIT;

    // Obtain & open the database
    p_database = QSqlDatabase::database("isre-result");

    if(!p_database.isValid()) {
        p_database = QSqlDatabase::addDatabase("QSQLITE", "isre-result");
        p_database.setDatabaseName("results.db");

    }

    if(!p_database.isOpen() && !p_database.open()) {
        qDebug() << "SimulationContext::prepare: failed to open database";
        throw;
    }

    // Initialize the JIT compiler
    mathpresso::Context context;

    context.addBuiltIns();

    switch(p_settings.dimension()) {
    case 1:
        context.addFunction(Functions::DefaultDensity, (void*)Settings::defaultDensityEnergyDependenceDimension1, mathpresso::kFunctionArg1);
        break;
    case 2:
        context.addFunction(Functions::DefaultDensity, (void*)Settings::defaultDensityEnergyDependenceDimension2, mathpresso::kFunctionArg1);
        break;
    case 3:
        context.addFunction(Functions::DefaultDensity, (void*)Settings::defaultDensityEnergyDependenceDimension3, mathpresso::kFunctionArg1);
        break;
    default:
        break;
    }

    // Create the mandatory functions
    p_compiledFunctions.push_back(new CompiledFunction(Functions::Frequency,        Index::Frequency,       p_simulation.code(Index::Frequency)));
    p_compiledFunctions.push_back(new CompiledFunction(Functions::Exchange,         Index::Exchange,        p_simulation.code(Index::Exchange)));
    p_compiledFunctions.push_back(new CompiledFunction(Functions::Damping,          Index::Damping,         p_simulation.code(Index::Damping)));
    p_compiledFunctions.push_back(new CompiledFunction(Functions::AtomLossFromF1,   Index::AtomLossFromF1,  p_simulation.code(Index::AtomLossFromF1)));
    p_compiledFunctions.push_back(new CompiledFunction(Functions::AtomLossFromF2,   Index::AtomLossFromF2,  p_simulation.code(Index::AtomLossFromF2)));

    p_compiledFunctions.push_back(new CompiledFunction(Functions::InitialSpinX,     Index::InitialSpinX,    p_simulation.code(Index::InitialSpinX)));
    p_compiledFunctions.push_back(new CompiledFunction(Functions::InitialSpinY,     Index::InitialSpinY,    p_simulation.code(Index::InitialSpinY)));
    p_compiledFunctions.push_back(new CompiledFunction(Functions::InitialSpinZ,     Index::InitialSpinZ,    p_simulation.code(Index::InitialSpinZ)));

    // Register the application-defined variables in the JIT
    context.addVariable(Functions::Frequency,       Index::Frequency        * sizeof(double));
    context.addVariable(Functions::Exchange,        Index::Exchange         * sizeof(double));
    context.addVariable(Functions::Damping,         Index::Damping          * sizeof(double));
    context.addVariable(Functions::AtomLossFromF1,  Index::AtomLossFromF1   * sizeof(double));
    context.addVariable(Functions::AtomLossFromF2,  Index::AtomLossFromF2   * sizeof(double));

    context.addVariable(Functions::InitialSpinX,    Index::InitialSpinX     * sizeof(double));
    context.addVariable(Functions::InitialSpinY,    Index::InitialSpinY     * sizeof(double));
    context.addVariable(Functions::InitialSpinZ,    Index::InitialSpinZ     * sizeof(double));

    context.addVariable(Variables::Duration,        Index::Duration         * sizeof(double));
    context.addVariable(Variables::Timestep,        Index::Timestep         * sizeof(double));
    context.addVariable(Variables::EnergyIndex,     Index::EnergyIndex      * sizeof(double));
    context.addVariable(Variables::Time,            Index::Time             * sizeof(double));
    context.addVariable(Variables::Energy,          Index::Energy           * sizeof(double));

    // User-defined variables index in the data[] start from Index::Start
    unsigned int variableIndex = Index::Start;

    // Register the function names in the JIT
    for(const Function & function : p_functions) {
        context.addVariable(function.name.toUtf8().data(), variableIndex * sizeof(double));

        p_compiledFunctions.push_back(new CompiledFunction(function.name, variableIndex, function.code));

        ++variableIndex;
    }

    // Register the parameters in the JIT
    for(Parameter & parameter : p_parameters) {
        context.addVariable(parameter.name.toUtf8().data(), variableIndex * sizeof(double));

        parameter.index = variableIndex;

        ++variableIndex;
    }

    // Number of registered variables/functions/paramter names in the JIT
    // The indexes start at 0, so we have to add 1 to get the actual size
    p_dataSize = variableIndex + 1;

    // Adding a variable duration (in case the simulation has moving operations)
    vector<double> durations;
    if(!p_allOperationsAreFixedTime) {
        for(double time = p_settings.saveInterval(); time < p_settings.duration(); time += p_settings.saveInterval()) {
            durations.push_back(time);
        }
    }
    durations.push_back(p_settings.duration());

    Parameter duration;
    duration.name   = Variables::Duration;
    duration.values = durations;
    duration.index  = Index::Duration;

    // Make sure the duration parameter is the front() of the p_parameters container
    p_parameters.push_back(duration);
    swap(p_parameters.front(), p_parameters.back());

    // Compile & cache functions
    prepareFunctions(context);

    // Saves the simulation descriptions in the database and stores the associated indexes
    prepareResults();
}

vector<ScriptObject::Result> SimulationContext::end()
{
    vector<ScriptObject::Result> savedResults;

    p_database.transaction();

    QSqlQuery query(p_database);
    query.prepare("INSERT INTO simulation (name, timestamp) VALUES (:n, :t);");

    for(unsigned int i = p_parameters.front().values.size() - 1; i < p_simulationCount; i += p_parameters.front().values.size()) {
        vector<shared_ptr<SimulationResult>> & results = p_simulationResults[i / p_parameters.front().values.size()];
        vector<shared_ptr<SimulationDescription>> & descriptions = p_simulationDescriptions[i / p_parameters.front().values.size()];

        query.bindValue(":n", p_simulation.name());
        query.bindValue(":t", QDateTime::currentDateTime().toString("yyyy-MM-dd HH:mm:ss:zzz"));

        if(!query.exec()) {
            qDebug() << "SimulationContext::end: failed to exec query" << query.lastQuery();
            qDebug() << "SimulationContext::end:" << query.lastError().text();
        }

        unsigned int mainSimulationDatabaseId = query.lastInsertId().toUInt();

        results.back()->write(p_database, mainSimulationDatabaseId);
        descriptions.back()->write(p_database, mainSimulationDatabaseId);

        savedResults.push_back(ScriptObject::Result(mainSimulationDatabaseId, results.back(), descriptions.back()));

        // Other shared_prt<SimulationResult> of results points to energy-bin results
        for(unsigned int j = 0; j < results.size() - 1; ++j) {
            descriptions[j]->setValue(Strings::Database::Description::MainSimulationID, QString::number(mainSimulationDatabaseId));

            query.bindValue(":n", p_simulation.name());
            query.bindValue(":t", QDateTime::currentDateTime().toString("yyyy-MM-dd HH:mm:ss:zzz"));

            if(!query.exec()) {
                qDebug() << "SimulationContext::end: failed to exec query" << query.lastQuery();
                qDebug() << "SimulationContext::end:" << query.lastError().text();
            }

            unsigned int databaseId = query.lastInsertId().toUInt();;

            results[j]->write(p_database, databaseId);
            descriptions[j]->write(p_database, databaseId);

            savedResults.push_back(ScriptObject::Result(databaseId, results[j], descriptions[j]));
        }
    }

    p_database.commit();

    return savedResults;
}

bool SimulationContext::next(unsigned int & simulationIndex)
{
    lock_guard<mutex> lock(p_mutex); (void) lock;

    if(p_simulationIndexes.empty()) {
        return false;
    }

    simulationIndex = p_simulationIndexes.back();

    p_simulationIndexes.pop_back();

    return true;
}

double * SimulationContext::data(int simulationIndex) const
{
    double * data = new double[p_dataSize];

    for(const Parameter & parameter : p_parameters) {
        data[parameter.index] = parameter.values.at(simulationIndex % parameter.values.size());
        simulationIndex /= parameter.values.size();
    }

    data[Index::TimestepCount] = data[Index::Duration] / p_settings.duration() * p_settings.timestepCount();
    data[Index::NextSave] = p_allOperationsAreFixedTime ? 0.0 : data[Index::TimestepCount];

    return data;
}

SimulationContext::Operations SimulationContext::operations(int simulationIndex) const
{
    // Duration of the simulation (the front() of p_parameters is the duration parameter)
    double duration = p_parameters.front().values.at(simulationIndex % p_parameters.front().values.size());

    Operations operations;
    operations.reserve(p_operations.size());

    for(const Operation & operation : p_operations) {
        // If the operation takes place before the simulation ends, add it to the list of operations
        if((operation.isFixedTime ? operation.time : operation.time * duration) <= duration) {
                operations.emplace_back(operation);
        }
    }

    // Sort in time-ascending order
    sort(operations.begin(), operations.end(), [=](const Operation & a, const Operation & b) {
        return (a.isFixedTime ? a.time : a.time * duration) < (b.isFixedTime ? b.time : b.time * duration);
    });

    // We are guaranted that the last operation is the indentity taking place at time t = duration (see constructor)
    return operations;
}

const vector<shared_ptr<SimulationResult>> & SimulationContext::results(unsigned int simulationIndex) const
{
    return p_simulationResults.at(simulationIndex / p_parameters.front().values.size());
}

/*************************************************************************************************************
 * Private Members
 ************************************************************************************************************/

void SimulationContext::prepareFunctions(mathpresso::Context & context)
{
    using namespace Strings::JIT;

    // Temporary parameters (for dependency detection)
    p_parameters.push_back(Parameter(Variables::Timestep,   vector<double>()));
    p_parameters.push_back(Parameter(Variables::Time,       vector<double>()));
    p_parameters.push_back(Parameter(Variables::EnergyIndex,vector<double>()));
    p_parameters.push_back(Parameter(Variables::Energy,     vector<double>()));

    // Build function dependency & reorder
    for(CompiledFunction * function : p_compiledFunctions) {
        if(!function->buildDependency(p_parameters, p_compiledFunctions)) {
            qDebug() << "SimulationContext::prepare : Failed to build function dependencies";
            throw;
        }
    }

    sort(p_compiledFunctions.begin(), p_compiledFunctions.end(),
        [](CompiledFunction * a, CompiledFunction * b) {
            if(b->dependsOn(a->name())) {
                return true;
            } else {
                return false;
            }
    });

    // Remove temporary parameters
    p_parameters.pop_back();
    p_parameters.pop_back();
    p_parameters.pop_back();
    p_parameters.pop_back();

    // Compile functions
    for(CompiledFunction * function : p_compiledFunctions) {
        if(!function->compile(context)) {
            qDebug() << "SimulationContext::prepare: Failed to compile user-defined function" << function->name().toUtf8().data();
            qDebug() << "SimulationContext::prepare: (code: " << function->code().toUtf8().data() << ")";
            qDebug() << "SimulationContext::prepare: (simulation name:" << p_simulation.name().toUtf8().data() << ")";
            qDebug() << "SimulationContext::prepare: Did you pass the full list of parameters ?";
            throw;
        }
    }

    // Set function flags (Compute at start, at each step or for each energy)
    for(CompiledFunction * function : p_compiledFunctions) {
        if(function->dependsOn(Variables::Energy) || function->dependsOn(Variables::EnergyIndex))  {
            function->flag() = Compute::OnceForEachEnergy;
        } else {
            if(function->dependsOn(Variables::Time) || function->dependsOn(Variables::Timestep)) {
                function->flag() = Compute::OnceForEachStep;
            } else {
                function->flag() = Compute::AtStart;
            }
        }
    }

    // Looking for caching-friendly functions
    for(CompiledFunction * function : p_compiledFunctions) {
        if(function->flag() == Compute::OnceForEachEnergy && !function->dependsOn(Variables::Time) && !function->dependsOn(Variables::Timestep)) {
            if(any_of(function->dependencies().begin(), function->dependencies().end(), [](const QString & name) {return QString::compare(Variables::Energy, name) != 0 && QString::compare(Variables::EnergyIndex, name) != 0;})) {
                //function->flag() = Compute::CacheAtStart;
            } else {
                //function->flag() = Compute::CacheNow;
            }
        }
    }

    // TODO: Cache functions
    for(CompiledFunction * function : p_compiledFunctions) {
        if(function->flag() == Compute::CacheNow) {
            (void) function;
        }
    }

    // Organize the compiled functions in lists depending on when they need evaluation (side effect: unused functions
    // won't be computed
    p_preparedFunctions.resize(ComputeFor::Max);

    p_preparedFunctions[ComputeFor::Initialization] = functionsRequiredToCompute(vector<QString>({Functions::Frequency, Functions::Exchange, Functions::Damping, Functions::AtomLossFromF2, Functions::AtomLossFromF2, Functions::InitialSpinX, Functions::InitialSpinY, Functions::InitialSpinZ}));
    p_preparedFunctions[ComputeFor::Step]           = functionsRequiredToCompute(vector<QString>({Functions::Frequency, Functions::Exchange, Functions::Damping, Functions::AtomLossFromF2, Functions::AtomLossFromF2}));
    p_preparedFunctions[ComputeFor::Operation]      = functionsRequiredToCompute(vector<QString>({}));
}

CompiledFunctions SimulationContext::functionsRequiredToCompute(std::vector<QString> functionNames) const
{
    CompiledFunctions result;

    for(const QString & name : functionNames) {
         for(CompiledFunction * function : functionsRequiredToCompute(name)) {
             if(find(result.begin(), result.end(), function) == result.end()) {
                 result.push_back(function);
             }
         }
    }

    return result;
}

CompiledFunctions SimulationContext::functionsRequiredToCompute(const QString & functionName) const
{
    auto it_theFunction = find_if(p_compiledFunctions.begin(), p_compiledFunctions.end(),
        [&](CompiledFunction * function) {
            return QString::compare(function->name(), functionName) == 0;
        });

    if(it_theFunction == p_compiledFunctions.end()) {
        qDebug() << "SimulationContext::prepare: failed to find function" << functionName.toUtf8().data();
        return CompiledFunctions();
    }

    CompiledFunction * theFunction = *it_theFunction;
    CompiledFunctions result;

    for(CompiledFunction * function : p_compiledFunctions) {
        if(theFunction->dependsOn(function->name())) {
            result.push_back(function);
        }
    }

    result.push_back(theFunction);

    return result;
}

void SimulationContext::prepareResults()
{
    // Number of simulation runs
    p_simulationCount = accumulate(p_parameters.begin(), p_parameters.end(), int(1),
        [](int value, const Parameter & parameter) {
            return value * parameter.values.size();
        });

    // The first parameter is the duration parameter. Simulation with all paramters identical but the duration
    // parameter share the same result table
    for(unsigned int i = p_parameters.front().values.size() - 1; i < p_simulationCount; i += p_parameters.front().values.size()) {
        p_simulationResults.push_back(vector<shared_ptr<SimulationResult>>());
        p_simulationDescriptions.push_back(vector<shared_ptr<SimulationDescription>>());

        vector<shared_ptr<SimulationResult>> & results = p_simulationResults.back();
        vector<shared_ptr<SimulationDescription>> & descriptions = p_simulationDescriptions.back();

        // Simulation result table for the energy bins
        for(Settings::EnergyBin bin : p_settings.energyBins()) {
            results.push_back(make_shared<SimulationResult>(p_settings.saveCount()));
            descriptions.push_back(make_shared<SimulationDescription>());

            descriptions.back()->setValue(QString(Strings::Database::Description::EnergyBinMinimum), QString::number(p_settings.energySamples().at(bin.first)));
            descriptions.back()->setValue(QString(Strings::Database::Description::EnergyBinMaximum),
                p_settings.energySamples().size() > bin.second ? QString::number(p_settings.energySamples().at(bin.second)) : QString("inf"));
        }

        // Main simulation result table
        results.push_back(make_shared<SimulationResult>(p_settings.saveCount()));
        descriptions.push_back(make_shared<SimulationDescription>());

        // Saving settings in the main result table metadata
        descriptions.back()->setValue(Strings::Database::Description::Dimension,              QString::number(p_settings.dimension()));
        descriptions.back()->setValue(Strings::Database::Description::Duration,               QString::number(p_settings.duration()));
        descriptions.back()->setValue(Strings::Database::Description::NumberOfTimesteps,      QString::number(p_settings.timestepCount()));
        descriptions.back()->setValue(Strings::Database::Description::NumberOfSavedPoints,    QString::number(p_settings.saveCount()));
        descriptions.back()->setValue(Strings::Database::Description::NumberOfEnergySamples,  QString::number(p_settings.energySamples().size()));

        // Saving functions and parameters
        unsigned int j = i;
        for(Parameter parameter : p_parameters) {
            if(QString::compare(parameter.name, Strings::JIT::Variables::Duration) != 0) {
                descriptions.back()->setValue(parameter.name, QString::number(parameter.values[j % parameter.values.size()]));
            }

            j /=  parameter.values.size();
        }

        for(CompiledFunction * function : p_compiledFunctions) {
            descriptions.back()->setValue(function->name(), function->code());
        }

        for(unsigned int j = p_parameters.front().values.size(); j != 0; --j) {
            p_simulationIndexes.push_back(i-j+1);
        }
    }
}
