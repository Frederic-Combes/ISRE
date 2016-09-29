#include "script.h"
#include "simulation.h"
#include "simulationresult.h"
#include "compute.h"
#include "strings.h"

// Multithreading
#include <thread>
#include <mutex>
#include <chrono>

// File support
#include <QFile>
#include <QTextStream>

// Script support
#include <QScriptEngine>
#include <QScriptValue>
#include "scriptobjects.h"

// Database support
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlError>

#include <QDebug>

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#pragma GCC diagnostic pop

using namespace std;
using namespace Eigen;
using namespace ScriptObject;

/*************************************************************************************************************
 * Class CompiledFunction
 ************************************************************************************************************/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

CompiledFunction::CompiledFunction(const QString & name, unsigned int index, const QString & code)
    : p_name(name), p_code(code), p_index(index), p_dependencyBuilt(false), p_valid(true), p_flag(0)
{
     QRegExp regexp(QString("\\b") + p_name + QString("\\b"));

     if(regexp.indexIn(p_code) != -1) {
         p_valid = false;
         qDebug() << "CompiledFunction: Function depends on itself!";
     }
}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

bool CompiledFunction::buildDependency(const vector<Parameter> & parameters, const vector<CompiledFunction *> & functions)
{
    if(p_dependencyBuilt || !p_valid) {
        return p_valid;
    }

    p_dependencyBuilt = true;

    for(const Parameter & parameter : parameters) {
        QRegExp regexp(QString("\\b") + parameter.name + QString("\\b"));

        if(regexp.indexIn(p_code) != -1) {
            p_dependencies.insert(parameter.name);
        }
    }

    for(CompiledFunction * function : functions) {
        QRegExp regexp(QString("\\b") + function->name() + QString("\\b"));

        if(regexp.indexIn(p_code) != -1) {
            function->buildDependency(parameters, functions);
            p_dependencies.insert(function->name());

            for(QString name : function->p_dependencies) {
                p_dependencies.insert(name);
            }
        }
    }

    if(p_dependencies.find(p_name) != p_dependencies.end()) {
        p_valid = false;
        qDebug() << "CompiledFunction::buildDependency : Function depends on itself!";
    }

    return p_valid;
}

bool CompiledFunction::compile(const mathpresso::Context & context)
{
    if(!p_valid) {
        return false;
    }

    mathpresso::Error error = p_compiledCode.compile(context, p_code.toUtf8().data(), mathpresso::kNoOptions);

    if(error != mathpresso::kErrorOk) {
        qDebug() << "CompiledFunction::compile : Error" << error;
        return false;
    }

    return true;
}

bool CompiledFunction::dependsOn(const QString & name) const
{
    if(!p_dependencyBuilt) {
        qDebug() << "CompiledFunction::dependsOn : called before buildDependency()";
        return true;
    }

    return p_dependencies.find(name) != p_dependencies.end();
}

/*************************************************************************************************************
 * Class ScriptLoader
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Members
 ************************************************************************************************************/

QScriptValue ScriptLoader::scriptFunctionLoadScript(QScriptContext * context, QScriptEngine *)
{
    // Check the argument count
    if(context->argumentCount() != 1) {
        qDebug() << "ScriptLoader::scriptFunction:loadScript: expecting exactly one argument.";
        return QScriptValue::UndefinedValue;
    }

    // Check the argument type
    if(!context->argument(0).isString()) {
        qDebug() << "ScriptLoader::scriptFunction:loadScript: argument is not a string";
        return QScriptValue::UndefinedValue;
    }

    // Locate the corresponding file
    QString filename = context->argument(0).toString();

    QFile file(filename);

    if(!file.open(QIODevice::ReadOnly)) {
        qDebug() << "ScriptLoader::scriptFunction:loadScript: failed to open the file" << filename.toUtf8().data();
        return QScriptValue::UndefinedValue;

    }

    // Load the script
    if(!ScriptLoader::load(QTextStream(&file).readAll())) {
        qDebug() << "ScriptLoader::scriptFunction:loadScript: failed to load the script";
    }

    return QScriptValue::UndefinedValue;
}

QScriptValue ScriptLoader::scriptFunctionSimulationConstructor(QScriptContext * context, QScriptEngine * engine)
{
    using namespace Strings::Script;

    QScriptValue thisObject;

    if(context->isCalledAsConstructor()) {
        thisObject = context->thisObject();
    } else {
        thisObject = engine->newObject();
        thisObject.setPrototype(context->callee().prototype());
    }

    // Bind C++ function to Simulation.run
    thisObject.setProperty(Simulation::Run, engine->newFunction(ScriptLoader::scriptFunctionSimulationRun));

    // TODO: move this to setup.js
    //Default values
    thisObject.setProperty(Simulation::Frequency, QScriptValue(engine, "0"));
    thisObject.setProperty(Simulation::Exchange, QScriptValue(engine, "0"));
    thisObject.setProperty(Simulation::Damping, QScriptValue(engine, "0"));

    QScriptValue atomLoss = engine->newObject();
    atomLoss.setProperty(AtomLoss::FromF1, QScriptValue(engine, "0"));
    atomLoss.setProperty(AtomLoss::FromF2, QScriptValue(engine, "0"));
    thisObject.setProperty(Simulation::AtomLoss, atomLoss);

    QScriptValue intialSpin = engine->newArray();
    intialSpin.setProperty(0, QScriptValue(engine, "1"));
    intialSpin.setProperty(1, QScriptValue(engine, "0"));
    intialSpin.setProperty(2, QScriptValue(engine, "0"));
    thisObject.setProperty(Simulation::InitialSpin, intialSpin);

    thisObject.setProperty(Simulation::Operations, engine->newArray());
    thisObject.setProperty(Simulation::Name, engine->undefinedValue());

    return thisObject;
}

QScriptValue ScriptLoader::scriptFunctionSimulationRun(QScriptContext * context, QScriptEngine * engine)
{
    namespace StrScr = Strings::Script;

    // Check the argument count
    if(context->argumentCount() > 0) {
        qDebug() << "ScriptLoader:Simulation.run called with an incorrect parameter count:" << context->argumentCount() << ", expected no argument.";
        return QScriptValue::UndefinedValue;
    }

    // TODO: Check the arguments type

    // TODO: Check the simulation validity

    QScriptValue sim = context->thisObject();

    Simulation simulation;

    simulation.setName(sim.property(StrScr::Simulation::Name).toString());

    simulation.setCode(Index::Frequency,        sim.property(Strings::Script::Simulation::Frequency).toString());
    simulation.setCode(Index::Exchange,         sim.property(Strings::Script::Simulation::Exchange).toString());
    simulation.setCode(Index::Damping,          sim.property(Strings::Script::Simulation::Damping).toString());
    simulation.setCode(Index::AtomLossFromF1,   sim.property(Strings::Script::Simulation::AtomLoss).property(Strings::Script::AtomLoss::FromF1).toString());
    simulation.setCode(Index::AtomLossFromF2,   sim.property(Strings::Script::Simulation::AtomLoss).property(Strings::Script::AtomLoss::FromF2).toString());
    simulation.setCode(Index::InitialSpinX,     sim.property(Strings::Script::Simulation::InitialSpin).property(0).toString());
    simulation.setCode(Index::InitialSpinY,     sim.property(Strings::Script::Simulation::InitialSpin).property(1).toString());
    simulation.setCode(Index::InitialSpinZ,     sim.property(Strings::Script::Simulation::InitialSpin).property(2).toString());

    SimulationContext simulationContext;

    simulationContext.setSettings(ScriptObject::Settings::fromScriptValue(sim.property(StrScr::Simulation::Settings)));
    simulationContext.setSimulation(simulation);

    for(quint32 i = 0; i < sim.property(StrScr::Simulation::Parameters).property("length").toUInt32(); ++i) {
        simulationContext.addParameter(ScriptObject::Parameter::fromScriptValue(sim.property(Strings::Script::Simulation::Parameters).property(i)));
    }

    for(quint32 i = 0; i < sim.property(StrScr::Simulation::Functions).property("length").toUInt32(); ++i) {
        simulationContext.addFunction(ScriptObject::Function::fromScriptValue(sim.property(Strings::Script::Simulation::Functions).property(i)));
    }

    for(quint32 j = 0; j < sim.property(StrScr::Simulation::Operations).property("length").toUInt32(); ++j) {
        simulationContext.addOperation(ScriptObject::Operation::fromScriptValue(sim.property(Strings::Script::Simulation::Operations).property(j)));
    }


    try {
        simulationContext.prepare();
    } catch(...) {
        qDebug() << "ScriptLoader::scriptFunction:Simulation.run: simulation is invalid";
        return engine->undefinedValue();
    }

    qDebug() << "ScriptLoader::scriptFunction:Simulation.run: starting simulation";
    auto startTime = chrono::high_resolution_clock::now();

    unsigned int threadCount = thread::hardware_concurrency()-1;
    thread * threadPool[threadCount == 0 ? 1 : threadCount];

    for(thread * & t : threadPool) {
        t = new thread(threadLoop, ref(simulationContext));
    }

    threadLoop(simulationContext);

    for(thread * & t : threadPool) {
        t->join();
        delete t;
    }

    vector<ScriptObject::Result> results = simulationContext.end();

    auto endTime = chrono::high_resolution_clock::now();
    qDebug() << "ScriptLoader::scriptFunction:Simulation.run: simulation ended succesfully ("
             << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << "ms)";

    QScriptValue returnValue = engine->newArray(results.size());

    for(unsigned int i = 0; i < results.size(); ++i) {
        QScriptValue resultScriptObject = engine->evaluate("new Result()");

        resultScriptObject.setProperty(Strings::Script::Result::Matches, engine->newFunction(ScriptObject::Result::scriptFunctionMatches));
        resultScriptObject.setProperty(Strings::Script::Result::Export, engine->newFunction(ScriptObject::Result::scriptFunctionExport));
        //resultScriptObject.setProperty(Strings::Script::Result::Erase, engine->newFunction(ScriptObject::Result::scriptFunctionErase));

        resultScriptObject.setData(engine->newVariant(QVariant::fromValue(results.at(i))));

        returnValue.setProperty(i, resultScriptObject);
    }

    return returnValue;
}

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

QScriptEngine * ScriptLoader::p_engine = NULL;

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

void ScriptLoader::init()
{
    if(p_engine == nullptr) {
        p_engine = new QScriptEngine();

        // Add additional functions to the script environment
        p_engine->globalObject().setProperty("loadScript",  p_engine->newFunction(ScriptLoader::scriptFunctionLoadScript));
        p_engine->globalObject().setProperty("Simulation",  p_engine->newFunction(ScriptLoader::scriptFunctionSimulationConstructor));

        // Load the setup file
        QFile setupFile(":/Setup/setup.js");
        setupFile.open(QIODevice::ReadOnly);

        QString setupScript = QTextStream(&setupFile).readAll();

        performSyntaxCheck(setupScript);
        evaluate(setupScript);
    }
}

bool ScriptLoader::load(QString script)
{
    if(performSyntaxCheck(script) && evaluate(script)) {
        return true;
    }

    return false;
}

/*************************************************************************************************************
 * Private Members
 ************************************************************************************************************/

bool ScriptLoader::performSyntaxCheck(const QString & script)
{
    QScriptSyntaxCheckResult syntaxCheck = p_engine->checkSyntax(script);

    if(syntaxCheck.state() != QScriptSyntaxCheckResult::Valid)
    {
        qDebug() << "ScriptLoader::load: Syntax error in the script";
        qDebug() << "Line" << syntaxCheck.errorLineNumber() << ": error:" << syntaxCheck.errorMessage();
        return false;
    }

    return true;
}

bool ScriptLoader::evaluate(const QString & script)
{
    p_engine->evaluate(script);

    if(p_engine->hasUncaughtException())
    {
        qDebug() << "ScriptLoader::load: Uncaught exception within the script:";
        qDebug() << "ScriptLoader::load:" << p_engine->uncaughtException().toString().toUtf8().data();
        qDebug() << "ScriptLoader::load: Line" << p_engine->uncaughtExceptionLineNumber();
        for(QString line : p_engine->uncaughtExceptionBacktrace())
            qDebug() << "ScriptLoader::load:" << line.toUtf8().data();
        return false;
    }

    return true;
}
