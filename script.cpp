#include "script.h"

#include <QFile>
#include <QTextStream>

#include <QScriptEngine>
#include <QScriptValue>

#include <QSqlQuery>
#include <QSqlError>

#include <QDebug>

#include <boost/math/special_functions/bessel.hpp>

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#pragma GCC diagnostic pop

using namespace std;

// Computes the bessel function of the second kind I_0(arg)
QScriptValue besselI_0(QScriptContext *context, QScriptEngine *engine)
{
    if(context->argumentCount() != 1)
        return QScriptValue(QScriptValue::UndefinedValue);

    if(!context->argument(0).isNumber())
        return QScriptValue(QScriptValue::UndefinedValue);

    return QScriptValue(engine, boost::math::cyl_bessel_i(0.0, context->argument(0).toNumber()));
}

QScriptValue loadScript(QScriptContext *context, QScriptEngine * engine)
{
    if(context->argumentCount() != 1) {
        qDebug() << "No argument given to loadScript()";
        return QScriptValue::UndefinedValue;
    }

    if(!context->argument(0).isString()) {
        qDebug() << "Invalid argument given to loadScript()";
        return QScriptValue::UndefinedValue;
    }

    QFile scriptFile(context->argument(0).toString());
    if(!scriptFile.open(QIODevice::ReadOnly)) {
        qDebug() << "Failed to open the file" << context->argument(0).toString();
        return QScriptValue::UndefinedValue;
    }

    QTextStream scriptStream(&scriptFile);
    QString script = scriptStream.readAll();

    QScriptSyntaxCheckResult syntaxCheck = engine->checkSyntax(script);

    if(syntaxCheck.state() != QScriptSyntaxCheckResult::Valid) {
        qDebug() << "Syntax check for script" << context->argument(0).toString() << "failed";
        qDebug() << "Line" << syntaxCheck.errorLineNumber() << ": error :" << syntaxCheck.errorMessage();
        return QScriptValue::UndefinedValue;
    }

    engine->evaluate(script);

    if(engine->hasUncaughtException())
    {
        qDebug() << "Uncaught exception within the script" << context->argument(0).toString();
        qDebug() << "Line" << engine->uncaughtExceptionLineNumber() << ": error :" << engine->uncaughtException().toString();
        qDebug() << "Backtrace :" << engine->uncaughtExceptionBacktrace();
    }

    return QScriptValue::UndefinedValue;
}

Eigen::Matrix3d toMatrix3d(const QScriptValue & value)
{
    Eigen::Matrix3d matrix;

    matrix << value.property(0).property(0).toNumber(), value.property(0).property(1).toNumber(), value.property(0).property(2).toNumber(),
              value.property(1).property(0).toNumber(), value.property(1).property(1).toNumber(), value.property(1).property(2).toNumber(),
              value.property(2).property(0).toNumber(), value.property(2).property(1).toNumber(), value.property(2).property(2).toNumber();

    qDebug() << "Matrix :";
    qDebug() << matrix(0,0) << matrix(0,1) << matrix(0,2) << endl
             << matrix(1,0) << matrix(1,1) << matrix(1,2) << endl
             << matrix(2,0) << matrix(2,1) << matrix(2,2);

    return matrix;
}

bool processFile(QString filename, SimulationTools &tools)
{
    // Opening "setup.js" to load default values

    QFile setupFile(":/Setup/setup.js");

    if(!setupFile.open(QIODevice::ReadOnly))
    {
        qWarning() << "The setup simulation file ""setup.js"" could not be opened.";
        return false;
    }

    // Opening the file defining this simulation parameters

    QFile simulationFile(filename);

    if(!simulationFile.open(QIODevice::ReadOnly))
    {
        qWarning() << "The file" << filename << "could not be open.";
        return false;
    }

    // Processing both files

    QScriptEngine engine;

    // Provide acces to specific functions

    // Script loading
    QScriptValue loadScriptFunction = engine.newFunction(loadScript, 1);
    engine.globalObject().setProperty("loadScript", loadScriptFunction);

    // Bessel I_0 function
    QScriptValue besselI_0Function = engine.newFunction(besselI_0, 1);
    engine.globalObject().setProperty("besselI_0", besselI_0Function);

    // Syntax checking for the setup file

    QTextStream setupStream(&setupFile);
    QString setup = setupStream.readAll();
    QScriptSyntaxCheckResult setupSyntaxCheck = engine.checkSyntax(setup);

    if(setupSyntaxCheck.state() != QScriptSyntaxCheckResult::Valid)
    {
        qDebug() << "Line" << setupSyntaxCheck.errorLineNumber() << ": error:" << setupSyntaxCheck.errorMessage();
        return false;
    }

    // Processing the setup file

    engine.evaluate(setup);

    if(engine.hasUncaughtException())
    {
        qDebug() << "Uncaught exception within the file ""setup.js""";
        qDebug() << engine.uncaughtException().toString();
        qDebug() << "Line" << engine.uncaughtExceptionLineNumber();
        qDebug() << engine.uncaughtExceptionBacktrace();
        return false;
    }

    // Syntax checking for the simulation file

    QTextStream simulationStream(&simulationFile);
    QString simulation = simulationStream.readAll();
    QScriptSyntaxCheckResult simulationSyntaxCheck = engine.checkSyntax(simulation);

    if(simulationSyntaxCheck.state() != QScriptSyntaxCheckResult::Valid) {
        qDebug() << "Line" << simulationSyntaxCheck.errorLineNumber() << ": Syntax error:" << simulationSyntaxCheck.errorMessage();
        return false;
    }

    // Processing

    engine.evaluate(simulation);

    if(engine.hasUncaughtException()) {
        qDebug() << "Uncaught exception within the file" << filename;
        qDebug() << engine.uncaughtException().toString();
        qDebug() << "Line" << engine.uncaughtExceptionLineNumber();
        qDebug() << engine.uncaughtExceptionBacktrace();
        return false;
    }

    QScriptValue simulationObject = engine.evaluate("simulation");

    // Check that all values are valid and turns single values to arrays
    if(!simulationObject.property("check").call(simulationObject).toBool()) {
        return false;
    }

    // Space dimension and energy sampling
    tools.sampleEnergy(simulationObject.property("dimension").toInteger(), simulationObject.property("energy").property("count").toInteger());

    // Time settings
    {
        QScriptValue timeSettings = simulationObject.property("time");
        tools.setTotalTime(timeSettings.property("length").toNumber());
        tools.setTimeStepsCount(timeSettings.property("steps").toInteger());
        tools.setTimePointsCount(timeSettings.property("points").toInteger());
    }

    // Reads the energy bins
    {
        QScriptValue energyBinSettings = simulationObject.property("energy").property("bins");
        int length = energyBinSettings.property("length").toNumber();
        for(int i = 0; i < length; ++i) {
            tools.addEnergyBin(energyBinSettings.property(i).property(0).toNumber(),
                               energyBinSettings.property(i).property(1).toNumber());
        }
    }

    // Adding operations (pi-pulses, ...)
    {
        qDebug() << "Adding operations...";
        QScriptValue operations = engine.evaluate("simulation.operations");
        int length = operations.property("length").toInteger();
        for(int i = 0; i < length; ++i) {
            QScriptValue scriptOperation = operations.property(i);

            Operation operation;
            operation.isFixedTime       = scriptOperation.property("isFixedTime").toNumber();
            operation.time              = scriptOperation.property("time").toNumber();
            operation.operationMatrix   = toMatrix3d(scriptOperation.property("operator"));

            if(abs(abs(operation.operationMatrix.determinant()) - 1.0) > 1e-6) {
                qDebug() << "The operation determinant is neither 1 nor -1! (" << operation.operationMatrix.determinant() << ")";
            }

            tools.addOperation(operation);
        }
    }

    // Time-dependence of the density
    vector<double> & densityTimeDependence = tools.densityTimeDependance();
    densityTimeDependence.resize(tools.timeStepsCount()+1);

    for(int i = 0; i <= tools.timeStepsCount(); ++i) {
        double timeInSeconds = double(i) / double(tools.timeStepsCount()) * tools.totalTime();
        densityTimeDependence[i] = engine.evaluate(QString("simulation.density.timeDependence(")+QString::number(timeInSeconds)+QString(")")).toNumber();
    }

    // Energy dependence of the density
    vector<double> & densityEnergyDependence = tools.densityEnergyDependance();
    densityEnergyDependence.resize(tools.energyCount());
    for(int i = 0; i < tools.energyCount(); ++i) {
        densityEnergyDependence[i] = engine.evaluate(QString("simulation.density.energyDependence(") + QString::number(tools.energy(i)) + QString(")")).toNumber();
    }

    // Simulation name
    tools.setSimulationName(simulationObject.property("name").toString().left(255));

    // Create the data queue
    vector<double> fl = makeList(&engine, "frequency.larmor");
    vector<double> il = makeList(&engine, "inhomogeneity.larmor");
    vector<double> id = makeList(&engine, "inhomogeneity.density");
    vector<double> fe = makeList(&engine, "frequency.exchange");
    vector<double> r  = makeList(&engine, "relaxation");
    vector<double> d  = makeList(&engine, "density.value");
    vector<double> l1 = makeList(&engine, "atomLoss.fromF1");
    vector<double> l2 = makeList(&engine, "atomLoss.fromF2");

    QSqlQuery query(tools.database());

    query.exec("BEGIN");

    for(double larmorFrequency : fl) {
    for(double larmorInhomogeneity : il) {
    for(double densityInhomogeneity : id) {
    for(double exchangeFrequency: fe) {
    for(double damping : r) {
    for(double lossFromF1 : l1) {
    for(double lossFromF2 : l2) {
    for(double density : d) {
        DataPoint point;

        point.larmor                = larmorFrequency;
        point.larmorInhomogeneity   = larmorInhomogeneity;
        point.densityInhomogeneity  = densityInhomogeneity * density;
        point.exchange              = exchangeFrequency * density;
        point.damping               = damping * density;
        point.lossFromF1            = lossFromF1 * density;
        point.lossFromF2            = lossFromF2 * density;

        // Inserts the simuation into the database
        query.prepare("INSERT INTO simulation"
                      " (spin_echo, larmor_frequency, larmor_inhomogeneity, density_inhomogeneity, exchange_frequency, damping_rate, loss_1, loss_2, name)"
                      " VALUES (:s, :l, :i, :d, :e, :r, :l1, :l2, :n);");
        query.bindValue(":s", int(tools.hasMovingOperation()));
        query.bindValue(":l", point.larmor);
        query.bindValue(":i", point.larmorInhomogeneity);
        query.bindValue(":d", point.densityInhomogeneity);
        query.bindValue(":e", point.exchange);
        query.bindValue(":r", point.damping);
        query.bindValue(":l1", point.lossFromF1);
        query.bindValue(":l2", point.lossFromF2);
        query.bindValue(":n", tools.simulationName());

        query.exec();

        // Obtains the simulation id
        point.simulationId = query.lastInsertId().toInt();

        // Inserts the energy bins into the database
        point.energyBinId = vector<int>();
        point.energyBinId.resize(tools.energyBinCount());
        for(int i = 0; i < tools.energyBinCount(); ++i) {
            query.prepare("INSERT INTO energy_bin (simulation_id, energy_min, energy_max) VALUES (:s, :a, :b);");
            query.bindValue(":s", point.simulationId);
            query.bindValue(":a", tools.energyBin(i).first);
            query.bindValue(":b", tools.energyBin(i).second);

            query.exec();

            // Obtains the energy bin id
            point.energyBinId[i] = query.lastInsertId().toInt();
        }

        if(tools.hasMovingOperation()) {
            for(int i = 1; i <= tools.timePointsCount(); ++i) {
                point.totalTime             = tools.totalTime() * ((double)i) / ((double)tools.timePointsCount());
                point.timeStepCount         = tools.timeStepsCount() * ((double)i) / ((double)tools.timePointsCount());

                // Adds the point to the queue
                tools.queue().push(point);
            }
        } else {
            point.totalTime     = tools.totalTime();
            point.timeStepCount = tools.timeStepsCount();

            tools.queue().push(point);
        }
    } // d
    } // l2
    } // l1
    } // r
    } // fe
    } // id
    } // il
    } // fl

    query.exec("COMMIT");

    // Initialize the tools object

    tools.init();

    return true;
}

vector<double> makeList(QScriptEngine * engine, QString pName)
{
    int length = engine->evaluate("simulation."+pName+".length").toNumber();

    vector<double> result;

    for(int i = 0; i < length; ++i) {
        result.push_back(engine->evaluate("simulation." + pName + "[" + QString::number(i) + "]").toNumber());
    }

    return result;
}
