#include "script.h"

#include <QFile>
#include <QTextStream>

#include <QScriptEngine>
#include <QScriptValue>

#include <QSqlQuery>
#include <QSqlError>

#include <QDebug>

using namespace std;

bool processFile(QString filename, ThreadTools & tools)
{
    /* Opening "setup.js" to load default values */

    QFile setupFile("setup.js");

    if(!setupFile.open(QIODevice::ReadOnly))
    {
        qWarning() << "The setup simulation file ""setup.js"" could not be open.";
        return false;
    }

    /* Opening the file defining this simulation parameters */

    QFile simulationFile(filename);

    if(!simulationFile.open(QIODevice::ReadOnly))
    {
        qWarning() << "The file" << filename << "could not be open.";
        return false;
    }

    /* Processing both files */

    QScriptEngine engine;

    QTextStream setupStream(&setupFile);
    QString setup = setupStream.readAll();
    QScriptSyntaxCheckResult setupSyntaxCheck = engine.checkSyntax(setup);

    if(setupSyntaxCheck.state() != QScriptSyntaxCheckResult::Valid)
    {
        qDebug() << "Line" << setupSyntaxCheck.errorLineNumber() << ": Syntax error:" << setupSyntaxCheck.errorMessage();
        return false;
    }

    engine.evaluate(setup);

    if(engine.hasUncaughtException())
    {
        qDebug() << "Uncaught exception within the file ""setup.js""";
        qDebug() << engine.uncaughtException().toString();
        qDebug() << "Line" << engine.uncaughtExceptionLineNumber();
        qDebug() << engine.uncaughtExceptionBacktrace();
        return false;
    }

    QTextStream simulationStream(&simulationFile);
    QString simulation = simulationStream.readAll();
    QScriptSyntaxCheckResult simulationSyntaxCheck = engine.checkSyntax(simulation);

    if(simulationSyntaxCheck.state() != QScriptSyntaxCheckResult::Valid)
    {
        qDebug() << "Line" << simulationSyntaxCheck.errorLineNumber() << ": Syntax error:" << simulationSyntaxCheck.errorMessage();
        return false;
    }

    engine.evaluate(simulation);

    if(engine.hasUncaughtException())
    {
        qDebug() << "Uncaught exception within the file" << filename;
        qDebug() << engine.uncaughtException().toString();
        qDebug() << "Line" << engine.uncaughtExceptionLineNumber();
        qDebug() << engine.uncaughtExceptionBacktrace();
        return false;
    }

    /* Check that all value are valid */
    if(!engine.evaluate("simulation.check()").toBool())
    {
        return false;
    }

    /* Get the values */

    tools.setSpaceDimension(engine.evaluate("simulation.dimension").toInteger());
    tools.setTotalTime(engine.evaluate("simulation.time.length").toNumber());
    tools.setTimeStepsCount(engine.evaluate("simulation.time.steps").toInteger());
    tools.setTimePointsCount(engine.evaluate("simulation.time.points").toInteger());
    tools.setEnergyCount(engine.evaluate("simulation.energy.count").toInteger());

    {
        QString pName = "energy.bins";
        int length = engine.evaluate("simulation."+pName+".length").toNumber();
        for(int i = 0; i < length; i++)
        {
            QString code1 = "simulation." + pName + "[" + QString::number(i) + "][0]";
            QString code2 = "simulation." + pName + "[" + QString::number(i) + "][1]";

            double e1 = engine.evaluate(code1).toNumber();
            double e2 = engine.evaluate(code2).toNumber();

            tools.addEnergyBin(e1,e2);
        }
    }

    tools.enableSpinEcho(engine.evaluate("simulation.spinEcho.enabled").toBool());

    if(tools.spinEcho())
    {
        tools.setSpinEchoTime(makeList(&engine, QString("spinEcho.time")));
    }

    double * densityFactor = new double[tools.timeStepsCount()];
    for(int i = 0; i < tools.timeStepsCount(); i++)
    {
        double timeInSeconds = double(i) / double(tools.timeStepsCount()) * tools.totalTime();
        densityFactor[i] = engine.evaluate("simulation.density.timeDependence("+QString::number(timeInSeconds)+")").toNumber();
    }

    tools.setDensityFactor(densityFactor);

    tools.setSimulationName(engine.evaluate("simulation.name").toString().left(255));

    /* Create the data queue */

    vector<double> fl = makeList(&engine, QString("frequency.larmor"));
    vector<double> fe = makeList(&engine, QString("frequency.exchange"));
    vector<double> il = makeList(&engine, QString("inhomogeneity.larmor"));
    vector<double> id = makeList(&engine, QString("inhomogeneity.density"));
    vector<double> r  = makeList(&engine, QString("relaxation"));
    vector<double> d  = makeList(&engine, QString("density.value"));

    QSqlQuery query(tools.database());

    query.exec("BEGIN");

    if(tools.spinEcho())
    {
        for(auto it1 = fl.begin(); it1 != fl.end(); it1++)
        {
            for(auto it2 = il.begin(); it2 != il.end(); it2++)
            {
                for(auto it3 = id.begin(); it3 != id.end(); it3++)
                {
                    for(auto it4 = fe.begin(); it4 != fe.end(); it4++)
                    {
                        for(auto it5 = r.begin(); it5 != r.end(); it5++)
                        {
                            for(auto it6 = d.begin(); it6 != d.end(); it6++)
                            {
                                DataPoint point;

                                point.larmor                = *it1;
                                point.larmorInhomogeneity   = *it2;
                                point.densityInhomogeneity  = *it3 * *it6;
                                point.exchange              = *it4 * *it6;
                                point.damping               = *it5 * *it6;

                                /* Inserts the simuation into the database */
                                query.prepare("INSERT INTO simulation"
                                              " (spin_echo, larmor_frequency, larmor_inhomogeneity, density_inhomogeneity, exchange_frequency, damping_rate, name)"
                                              " VALUES (1, :l, :i, :d, :e, :r, :n);");
                                query.bindValue(":l", point.larmor);
                                query.bindValue(":i", point.larmorInhomogeneity);
                                query.bindValue(":d", point.densityInhomogeneity);
                                query.bindValue(":e", point.exchange);
                                query.bindValue(":r", point.damping);
                                query.bindValue(":n", tools.simulationName());

                                query.exec();

                                /* Obtains the simulation id */
                                point.simulationId = query.lastInsertId().toInt();

                                /* Inserts the energy bins into the database */
                                point.energyBinId = vector<int>();
                                point.energyBinId.resize(tools.energyBinCount());
                                for(int i = 0; i < tools.energyBinCount(); i++)
                                {
                                    query.prepare("INSERT INTO energy_bin (simulation_id, energy_min, energy_max) VALUES (:s, :a, :b);");
                                    query.bindValue(":s", point.simulationId);
                                    query.bindValue(":a", tools.energyBin(i).first);
                                    query.bindValue(":b", tools.energyBin(i).second);

                                    query.exec();

                                    /* Obtains the energy bin id */
                                    point.energyBinId[i] = query.lastInsertId().toInt();
                                }

                                for(int i = 1; i <= tools.timePointsCount(); i++)
                                {
                                    point.totalTime             = tools.totalTime() * ((double)i) / ((double)tools.timePointsCount());
                                    point.timeStepCount         = tools.timeStepsCount() * ((double)i) / ((double)tools.timePointsCount());

                                    /* Adds the point to the queue */
                                    tools.queue().push(point);
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    else
    {
        for(auto it1 = fl.begin(); it1 != fl.end(); it1++)
        {
            for(auto it2 = il.begin(); it2 != il.end(); it2++)
            {
                for(auto it3 = id.begin(); it3 != id.end(); it3++)
                {
                    for(auto it4 = fe.begin(); it4 != fe.end(); it4++)
                    {
                        for(auto it5 = r.begin(); it5 != r.end(); it5++)
                        {
                            for(auto it6 = d.begin(); it6 != d.end(); it6++)
                            {
                                DataPoint point;

                                point.totalTime             = tools.totalTime();
                                point.timeStepCount         = tools.timeStepsCount();

                                point.larmor                = *it1;
                                point.larmorInhomogeneity   = *it2;
                                point.densityInhomogeneity  = *it3 * *it6;
                                point.exchange              = *it4 * *it6;
                                point.damping               = *it5 * *it6;

                                /* Inserts the simuation into the database */
                                query.prepare("INSERT INTO simulation"
                                              " (spin_echo, larmor_frequency, larmor_inhomogeneity, density_inhomogeneity, exchange_frequency, damping_rate, name)"
                                              " VALUES (0, :l, :i, :d, :e, :r, :n);");
                                query.bindValue(":l", point.larmor);
                                query.bindValue(":i", point.larmorInhomogeneity);
                                query.bindValue(":d", point.densityInhomogeneity);
                                query.bindValue(":e", point.exchange);
                                query.bindValue(":r", point.damping);
                                query.bindValue(":n", tools.simulationName());

                                query.exec();

                                /* Obtains the simulation id */
                                point.simulationId = query.lastInsertId().toInt();

                                /* Inserts the energy bins into the database */
                                point.energyBinId = vector<int>();
                                point.energyBinId.resize(tools.energyBinCount());
                                for(int i = 0; i < tools.energyBinCount(); i++)
                                {
                                    query.prepare("INSERT INTO energy_bin (simulation_id, energy_min, energy_max) VALUES (:s, :a, :b);");
                                    query.bindValue(":s", point.simulationId);
                                    query.bindValue(":a", tools.energyBin(i).first);
                                    query.bindValue(":b", tools.energyBin(i).second);

                                    query.exec();

                                    /* Obtains the energy bin id */
                                    point.energyBinId.at(i) = query.lastInsertId().toInt();
                                }

                                /* Adds the point to the queue */
                                tools.queue().push(point);
                            }
                        }
                    }
                }
            }
        }
    }

    query.exec("COMMIT");

    /* Initialize the tools object */

    tools.init();

    return true;
}

vector<double> makeList(QScriptEngine * engine, QString pName)
{
    int length = engine->evaluate("simulation."+pName+".length").toNumber();

    vector<double> result;

    for(int i = 0; i < length; i++)
    {
        result.push_back(engine->evaluate("simulation." + pName + "[" + QString::number(i) + "]").toNumber());
    }

    return result;
}


