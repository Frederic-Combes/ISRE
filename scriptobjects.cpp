#include "scriptobjects.h"

#include "strings.h"
#include "result.h"

#include <QScriptContext>
#include <QScriptEngine>

#include <QFile>
#include <QTextStream>

#include <boost/math/special_functions/bessel.hpp>

#include <QDebug>

using namespace std;
namespace ScriptObject {

/*************************************************************************************************************
 * Class Operation
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Public Members
 ************************************************************************************************************/

Operation Operation::fromScriptValue(QScriptValue value)
{
    Operation operation;

    operation.time        = value.property(Strings::Script::Operation::Time).toNumber();
    operation.isFixedTime = value.property(Strings::Script::Operation::IsFixedTime).toNumber();

    QScriptValue matrix = value.property(Strings::Script::Operation::Operator);
    operation.matrix <<
            matrix.property(0).property(0).toNumber(), matrix.property(0).property(1).toNumber(), matrix.property(0).property(2).toNumber(),
            matrix.property(1).property(0).toNumber(), matrix.property(1).property(1).toNumber(), matrix.property(1).property(2).toNumber(),
            matrix.property(2).property(0).toNumber(), matrix.property(2).property(1).toNumber(), matrix.property(2).property(2).toNumber();

    return operation;
}

/*************************************************************************************************************
 * Class Parameter
 ************************************************************************************************************/

Parameter Parameter::fromScriptValue(QScriptValue value)
{
    Parameter parameter;

    parameter.name = value.property(Strings::Script::Parameter::Name).toString();

    for(unsigned int i = 0; i < value.property(Strings::Script::Parameter::Values).property("length").toUInt32(); ++i) {
        parameter.values.push_back(value.property(Strings::Script::Parameter::Values).property(i).toNumber());
    }

    parameter.index = 0;

    return parameter;
}

/*************************************************************************************************************
 * Class Function
 ************************************************************************************************************/

Function Function::fromScriptValue(QScriptValue value)
{
    Function function;

    function.name = value.property(Strings::Script::Function::Name).toString();
    function.code = value.property(Strings::Script::Function::Code).toString();

    return function;
}


/*************************************************************************************************************
 * Class Settings
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Public Members
 ************************************************************************************************************/

double Settings::defaultDensityEnergyDependenceDimension1(double energy)
{
    return exp(-0.5*energy) * boost::math::cyl_bessel_i(0, 0.5*energy);
}

double Settings::defaultDensityEnergyDependenceDimension2(double energy)
{
    return exp(-0.5*energy) * boost::math::cyl_bessel_i(0, 0.426119 * energy);
}

double Settings::defaultDensityEnergyDependenceDimension3(double energy)
{
    return exp(-0.5*energy) * boost::math::cyl_bessel_i(0, 0.375126 * energy);
}

Settings Settings::fromScriptValue(QScriptValue value)
{
    Settings settings;

    settings.p_dimension    = value.property(Strings::Script::Settings::Dimension).toInt32();
    settings.p_duration     = value.property(Strings::Script::Settings::Duration).toNumber();
    settings.p_timestepCount= value.property(Strings::Script::Settings::Timesteps).toInt32();
    settings.p_saveInterval = value.property(Strings::Script::Settings::SaveInterval).toNumber();

    settings.sampleEnergy(value.property(Strings::Script::Settings::EnergySamples).toInt32());

    for(unsigned int i = 0; i < value.property(Strings::Script::Settings::EnergyBins).property("length").toUInt32(); ++i) {
        QScriptValue energyBin = value.property(Strings::Script::Settings::EnergyBins).property(i);
        settings.addEnergyBin(energyBin.property(0).toNumber(), energyBin.property(1).toNumber());
    }

    return settings;
}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

void Settings::sampleEnergy(int sampleCount)
{
    p_energySamples.resize(sampleCount);
    fill(p_energySamples.begin(), p_energySamples.end(), 0);

    if(sampleCount <= 0) {
        return;
    }

    switch(p_dimension)
    {
    case 1:
        for(int i = 0; i < sampleCount; ++i) {
            p_energySamples[i] = -log(1.0-double(i)/double(sampleCount));
        }
        break;
    case 2:
        p_energySamples[0] = 0;
        p_energySamples[1] = pow(2.0 / double(sampleCount), 0.5);
        for(int i = 2; i < sampleCount; ++i) {
            p_energySamples[i] = p_energySamples[i-1] + exp(p_energySamples[i-1])/(p_energySamples[i-1]*double(sampleCount));
        }
        break;
    case 3:
        p_energySamples[0] = 0;
        p_energySamples[1] = pow(8.0 / double(sampleCount), 1.0/3.0);
        for(int i = 2; i < sampleCount; ++i) {
            p_energySamples[i] = p_energySamples[i-1] + 2.0 * exp(p_energySamples[i-1]) / (p_energySamples[i-1]*p_energySamples[i-1] * double(sampleCount));
        }
        break;
    default:
        qDebug() << "Supported space dimensions are 1, 2 and 3 - this message should never appear.";
        return;
    }
}

void Settings::addEnergyBin(double energyMin, double energyMax)
{
    // Reorder the energy range
    if(energyMin > energyMax) {
        swap(energyMin, energyMax);
    }

    vector<double>::const_iterator begin = find_if(p_energySamples.begin(), p_energySamples.end(),
    [=](double energy) {
        return energy >= energyMin;
    });

    vector<double>::const_iterator end = find_if(p_energySamples.begin(), p_energySamples.end(),
    [=](double energy) {
        return energy > energyMax;
    });

    if(begin != end) {
        p_energyBins.push_back(make_pair(distance(p_energySamples.cbegin(), begin), distance(p_energySamples.cbegin(), end)));
    } else {
        qDebug() << "Energy bins: bin (" << energyMin << "," << energyMax << ") not added: no energy samples in this range";
    }
}

/*************************************************************************************************************
 * Class Result
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Public Members
 ************************************************************************************************************/

QScriptValue Result::scriptFunctionMatches(QScriptContext * context, QScriptEngine * engine)
{
    // Check argument count
    if(context->argumentCount() != 2) {
        qDebug() << "ScriptObject::Result.matches called with an incorrect parameter count:" << context->argumentCount() << ", expected 2 arguments.";
        return QScriptValue::UndefinedValue;
    }

    // Check argument type
    if(!context->argument(0).isString()) {
        qDebug() << "ScriptObject::Result.matches: argument 0 is not a string";
        return QScriptValue::UndefinedValue;
    }

    if(!context->argument(1).isArray()) {
        qDebug() << "ScriptObject::Result.matches: argument 1 is not an array";
        return QScriptValue::UndefinedValue;
    }

    Result result = context->thisObject().data().toVariant().value<Result>();

    // Check whether the object metadata contains the field given as argument 0
    auto pair = result.metadata(context->argument(0).toString());

    if(!pair.first) {
        return false;
    }

    QString & value = pair.second;

    // Check whether the value matches
    for(qint32 i = 0; i < context->argument(1).property("length").toInt32(); ++i) {
        QScriptValue compareTo = context->argument(1).property(i);

        if(compareTo.isString()) {
            if(QString::compare(value, compareTo.toString()) == 0) {
                return true;
            }
        } else if(compareTo.isNumber()) {
            if(value.toDouble() == compareTo.toNumber()) {
                // TODO : precision / conversion issues
                return true;
            }
        } else if(compareTo.isFunction()){
            if(compareTo.call(QScriptValue(), QScriptValue(engine, value)).toBool()) {
                return true;
            }
        } else {
            qDebug() << "ScriptObject::Result.matches: value is neither a string nor a number nor a function";
        }
    }

    return false;
}

QScriptValue Result::scriptFunctionExport(QScriptContext * context, QScriptEngine * /*engine*/)
{
    // Check argument count
    if(context->argumentCount() != 3) {
        qDebug() << "ScriptObject::Result.export called with an incorrect parameter count:" << context->argumentCount() << ", expected 3 arguments.";
        return QScriptValue::UndefinedValue;
    }

    // Check argument type
    if(!context->argument(0).isString()) {
        qDebug() << "ScriptObject::Result.export: argument 0 is not a string";
        return QScriptValue::UndefinedValue;
    }

    if(!context->argument(1).isString()) {
        qDebug() << "ScriptObject::Result.export: argument 1 is not a string";
        return QScriptValue::UndefinedValue;
    }

    if(!context->argument(2).isArray()) {
        qDebug() << "ScriptObject::Result.export: argument 2 is not an array of string";
        return QScriptValue::UndefinedValue;
    }

    Result result = context->thisObject().data().toVariant().value<Result>();

    QString filename = context->argument(0).toString();
    QString options = context->argument(1).toString();

    // Converts column names to indexes
    vector<SimulationResult::Index> columnIndexes;
    for(quint32 i = 0; i < context->argument(2).property("length").toUInt32(); ++i) {
        SimulationResult::Index columnIndex;

        if(SimulationResult::convertToColumnIndex(context->argument(2).property(i).toString(), columnIndex)) {
            columnIndexes.push_back(columnIndex);
        }
    }

    qDebug() << "ScriptObject::Result.write: opening file" << filename;

    QFile file(filename);

    QIODevice::OpenMode openMode = QIODevice::WriteOnly;
    openMode |= options.contains("a") ? QIODevice::Append : QIODevice::OpenModeFlag(0);
    openMode |= options.contains("o") ? QIODevice::Truncate : QIODevice::OpenModeFlag(0);

    if(!file.open(openMode)) {
        qDebug() << "ScriptObject::Result.export: failed to open the file" << filename;
        return QScriptValue::UndefinedValue;
    }

    QTextStream stream(&file);

    stream.setFieldWidth(12);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    for(double time : *result.p_simulationResult.get()) {
        for(auto columnIndex : columnIndexes) {
            stream << result.p_simulationResult->interpolate(time, columnIndex);
        }

        stream.setFieldWidth(0);
        stream << endl;
        stream.setFieldWidth(12);
    }
    stream.setFieldWidth(0);
    stream << endl;

    file.close();
    return QScriptValue::UndefinedValue;

}

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

Result::Result()
{}

Result::Result(unsigned int id, shared_ptr<SimulationResult> simulationResult, shared_ptr<SimulationDescription> simulationDescription)
    : p_id(id), p_simulationResult(simulationResult), p_simulationDescription(simulationDescription) {}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

void Result::setMetadata(QString name, QString value)
{
    if(p_simulationDescription == nullptr) {
        p_simulationDescription = make_shared<SimulationDescription>();
    }

    p_simulationDescription->setValue(name, value);
}

pair<bool, QString> Result::metadata(QString name) const
{
    if(p_simulationDescription == nullptr) {
        return {false, QString()};
    }

    return p_simulationDescription->value(name);
}

} // namespace ScriptObject
