#include "scriptobjects.h"
#include "strings.h"

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

}
