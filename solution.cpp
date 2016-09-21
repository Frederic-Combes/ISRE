#include "solution.h"

#include "result.h"

#include <QDebug>

using namespace std;
using namespace Eigen;
using namespace ScriptObject;

Solution::Solution(const Settings & settings)
{
    unsigned int numberOfEnergySamples   = settings.energySamples().size();
    unsigned int numberOfEnergyBins      = settings.energyBins().size();

    current.spin.resize(numberOfEnergySamples);
    current.partial.resize(numberOfEnergyBins);

    previous.partial.resize(numberOfEnergyBins);
}

Solution::~Solution()
{
}

void Solution::save(double time, const std::vector<ResultPtr> &results)
{
    if(results.size() > 1) {
        for(unsigned int i = 0; i < (unsigned int)(results.size()-1); ++i) {
            lock_guard<mutex> lock(results.at(i)->lock()); (void) lock;
            vector<double> & values = results[i]->at(time);

            values[Result::Time]            = time;
            values[Result::SpinX]           = current.partial(i).spin(0);
            values[Result::SpinY]           = current.partial(i).spin(1);
            values[Result::SpinZ]           = current.partial(i).spin(2);
            values[Result::Contrast]        = current.partial(i).spin.head<3>().norm();
            values[Result::AtomLossFromF1]  = 0;
            values[Result::AtomLossFromF2]  = 0;
            values[Result::Phase]           = current.partial(i).phase;
        }
    }

    lock_guard<mutex> lock(results.back()->lock()); (void) lock;
    vector<double> & values = results.back()->at(time);

    values[Result::Time]            = time;
    values[Result::SpinX]           = current.total.spin(0);
    values[Result::SpinY]           = current.total.spin(1);
    values[Result::SpinZ]           = current.total.spin(2);
    values[Result::Contrast]        = current.total.spin.head<3>().norm();
    values[Result::AtomLossFromF1]  = current.total.atomLosses(1);
    values[Result::AtomLossFromF2]  = current.total.atomLosses(2);
    values[Result::Phase]           = current.total.phase;
}

void Solution::clear()
{
    this->timeStep = 0;

    current.clear();
    previous.clear();
    //save.clear();
}

void Solution::Current::clear()
{
    fill(spin.begin(), spin.end(), Vector4d::Zero().eval());

    total.clear();
    partial.clear();
}

void Solution::Current::Total::clear()
{
    spin  = Vector4d::Zero();
    phase = 0;

    atomLosses.clear();
}

void Solution::Current::Total::AtomLosses::clear()
{
    p_fromF1 = 0;
    p_fromF2 = 0;
}

void Solution::Current::Partial::clear()
{
    for_each(p_internal.begin(), p_internal.end(), [](Internal & internal){internal.clear();});
}

void Solution::Current::Partial::Internal::clear()
{
    spin  = Vector4d::Zero();
    phase = 0;
}

void Solution::Previous::clear()
{
    total.clear();
    partial.clear();
}

void Solution::Previous::Total::clear()
{
    spin  = Vector4d::Zero();
    phase = 0;
}

void Solution::Previous::Partial::clear()
{
    for_each(p_internal.begin(), p_internal.end(), [](Internal & internal){internal.clear();});
}

void Solution::Previous::Partial::Internal::clear()
{
    spin  = Vector4d::Zero();
    phase = 0;
}

/*void Solution::Save::clear()
{
    p_internal.clear();
}

unsigned int Solution::Save::Internal::Partial::size() const
{
    return p_internal.size();
}

void Solution::Save::Internal::clear()
{
    timestep = 0;
    total.clear();
    partial.clear();
}

void Solution::Save::Internal::Total::clear()
{
    spin  = Vector4d::Zero();
    phase = 0;

    atomLosses.clear();
}

void Solution::Save::Internal::Total::AtomLosses::clear()
{
    p_fromF1 = 0;
    p_fromF2 = 0;
}

void Solution::Save::Internal::Partial::clear()
{
    for_each(p_internal.begin(), p_internal.end(), [](Internal & internal) {internal.clear();});
}

void Solution::Save::Internal::Partial::Internal::clear()
{
    spin  = Vector4d::Zero();
    phase = 0;
}*/
