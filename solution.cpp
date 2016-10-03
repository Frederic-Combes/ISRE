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

void Solution::save(double time, const vector<shared_ptr<SimulationResult>> & results)
{
    if(results.size() > 1) {
        for(unsigned int i = 0; i < (unsigned int)(results.size()-1); ++i) {
            SimulationResult::lock_type lock(results.at(i)->lock()); (void) lock;

            results.at(i)->value(time, SimulationResult::Time)           = time;
            results.at(i)->value(time, SimulationResult::SpinX)          = current.partial(i).spin(0);
            results.at(i)->value(time, SimulationResult::SpinY)          = current.partial(i).spin(1);
            results.at(i)->value(time, SimulationResult::SpinZ)          = current.partial(i).spin(2);
            results.at(i)->value(time, SimulationResult::Contrast)       = current.partial(i).spin.head<3>().norm();
            results.at(i)->value(time, SimulationResult::AtomLossFromF1) = 0;
            results.at(i)->value(time, SimulationResult::AtomLossFromF2) = 0;
            results.at(i)->value(time, SimulationResult::Phase)          = current.partial(i).phase;
        }
    }

    SimulationResult::lock_type lock(results.back()->lock()); (void) lock;

    results.back()->value(time, SimulationResult::Time)           = time;
    results.back()->value(time, SimulationResult::SpinX)          = current.total.spin(0);
    results.back()->value(time, SimulationResult::SpinY)          = current.total.spin(1);
    results.back()->value(time, SimulationResult::SpinZ)          = current.total.spin(2);
    results.back()->value(time, SimulationResult::Contrast)       = current.total.spin.head<3>().norm();
    results.back()->value(time, SimulationResult::AtomLossFromF1) = current.total.atomLosses(1);
    results.back()->value(time, SimulationResult::AtomLossFromF2) = current.total.atomLosses(2);
    results.back()->value(time, SimulationResult::Phase)          = current.total.phase;
}

void Solution::clear()
{
    this->timeStep = 0;

    current.clear();
    previous.clear();
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
    phase = 0;
}
