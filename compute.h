#ifndef COMPUTE_H
#define COMPUTE_H

#include "simulation.h"
#include "solution.h"
class SimulationContext;

namespace ScriptObject {
class Settings;
class Operation;
}

class Result;

// Thread loop for multithreading
void threadLoop(SimulationContext & context);

// Initial conditions to allow step-by-step solving
void initialize(double * data, Solution & solution, SimulationContext & context, const std::vector<std::shared_ptr<SimulationResult>> & results);

// Performs a time step
void step(double * data, Solution & solution, const SimulationContext & context, const std::vector<std::shared_ptr<SimulationResult>> & results);

// Sate the data if necessary
void save(double * data, Solution & solution, const ScriptObject::Settings & settings, const std::vector<std::shared_ptr<SimulationResult>> & results);

// Computes the phase
void phase(Solution & solution, const ScriptObject::Settings & settings);

// Applies an operation (pi-pulse, ...)
void applyOperation(const ScriptObject::Operation & operation, Solution & solution, const ScriptObject::Settings & settings);

#endif // COMPUTE_H
