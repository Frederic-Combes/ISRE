#ifndef SCRIPT_H
#define SCRIPT_H

#include "tools.h"

#include <QString>
class QScriptEngine;

#include <queue>
#include <vector>

typedef std::queue<DataPoint> DataQueue;

bool processFile(QString filename, SimulationTools & tools);

std::vector<double> makeList(QScriptEngine * engine, QString pName);

#endif // SCRIPT_H
