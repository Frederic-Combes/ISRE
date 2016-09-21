#include "solution.h"
#include "script.h"
#include "simulation.h"

#include <QCoreApplication>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <QDebug>

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include <queue>
#include <map>
#include <vector>
#include <thread>
#include <mutex>

#include <QtSql>

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#pragma GCC diagnostic pop
using namespace Eigen;

constexpr double PI = M_PI;

using namespace std;

int main(int argc, char * argv[])
{
    // Necessary to use QScriptEngine and QSqlDatabase
    QCoreApplication app(argc, argv); Q_UNUSED(app);

    // Check that a script name was passed to the application
    if(argc < 2) {
        qDebug() << "main: No script provided.";
        return -1;
    }

    QFile scriptFile(argv[1]);

    if(!scriptFile.open(QIODevice::ReadOnly)) {
        qDebug() << "main: failed to open the file" << argv[1];
        return -1;
    }

    // Load the script
    ScriptLoader::init();
    ScriptLoader::load(QTextStream(&scriptFile).readAll());

    return 0;
}
