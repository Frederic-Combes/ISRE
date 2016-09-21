#ifndef SCRIPT_H
#define SCRIPT_H

#include <vector>
#include <set>

#include <QString>

// Script support
#include <QObject>
#include <QScriptable>
#include <QScriptEngine>
#include <QScriptValue>

// Mathematical expressions JIT
#include <mathpresso/mathpresso.h>

class SimulationContext;
class CompiledFunction;
class ScriptLoader;

namespace ScriptObject {
class Parameter;
}

struct FunctionCode
{
    FunctionCode(QString n, QString c) : name(n), code(c) {}

    QString name;
    QString code;
};

class CompiledFunction
{
public:
    CompiledFunction(const QString & name, unsigned int index, const QString & code);

    inline const QString & name() const {return p_name;}
    inline const QString & code() const {return p_code;}
    inline unsigned int index() const {return p_index;}

    inline int flag() const {return p_flag;}
    inline int & flag() {return p_flag;}

    bool buildDependency(const std::vector<ScriptObject::Parameter> & parameters, const std::vector<CompiledFunction *> & functions);
    bool compile(const mathpresso::Context & context);

    bool dependsOn(const QString & name) const;
    const std::set<QString> & dependencies() const {return p_dependencies;}

    inline double eval(double * values) const {return p_compiledCode.evaluate(values);}

private:
    QString                         p_name;
    QString                         p_code;

    unsigned int                    p_index;
    mathpresso::Expression          p_compiledCode;

    std::set<QString>               p_dependencies;

    bool                            p_dependencyBuilt;
    bool                            p_valid;

    int                             p_flag;
};

// TODO: use unique ptr
typedef std::vector<CompiledFunction*> CompiledFunctions;

class ScriptLoader
{
public: // Script-C++ bindings
    static QScriptValue scriptFunctionLoadScript(QScriptContext * context, QScriptEngine * engine);

    static QScriptValue scriptFunctionSimulationConstructor(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionSimulationRun(QScriptContext * context, QScriptEngine * engine);

public:
    static void init();
    static bool load(QString script);

private:
    static bool performSyntaxCheck(const QString & script);
    static bool evaluate(const QString & script);

private:
    static QScriptEngine * p_engine;
};

namespace Index
{
    enum {Frequency = 0, Exchange = 1, Damping = 2, AtomLossFromF1 = 3, AtomLossFromF2 = 4,
          Timestep = 5, EnergyIndex = 6, Time = 7, Energy = 8, Duration = 9, TimestepCount = 10,
          InitialSpinX = 11, InitialSpinY = 12, InitialSpinZ = 13, NextSave = 14,
          Start = 15};
}

namespace Compute
{
    enum {Never, AtStart, OnceForEachStep, OnceForEachEnergy, CacheNow, CacheAtStart};
}

#endif // SCRIPT_H
