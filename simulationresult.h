#ifndef SIMULATIONRESULT_H
#define SIMULATIONRESULT_H

// STD helpers
#include <vector>

// Basic types
#include <QString>
#include <QFlags>

// File support
#include <QTextStream>

// SQL support
#include <QSqlDatabase>
#include <QSqlQuery>

// Script support
#include <QScriptValue>
class QScriptEngine;
class QScriptContext;
#include <QSharedPointer>

// TODO: Rename to TableScriptObject (?)
class SimulationResult
{
public:
    enum Table {Result, Description};
    enum WriteFlag {Overwrite = 0x01, NewLine = 0x02, SeparateIds = 0x04};
    Q_DECLARE_FLAGS(WriteFlags, WriteFlag)

public:
    static QScriptValue scriptFunctionTableConstructor(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableSelect(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableMergeWith(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableJoinWith(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableOrder(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableList(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableWrite(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionTableErase(QScriptContext * context, QScriptEngine * engine);

public:
    SimulationResult(std::vector<int> ids);
    SimulationResult(const SimulationResult & other);

    void select(QString name, QScriptValue function,      Table table = Description);
    void select(QString name, QString value,              Table table = Description);
    void select(QString name, int value,                  Table table = Description);
    void select(QString name, std::vector<int> values,    Table table = Description);
    void select(QString name, double value,               Table table = Description);
    void select(QString name, std::vector<double> values, Table table = Description);

    // Merges the id list and removes the duplicates
    void merge(const SimulationResult & other);
    // Joins the id list, keeps the duplicates
    void join(const SimulationResult & other);

    void order(QString name, bool asc = true);
    void clearOrder();

    std::vector<SimulationResult> list() const;

    void write(QString filename, WriteFlags flags, std::vector<QString> columns) const;

    void erase() const;

private:
    void selectFromDescription(QString whereClause);
    void openDatatabase() const;

    void write(std::vector<int> ids, QTextStream & stream, const std::vector<QString> & columns) const;

private:
    QSqlDatabase            p_database;

    std::vector<int>        p_ids;

    std::vector<QString>    p_SQLWhereClauses;
    std::vector<QString>    p_SQLOrderClauses;

    QSharedPointer<QSqlQuery> p_query;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(SimulationResult::WriteFlags)

// TODO: use std::shared_ptr (?) or better, unique_ptr
typedef QSharedPointer<SimulationResult> SimulationResultPtr;
Q_DECLARE_METATYPE(SimulationResultPtr)

// TODO: Rename to LineScriptObject (?) & implement
class SimulationResultLine
{
public:
    static QScriptValue scriptFunctionConstructor(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionValue(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionSet(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionPrevious(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionNext(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionExists(QScriptContext * context, QScriptEngine * engine);
};

// TODO: use std::shared_ptr (?)
typedef QSharedPointer<QSqlQuery> QSQLQueryPtr;
Q_DECLARE_METATYPE(QSQLQueryPtr)

#endif // SIMULATIONRESULT_H
