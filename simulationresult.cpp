#include "simulationresult.h"

#include "strings.h"

// STD helpers
#include <functional>

// Basic types
#include <QString>
#include <QRegularExpression>
#include <QDate>

// SQL support
#include <QSqlQuery>
#include <QSqlRecord>
#include <QSqlError>

// File support
#include <QFile>
#include <QTextStream>

// Script Support
#include <QScriptEngine>
#include <QScriptContext>
#include "script.h"

// Debug
#include <QDebug>

using namespace std;

/*************************************************************************************************************
 * SimulationResult
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Public Members
 ************************************************************************************************************/

/*QScriptValue SimulationResult::scriptFunctionTableConstructor(QScriptContext * context, QScriptEngine * engine)
{
    // C++ object that is wrapped by the script
    SimulationResultPtr result = SimulationResultPtr::create(vector<int>());

    // Script object
    QScriptValue thisObject;

    if(context->isCalledAsConstructor()) {
        thisObject = context->thisObject();
    } else {
        thisObject = engine->newObject();
        thisObject.setPrototype(context->callee().property("prototype"));
    }

    // Adds functions to the script object
    thisObject.setProperty("select",    engine->newFunction(scriptFunctionTableSelect));
    thisObject.setProperty("mergeWith", engine->newFunction(scriptFunctionTableMergeWith));
    thisObject.setProperty("joinWith",  engine->newFunction(scriptFunctionTableJoinWith));
    thisObject.setProperty("order",     engine->newFunction(scriptFunctionTableOrder));
    thisObject.setProperty("list",      engine->newFunction(scriptFunctionTableList));
    thisObject.setProperty("write",     engine->newFunction(scriptFunctionTableWrite));
    thisObject.setProperty("erase",     engine->newFunction(scriptFunctionTableErase));

    // Store the C++ object into the script object data
    thisObject.setData(engine->newVariant(QVariant::fromValue(result)));

    // Perfroms the SQL table insertion
    QSqlQuery query(result->p_database);

    query.prepare("INSERT INTO simulation (name, timestamp) VALUES (:n,:t);");

    if(context->argumentCount() > 0 && context->argument(0).isString()) {
        query.bindValue(":n", context->argument(0).toString());
    } else {
         query.bindValue(":n", "no-name");
    }

    query.bindValue(":t", QDate::currentDate().toString("yyyy-MM-dd"));

    if(!query.exec()) {
        qDebug() << "SimulationResult::scriptFunction:TableConstructor: failed to execute the query:" << query.lastQuery().toUtf8().data();
        qDebug() << "SimulationResult::scriptFunction:TableConstructor: error:" << query.lastError().text().toUtf8().data();
    } else {
        qDebug() << "SimulationResult::scriptFunction:TableConstructor: executed the query:" << query.executedQuery().toUtf8().data();
    }

    result->p_ids = {query.lastInsertId().toInt()};

    query.prepare("INSERT INTO description (id, namen value) VALUES (:id, :name, :value)");
    query.bindValue(":id",      query.lastInsertId().toInt());
    query.bindValue(":name",    Strings::Database::Description::FromScript);
    query.bindValue(":value",   true);

    if(!query.exec()) {
        qDebug() << "SimulationResult::scriptFunction:TableConstructor: failed to execute the query:" << query.lastQuery().toUtf8().data();
        qDebug() << "SimulationResult::scriptFunction:TableConstructor: error:" << query.lastError().text().toUtf8().data();
    } else {
        qDebug() << "SimulationResult::scriptFunction:TableConstructor: executed the query:" << query.executedQuery().toUtf8().data();
    }

    if(context->argumentCount() > 1 && context->argument(1).isString()) {
        query.bindValue(":name",    "Comment");
        query.bindValue(":value",   context->argument(1).toString());

        if(!query.exec()) {
            qDebug() << "SimulationResult::scriptFunction:TableConstructor: failed to execute the query:" << query.lastQuery().toUtf8().data();
            qDebug() << "SimulationResult::scriptFunction:TableConstructor: error:" << query.lastError().text().toUtf8().data();
        } else {
            qDebug() << "SimulationResult::scriptFunction:TableConstructor: executed the query:" << query.executedQuery().toUtf8().data();
        }
    }

    return thisObject;
}

QScriptValue SimulationResult::scriptFunctionTableSelect(QScriptContext * context, QScriptEngine * engine)
{
    using namespace Strings::Database;

    // Argument checking
    if(context->argumentCount() < 2) {
        return context->thisObject();
    }

    // Copies the C++ object that is wrapped by the script
    SimulationResultPtr result = SimulationResultPtr::create(*context->thisObject().data().toVariant().value<SimulationResultPtr>());

    if(result.isNull()) {
        return QScriptValue::UndefinedValue;
    }

    // Create a new script object (a table)
    QScriptValue thisObject = engine->newObject();
    thisObject.setPrototype(context->thisObject().prototype());

    // Adds functions to the script object
    // TODO: aren't these functions added with setPrototype?
    thisObject.setProperty("select",    engine->newFunction(scriptFunctionTableSelect));
    thisObject.setProperty("mergeWith", engine->newFunction(scriptFunctionTableMergeWith));
    thisObject.setProperty("joinWith",  engine->newFunction(scriptFunctionTableJoinWith));
    thisObject.setProperty("order",     engine->newFunction(scriptFunctionTableOrder));
    thisObject.setProperty("list",      engine->newFunction(scriptFunctionTableList));
    thisObject.setProperty("write",     engine->newFunction(scriptFunctionTableWrite));
    thisObject.setProperty("erase",     engine->newFunction(scriptFunctionTableErase));

    // Store the C++ object into the script object's data
    thisObject.setData(engine->newVariant(QVariant::fromValue(result)));

    QString name;
    Table table = Description;

    if(context->argument(0).isObject() && context->argument(0).property("name").isString()) {
        // Select by parameter value
        name = context->argument(0).property("name").toString();
    } else if(context->argument(0).isString()) {
        // Select by name,value pair
        name = context->argument(0).toString();

        for(QString columnName : {Result::Time, Result::SX, Result::SY, Result::SZ, Result::Constrast, Result::Phase, Result::AtomLossFromF1, Result::AtomLossFromF2}) {
            if(QString::compare(name, columnName) == 0) {
                // name is the name of a column in the sql table result
                table = Result;
            }
        }
    }

    if(context->argument(1).isFunction()) {
        // Select with a function
        result->select(name, context->argument(1), table);
    } else if(context->argument(1).isArray()) {
        // Select a list of values
        vector<double> values;
        for(unsigned int i = 0; i < context->argument(1).property("length").toUInt32(); ++i) {
            values.push_back(context->argument(1).property(i).toNumber());
        }
        result->select(name, values, table);
    } else if(context->argument(1).isNumber()){
        // Select a single numeric value
        result->select(name, context->argument(1).toNumber(), table);
    } else if(context->argument(1).isString()) {
        // Select a single string value
        result->select(name, context->argument(1).toString(), table);
    }

    return thisObject;
}

QScriptValue SimulationResult::scriptFunctionTableMergeWith(QScriptContext * context, QScriptEngine *)
{
    // Result that is wrapped by the script "this"
    SimulationResultPtr result = context->thisObject().data().toVariant().value<SimulationResultPtr>();

    if(result.isNull()) {
        return QScriptValue::UndefinedValue;
    }

    for(int i = 0; i < context->argumentCount(); ++i) {
        SimulationResultPtr mergeWith = context->argument(i).data().toVariant().value<SimulationResultPtr>();

        if(mergeWith.isNull()) {
            continue;
        }

        result->merge(*mergeWith);
    }

    return QScriptValue::UndefinedValue;
}

QScriptValue SimulationResult::scriptFunctionTableJoinWith(QScriptContext * context, QScriptEngine *)
{
    // C++ object wrapped by the script object
    SimulationResultPtr result = context->thisObject().data().toVariant().value<SimulationResultPtr>();

    if(result.isNull()) {
        return QScriptValue::UndefinedValue;
    }

    for(int i = 0; i < context->argumentCount(); ++i) {
        SimulationResultPtr joinWith = context->argument(i).data().toVariant().value<SimulationResultPtr>();

        if(joinWith.isNull()) {
            continue;
        }

        result->join(*joinWith);
    }

    return QScriptValue::UndefinedValue;
}

QScriptValue SimulationResult::scriptFunctionTableOrder(QScriptContext * context, QScriptEngine *)
{
    // C++ object wrapped by the script object
    SimulationResultPtr result = context->thisObject().data().toVariant().value<SimulationResultPtr>();

    if(result.isNull()) {
        return QScriptValue::UndefinedValue;
    }

    if(context->argumentCount() == 0) {
        // No argument given: all ordering are to be removed
        result->clearOrder();
    } else if(context->argumentCount() == 1) {
        // Only one argument given, default ordering is "ASC"
        result->order(context->argument(0).toString());
    } else {
        result->order(context->argument(0).toString(), QString::compare(context->argument(1).toString(), "ASC") == 0);
    }

    return QScriptValue::UndefinedValue;
}

QScriptValue SimulationResult::scriptFunctionTableList(QScriptContext * context, QScriptEngine * engine)
{
    // C++ object wrapped by the script object
    SimulationResultPtr result = context->thisObject().data().toVariant().value<SimulationResultPtr>();

    if(result.isNull()) {
        return engine->newArray();
    }

    vector<SimulationResult> results = result->list();
    QScriptValue array = engine->newArray(results.size());

    for(unsigned int i = 0; i < results.size(); ++i) {
        // Allocate a new SimulationResult
        SimulationResultPtr singleResult = SimulationResultPtr::create(cref(results.at(i)));

        // Create a new script object
        QScriptValue thisObject = engine->newObject();
        thisObject.setPrototype(context->thisObject().prototype());

        // Adds the function to the "this" object
        thisObject.setProperty("select",    engine->newFunction(scriptFunctionTableSelect));
        thisObject.setProperty("mergeWith", engine->newFunction(scriptFunctionTableMergeWith));
        thisObject.setProperty("joinWith",  engine->newFunction(scriptFunctionTableJoinWith));
        thisObject.setProperty("order",     engine->newFunction(scriptFunctionTableOrder));
        thisObject.setProperty("list",      engine->newFunction(scriptFunctionTableList));
        thisObject.setProperty("write",     engine->newFunction(scriptFunctionTableWrite));
        thisObject.setProperty("erase",     engine->newFunction(scriptFunctionTableErase));

        // Store the result into this.internal_pointer
        thisObject.setData(engine->newVariant(QVariant::fromValue(singleResult)));

        array.setProperty(i, thisObject);
    }

    return array;
}

QScriptValue SimulationResult::scriptFunctionTableWrite(QScriptContext * context, QScriptEngine *)
{
    using namespace Strings::Database;

    // C++ object wrapped by the script object
    SimulationResultPtr result = context->thisObject().data().toVariant().value<SimulationResultPtr>();

    if(result.isNull()) {
        qDebug() << "SimulationResult::scriptFunction:Table.write: unable to retrieve the object data.";
        return QScriptValue::UndefinedValue;
    }

    if(context->argumentCount() < 1 || !context->argument(0).isString()) {
        qDebug() << "SimulationResult::scriptFunction:Table.write: missing first argument (filename)";
        return QScriptValue::UndefinedValue;
    }

    vector<QString> columnNames;

    if(context->argumentCount() < 3) {
        // No column names provided: we save them all
        columnNames = {Result::Time, Result::SX, Result::SY, Result::SZ, Result::Constrast, Result::Phase, Result::AtomLossFromF1, Result::AtomLossFromF2};
    } else if(context->argument(2).isArray()) {
        // Column names provided as an array of string
        for(unsigned int i = 0; i < context->argument(2).property("length").toUInt32(); ++i) {
            if(context->argument(2).property(i).isString()) {
                columnNames.push_back(context->argument(2).property(i).toString());
            }
        }
    } else {
        for(int i = 2; i < context->argumentCount(); ++i) {
            if(context->argument(i).isString()) {
                columnNames.push_back(context->argument(i).toString());
            }
        }
    }

    QString options  = context->argument(1).toString();

    WriteFlags flags = 0;
    if(options.contains(QRegularExpression("[sS]"))) {
        flags |= SeparateIds;
    }
    if(options.contains(QRegularExpression("[nN]"))) {
        flags |= NewLine;
    }
    if(options.contains(QRegularExpression("[oO]"))) {
        flags |= Overwrite;
    }
    if(options.contains(QRegularExpression("[aA]"))) {
        flags &= ~Overwrite;
    }

    result->write(context->argument(0).toString(), flags, columnNames);

    return QScriptValue::UndefinedValue;
}

QScriptValue SimulationResult::scriptFunctionTableErase(QScriptContext * context, QScriptEngine *)
{
    // C++ object wrapped by the script object
    SimulationResultPtr result = context->thisObject().data().toVariant().value<SimulationResultPtr>();

    if(result.isNull()) {
        return QScriptValue::UndefinedValue;
    }

    result->erase();

    return QScriptValue::UndefinedValue;
}*/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

/*SimulationResult::SimulationResult(vector<int> ids)
    : p_ids(ids)
{
    // Obtain & open the database
    p_database = QSqlDatabase::database("isre-result");

    if(!p_database.isValid()) {
        p_database = QSqlDatabase::addDatabase("QSQLITE", "isre-result");
        p_database.setDatabaseName("results.db");

    }

    openDatatabase();
}

SimulationResult::SimulationResult(const SimulationResult & other)
    : p_ids(other.p_ids)
{
    // Obtain & open the database
    p_database = QSqlDatabase::database("isre-result");

    if(!p_database.isValid()) {
        p_database = QSqlDatabase::addDatabase("QSQLITE", "isre-result");
        p_database.setDatabaseName("results.db");
    }

    openDatatabase();
}*/

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

/*void SimulationResult::select(QString name, QScriptValue function, Table table)
{
    switch(table) {
    case Description:
    {
        vector<int> newIds;

        // Ensures the database is open
        openDatatabase();

        // Selects all values associated to name
        QSqlQuery query(p_database);

        query.prepare("SELECT id, value FROM description WHERE name = :n;");
        query.bindValue(":n", name);

        if(!query.exec()) {
            qDebug() << "SimulationResult::select: failed to execute query:" << query.lastQuery().toUtf8().data();
            qDebug() << "SimulationResult::select: error:" << query.lastError().text().toUtf8().data();
            break;
        } else {
            qDebug() << "SimulationResult::select: executed the query:" << query.executedQuery().toUtf8().data();
        }

        // Iterates over the values
        while(query.next()) {
            int id          = query.record().value("id").toInt();
            QString value   = query.record().value("value").toString();

            // If the id is not in the list of already selected ids, we don't select it
            if(find(p_ids.begin(), p_ids.end(), id) == p_ids.end()) {
                continue;
            }

            // Try to cast the value to a number
            bool isNumber = false;
            double number = value.toDouble(&isNumber);

            QScriptValue isValid;

            if(isNumber) {
                isValid = function.call(QScriptValue(), QScriptValueList() << number);
            } else {
                isValid = function.call(QScriptValue(), QScriptValueList() << value);
            }

            if(isValid.isBool() && isValid.toBool()) {
                newIds.push_back(id);
            }
        }

        // Update the selected ids
        p_ids = newIds;

        break;
    }
    case Result:
        qDebug() << "SimulationResult::select: called with table = Result : not implemented";
        break;
    }
}

void SimulationResult::select(QString name, QString value, Table table)
{
    if(table == Description) {
        selectFromDescription("name = \"" + name + "\" AND value = \"" + value + "\";");
    } else if(table == Result) {

    }
}

void SimulationResult::select(QString name, int value, Table table)
{
    if(table == Description) {
        selectFromDescription("name = \"" + name + "\" AND value = \"" + QString::number(value) + "\";");
    } else if(table == Result) {
        p_SQLWhereClauses.push_back(name + " = " + QString::number(value));
        qDebug() << "SimulationResult::select: added where clause" << QString(name + " = " + QString::number(value)).toUtf8().data();
    }
}

void SimulationResult::select(QString name, vector<int> values, Table table)
{
    QString valueList;

    for(int value : values) {
        valueList += QString::number(value) + ", ";
    }

    valueList.chop(2);

    if(table == Description) {
        selectFromDescription("name = \"" + name + "\" AND value IN (" + valueList + ");");
    } else if(table == Result) {
        p_SQLWhereClauses.push_back(name + " IN (" + valueList + ")");
        qDebug() << "SimulationResult::select: added where clause" << QString(name + " IN (" + valueList + ")").toUtf8().data();
    }
}

void SimulationResult::select(QString name, double value, Table table)
{
    if(table == Description) {
        selectFromDescription("name = \"" + name + "\" AND value = \"" + QString::number(value) + "\";");
    } else if(table == Result) {
        p_SQLWhereClauses.push_back(name + " = " + QString::number(value));
        qDebug() << "SimulationResult::select: added where clause" << QString(name + " = " + QString::number(value)).toUtf8().data();
    }
}

void SimulationResult::select(QString name, vector<double> values, Table table)
{
    QString valueList;

    for(double value : values) {
        valueList += QString::number(value) + ", ";
    }

    valueList.chop(2);

    if(table == Description) {
        selectFromDescription("name = \"" + name + "\" AND value IN (" + valueList + ");");
    } else if(table == Result) {
        p_SQLWhereClauses.push_back(name + " IN (" + valueList + ")");
        qDebug() << "SimulationResult::select: added where clause" << QString(name + " IN (" + valueList + ")").toUtf8().data();
    }
}

void SimulationResult::merge(const SimulationResult & other)
{
    for(int id : other.p_ids) {
        auto it = find(p_ids.begin(), p_ids.end(), id);

        if(it == p_ids.end()) {
            p_ids.push_back(id);
        }
    }

    for(QString whereClause : other.p_SQLWhereClauses) {
        p_SQLWhereClauses.push_back(whereClause);
    }
}

void SimulationResult::join(const SimulationResult & other)
{
    for(int id : other.p_ids) {
        p_ids.push_back(id);
    }

    for(QString whereClause : other.p_SQLWhereClauses) {
        p_SQLWhereClauses.push_back(whereClause);
    }
}

void SimulationResult::order(QString name, bool asc)
{
    p_SQLOrderClauses.push_back(name + (asc ? " ASC" : " DESC"));
}

void SimulationResult::clearOrder()
{
    p_SQLOrderClauses.clear();
}

vector<SimulationResult> SimulationResult::list() const {
    vector<SimulationResult> ret;

    for(int id : p_ids) {
        SimulationResult result({id});
        result.p_SQLWhereClauses = p_SQLWhereClauses;

        ret.push_back(result);
    }

    return ret;
}

void SimulationResult::write(QString filename, WriteFlags flags, vector<QString> columns) const
{
    if(p_ids.empty()) {
        return;
    }

    // Opens the file
    QFile file(filename);

    QIODevice::OpenMode openMode = QIODevice::ReadWrite;
    if(flags & Overwrite) {
        openMode |= QIODevice::Truncate;
    } else {
        openMode |= QIODevice::Append;
    }

    if(!file.open(openMode)) {
        qDebug() << "SimulationResult::write: failed to open the file" << filename.toUtf8().data();
        return;
    }

    // Prepares the stream
    QTextStream stream(&file);

    stream.setFieldWidth(20);
    stream.setRealNumberPrecision(12);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    // Appends a new line at the end of the file if the file is not empty
    if(file.size() != 0 && flags & NewLine) {
        stream << endl;
    }

    if(flags & SeparateIds) {
        for(int id : p_ids) {
            write({id}, ref(stream), columns);

            if(flags & NewLine) {
                stream << endl;
            }
        }
    } else {
        write(p_ids, ref(stream), columns);

        if(flags & NewLine) {
            stream << endl;
        }
    }

    file.close();
}

void SimulationResult::erase() const
{
    if(p_ids.empty()) {
        return;
    }

    // Comma-separated list of ids
    QString idList;

    for(int id : p_ids) {
        idList += QString::number(id) + ", ";
    }
    idList.chop(2);

    // WHERE statement
    QString whereClauses;
    whereClauses = "WHERE id IN (" + idList + ")";

    for(const QString & whereClause : p_SQLWhereClauses) {
        whereClauses += " AND " + whereClause;
    }

    // Ensure the database is open
    openDatatabase();

    QSqlQuery query(p_database);
    query.prepare("DELETE FROM result " + whereClauses + ";");

    if(!query.exec()) {
        qDebug() << "SimulationResult::erase: error in the query:" << query.lastQuery().toUtf8().data();
        qDebug() << "SimulationResult::erase:" << query.lastError().text().toUtf8().data();
    } else {
        qDebug() << "SimulationResult::erase: executed the query:" << query.executedQuery().toUtf8().data();
    }

}*/

/*************************************************************************************************************
 * Private Members
 ************************************************************************************************************/

/*void SimulationResult::openDatatabase() const
{
    if(p_database.isOpen()) {
        return;
    }

    if(!const_cast<SimulationResult*>(this)->p_database.open()) {
        qDebug() << "SimulationResult::openDatabase: Failed to open the database. Queries will fail";
    }
}

void SimulationResult::selectFromDescription(QString whereClause)
{
    // Will contains the id that are in p_ids and in the list returned by the query
    vector<int> newIds;

    // Ensure the database is open
    openDatatabase();

    QSqlQuery query(p_database);
    query.prepare("SELECT DISTINCT id FROM description WHERE " + whereClause);


    if(!query.exec()) {
        qDebug() << "SimulationResult::select: error in the query:" << query.lastQuery().toUtf8().data();
        qDebug() << "SimulationResult::select:" << query.lastError().text().toUtf8().data();
        return;
    } else {
        qDebug() << "SimulationResult::select: executed the query:" << query.executedQuery().toUtf8().data();
    }

    while(query.next()) {

        auto it = find(p_ids.begin(), p_ids.end(), query.record().value("id").toInt());

        if(it != p_ids.end()) {
            newIds.push_back(query.record().value("id").toInt());
        }
    }

    p_ids = newIds;
}

void SimulationResult::write(std::vector<int> ids, QTextStream & stream, const vector<QString> & columns) const
{
    // WHERE id IN (...) clause
    QString idList;

    for(int id : ids) {
        idList += QString::number(id) + ", ";
    }
    idList.chop(2);

    // WHERE clauses
    QString whereClauses;
    for(const QString & whereClause : p_SQLWhereClauses) {
        whereClauses += " AND " + whereClause;
    }

    // ORDER BY clause
    QString OrderByClauses;
    if(!p_SQLOrderClauses.empty()) {
        OrderByClauses += " ORDER BY ";
        for(const QString & orderClause : p_SQLOrderClauses) {
            OrderByClauses += orderClause + ", ";
        }
        OrderByClauses.chop(2);
    }

    // Ensures the database is open
    openDatatabase();

    QSqlQuery query(p_database);
    query.prepare("SELECT * FROM result WHERE id IN (" + idList + ")" + whereClauses + OrderByClauses + ";");

    if(!query.exec()) {
        qDebug() << "SimulationResult::write: error in the query:" << query.lastQuery().toUtf8().data();
        qDebug() << "SimulationResult::write:" << query.lastError().text().toUtf8().data();
    } else {
        qDebug() << "SimulationResult::write: executed the query:" << query.executedQuery().toUtf8().data();
    }

    int fieldWidth = stream.fieldWidth();

    while(query.next()) {
        for(const QString & column : columns) {
            stream << query.record().value(column).toDouble();
        }
        stream << qSetFieldWidth(0) << endl << qSetFieldWidth(fieldWidth);
    }
}*/
