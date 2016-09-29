#include "result.h"
using namespace std;

#include "strings.h"

#include <QFile>
#include <QTextStream>

#include <QSqlQuery>
#include <QSqlError>
#include <QDateTime>

#include <QScriptValue>
#include <QScriptContext>
#include <QScriptEngine>

#include <QDebug>

/*************************************************************************************************************
 * Class Result
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Public Members
 ************************************************************************************************************/

/*QScriptValue Result::scriptFunctionConstructor(QScriptContext * context, QScriptEngine * engine)
{
    QScriptValue thisObject;

    if(context->isCalledAsConstructor()) {
        thisObject = context->thisObject();
    } else {
        thisObject = engine->newObject();
        thisObject.setPrototype(context->callee().prototype());
    }

    thisObject.setProperty(Strings::Script::Result::Matches, engine->newFunction(scriptFunctionMatches));
    thisObject.setProperty(Strings::Script::Result::Export, engine->newFunction(scriptFunctionExport));
    //thisObject.setProperty(Strings::Script::Result::Erase, engine->newFunction(scriptFunctionErase));

    thisObject.setData(engine->newVariant(QVariant::fromValue(new Result())));

    return thisObject;
}





QScriptValue Result::scriptFunctionErase(QScriptContext * context, QScriptEngine * engine)
{
    // Check argument count
    if(context->argumentCount() > 0) {
        qDebug() << "ScriptObject::Result.erase called with an incorrect parameter count:" << context->argumentCount() << ", expected 0 argument.";
        return QScriptValue::UndefinedValue;
    }

    Result * result = context->thisObject().data().toVariant().value<Result*>();

    if(result == nullptr || result->end() == result->begin()) {
        return QScriptValue::UndefinedValue;
    }

    result->erase();
    return QScriptValue::UndefinedValue;
}

*/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

/*Result::Result(QSqlDatabase & db, unsigned int id)
    : p_saved(true), p_id(id)
{
    (void) db;
}

Result::Result()
    : p_saved(false), p_id(0)
{

}*/

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

/*mutex & Result::lock()
{
    return p_mutex;
}

const map<QString, QString> & Result::metadata() const
{
    return p_metadata;
}

map<QString, QString> & Result::metadata()
{
    return p_metadata;
}

double Result::interpolate(double time, Index index) const
{
    map<double,vector<double>>::const_iterator it = p_data.lower_bound(time);

    if(it == p_data.end()) {
        return 0;
    } else if (it->first == time || it == p_data.begin()){
        return it->second[index];
    } else {
        double t2 = it->first, v2 = it->second[index];
        --it;
        double t1 = it->first, v1 = it->second[index];
        double t = time;

        return ((t2-t) * v1 + (t-t1) * v2) / (t2-t1);
    }
}

vector<double> & Result::at(double time)
{
    auto it = p_data.find(time);

    if(it == p_data.end()) {
        it = p_data.insert({time, vector<double>(size)}).first;
        it->second.at(Time) = time;
    }

    return it->second;
}

vector<double> & Result::operator[](double time)
{
    return at(time);
}

unsigned int Result::write(QSqlDatabase & db, const QString & name)
{
    QSqlQuery query(db);

    if(p_saved) {
        query.prepare("DELETE FROM result WHERE id = :id");
        query.bindValue(":id", p_id);

        if(!query.exec()) {
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec:" << query.lastError().text();
        }

        query.prepare("DELETE from description WHERE id = :id");
        query.bindValue(":id", p_id);

        if(!query.exec()) {
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec:" << query.lastError().text();
        }

    } else {
        query.prepare("INSERT INTO simulation (name, timestamp) VALUES (:n, :t);");
        query.bindValue(":n", name);
        query.bindValue(":t", QDateTime::currentDateTime());

        if(!query.exec()) {
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec:" << query.lastError().text();
        }

        p_saved = true;
        p_id = query.lastInsertId().toUInt();
    }

    query.prepare("INSERT INTO description (id, name, value) VALUES (:id, :n, :v);");

    for(map<QString,QString>::value_type value : p_metadata) {
        query.bindValue(":id", p_id);
        query.bindValue(":n", value.first);
        query.bindValue(":v", value.second);

        if(!query.exec()) {
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec:" << query.lastError().text();
        }
    }

    query.prepare("INSERT INTO result (id, time, s_x, s_y, s_z, contrast, lost_1, lost_2, phase) VALUES (:id, :t, :x, :y, :z, :c, :l1, :l2, :p);");

    for(map<double, std::vector<double>>::value_type value: p_data) {
        query.bindValue(":id", p_id);
        query.bindValue(":t", value.second.at(Result::Time));
        query.bindValue(":x", value.second.at(Result::SpinX));
        query.bindValue(":y", value.second.at(Result::SpinY));
        query.bindValue(":z", value.second.at(Result::SpinZ));
        query.bindValue(":c", value.second.at(Result::Contrast));
        query.bindValue(":l1", value.second.at(Result::AtomLossFromF1));
        query.bindValue(":l2", value.second.at(Result::AtomLossFromF2));
        query.bindValue(":p", value.second.at(Result::Phase));

        if(!query.exec()) {
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec.";
            qDebug() << "SimulationDescription::write:" << query.lastError().text();
        }
    }

    return p_id;
}*/

/*************************************************************************************************************
 * Class SimulationResult
 ************************************************************************************************************/

/*************************************************************************************************************
 * Static Public Members
 ************************************************************************************************************/

bool SimulationResult::convertToColumnIndex(QString columnName, SimulationResult::Index & columnIndex)
{
    if(QString::compare(columnName, Strings::Database::Result::Time, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::Time;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::SX, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::SpinX;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::SY, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::SpinY;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::SZ, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::SpinZ;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::Constrast, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::Contrast;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::AtomLossFromF1, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::AtomLossFromF1;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::AtomLossFromF2, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::AtomLossFromF2;
        return true;
    }

    if(QString::compare(columnName, Strings::Database::Result::Phase, Qt::CaseInsensitive) == 0) {
        columnIndex = SimulationResult::Phase;
        return true;
    }

    qDebug() << "ScriptObject::Result: unrecognized column name:" << columnName;

    return false;
}

/*************************************************************************************************************
 * Static Private Members
 ************************************************************************************************************/

SimulationResult::size_type SimulationResult::toSizeType(Index index)
{
    switch (index) {
    case Time: return 0;
    case SpinX: return 1;
    case SpinY: return 2;
    case SpinZ: return 3;
    case Contrast: return 4;
    case AtomLossFromF1: return 5;
    case AtomLossFromF2: return 6;
    case Phase: return 7;
    }

    // Warning supression
    return 0;
}

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

SimulationResult::SimulationResult(unsigned int numberOfSavedPoints)
{
    p_data.reserve(numberOfSavedPoints * IndexCount);
}

/*************************************************************************************************************
 * Public members
 ************************************************************************************************************/

SimulationResult::mutex_type & SimulationResult::lock()
{
    return p_mutex;
}

double & SimulationResult::value(double time, Index index)
{
    map<double, size_type>::const_iterator it = p_mapTimeToData.find(time);

    if(it == p_mapTimeToData.end()) {
        it = p_mapTimeToData.insert({time, p_data.size()}).first;
        p_data.insert(p_data.end(), IndexCount, 0);
    }

    return p_data[it->second + toSizeType(index)];
}

double SimulationResult::interpolate(double time, Index index) const
{
    if(p_mapTimeToData.empty()) {
        return 0;
    }

    map<double,size_type>::const_iterator it = p_mapTimeToData.lower_bound(time);

    if(it == p_mapTimeToData.end()) {
        it = p_mapTimeToData.end()--;
        return p_data.at(it->second + toSizeType(index));
    } else if (it->first == time || it == p_mapTimeToData.begin()){
        return p_data.at(it->second + toSizeType(index));
    } else {
        double t2 = it->first;
        double v2 = p_data.at(it->second + toSizeType(index));

        (--it);

        double t1 = it->first;
        double v1 = p_data.at(it->second + toSizeType(index));

        double t = time;

        return ((t2-t) * v1 + (t-t1) * v2) / (t2-t1);
    }
}

SimulationResult::time_iterator SimulationResult::begin() const
{
    return time_iterator(p_mapTimeToData.cbegin());
}

SimulationResult::time_iterator SimulationResult::end() const
{
    return time_iterator(p_mapTimeToData.cend());
}

void SimulationResult::write(QSqlDatabase & db, unsigned int databaseID) const
{
    QSqlQuery query(db);

    query.prepare("INSERT INTO result (id, time, s_x, s_y, s_z, contrast, lost_1, lost_2, phase) VALUES (:id, :t, :x, :y, :z, :c, :l1, :l2, :p);");

    for(map<double, size_type>::value_type value: p_mapTimeToData) {
        query.bindValue(":id", databaseID);
        query.bindValue(":t",  p_data.at(value.second + toSizeType(Time)));
        query.bindValue(":x",  p_data.at(value.second + toSizeType(SpinX)));
        query.bindValue(":y",  p_data.at(value.second + toSizeType(SpinY)));
        query.bindValue(":z",  p_data.at(value.second + toSizeType(SpinZ)));
        query.bindValue(":c",  p_data.at(value.second + toSizeType(Contrast)));
        query.bindValue(":l1", p_data.at(value.second + toSizeType(AtomLossFromF1)));
        query.bindValue(":l2", p_data.at(value.second + toSizeType(AtomLossFromF2)));
        query.bindValue(":p",  p_data.at(value.second + toSizeType(Phase)));

        if(!query.exec()) {
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec.";
            qDebug() << "SimulationDescription::write:" << query.lastError().text();
        }
    }
}


/*************************************************************************************************************
 * Class SimulationResult::time_iterator
 ************************************************************************************************************/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

SimulationResult::time_iterator::time_iterator(const_iterator it)
    : p_it(it)
{}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

SimulationResult::time_iterator::value_type SimulationResult::time_iterator::operator*()
{
    return p_it->first;
}

/*SimulationResult::time_iterator::reference SimulationResult::time_iterator::operator[](Index index)
{
    return *(p_it->second + toSizeType(index));
}*/

SimulationResult::time_iterator & SimulationResult::time_iterator::operator++()
{
    ++p_it;

    return *this;
}

SimulationResult::time_iterator SimulationResult::time_iterator::operator++(int)
{
    time_iterator it(p_it);

    ++p_it;

    return it;
}

SimulationResult::time_iterator & SimulationResult::time_iterator::operator--()
{
    --p_it;

    return *this;
}

SimulationResult::time_iterator SimulationResult::time_iterator::operator--(int)
{
    time_iterator it(p_it);

    --p_it;

    return it;
}

bool SimulationResult::time_iterator::operator==(const time_iterator & other) const
{
    return p_it == other.p_it;
}

bool SimulationResult::time_iterator::operator!=(const time_iterator & other) const
{
    return p_it != other.p_it;
}

/*************************************************************************************************************
 * Class SimulationDescription
 ************************************************************************************************************/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

SimulationDescription::SimulationDescription()
{}

/*************************************************************************************************************
 * Public members
 ************************************************************************************************************/

void SimulationDescription::setValue(QString name, QString value)
{
    auto it = p_description.find(name);

    if(it == p_description.end()) {
        p_description.insert({name, value});
    } else {
        it->second = value;
    }
}

pair<bool, QString> SimulationDescription::value(QString name) const
{
    auto it = p_description.find(name);

    if(it == p_description.end()) {
        return {false, QString()};
    } else {
        return {true, it->second};
    }
}

void SimulationDescription::write(QSqlDatabase & db, unsigned int databaseID) const
{
    QSqlQuery query(db);

    query.prepare("INSERT INTO description (id, name, value) VALUES (:i, :n, :v)");

    for(const pair<QString, QString> & pairs : p_description) {
        query.bindValue(":i", databaseID);
        query.bindValue(":n", pairs.first);
        query.bindValue(":v", pairs.second);

        if(!query.exec()) {
            qDebug() << "SimulationDescription::write: query" << query.lastQuery() << "failed to exec.";
            qDebug() << "SimulationDescription::write:" << query.lastError().text();
        }
    }
}

