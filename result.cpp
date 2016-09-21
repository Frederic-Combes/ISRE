#include "result.h"
using namespace std;

#include <QSqlQuery>
#include <QSqlError>
#include <QDateTime>

#include <QDebug>

/*************************************************************************************************************
 * Static Members
 ************************************************************************************************************/

/*************************************************************************************************************
 * Constructors
 ************************************************************************************************************/

Result::Result(QSqlDatabase & db, unsigned int id)
    : p_saved(true), p_id(id)
{
    (void) db;
}

Result::Result()
    : p_saved(false), p_id(0)
{

}

/*************************************************************************************************************
 * Public Members
 ************************************************************************************************************/

mutex & Result::lock()
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
            qDebug() << "Result::write: query" << query.lastQuery() << "failed to exec:" << query.lastError().text();
        }
    }

    return p_id;
}


