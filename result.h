#ifndef RESULT_H
#define RESULT_H

#include <vector>
#include <map>
#include <memory>

#include <mutex>

#include <QSqlDatabase>
#include <QString>

#include <QScriptValue>
class QScriptContext;
class QScriptEngine;

/*class Result
{    
public:
    enum Index {Time = 0, SpinX = 1, SpinY = 2, SpinZ = 3, Contrast = 4, AtomLossFromF1 = 5, AtomLossFromF2 = 6, Phase = 7};

    static QScriptValue scriptFunctionConstructor(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionMatches(QScriptContext * context, QScriptEngine * engine);
    static QScriptValue scriptFunctionExport(QScriptContext * context, QScriptEngine * engine);
    //static QScriptValue scriptFunctionErase(QScriptContext * context, QScriptEngine * engine);

private :
    static constexpr unsigned int size = 8;
    static bool convertToColumnIndex(QString columnName, Index & columnIndex);

public:
    Result();
    Result(QSqlDatabase & db, unsigned int id);

    std::mutex & lock();

    const std::map<QString,QString> & metadata() const;
    std::map<QString,QString> & metadata();

    double interpolate(double time, Index index) const;

    std::vector<double> & operator[](double time);
    std::vector<double> & at(double time);

    std::map<double, std::vector<double>>::const_iterator begin() const;
    std::map<double, std::vector<double>>::const_iterator end() const;

    //
     * @brief Writes (or updates) the results to the database.
     *
     * @warning The database must be open (and preferably the call must be enclosed in a transaction)
     *
     * @param db database
     * @param name name of the simulation
    //
    unsigned int write(QSqlDatabase & db, const QString & name);
    //void erase(QSqlDatabase & db) {}

private:
    std::mutex p_mutex;

    bool p_saved;
    unsigned int p_id;
    std::map<QString, QString> p_metadata;

    std::map<double, std::vector<double>> p_data;
};

typedef std::unique_ptr<Result> ResultPtr;*/

/**
 * @brief Holds the result of a simulation (C++ side)
 */
class SimulationResult
{
    static constexpr unsigned int IndexCount = 8;
    typedef std::vector<double>::size_type size_type;
public:
    enum Index {Time, SpinX, SpinY, SpinZ, Contrast, AtomLossFromF1, AtomLossFromF2, Phase};
    static bool convertToColumnIndex(QString columnName, SimulationResult::Index & columnIndex);

    typedef std::mutex mutex_type;
    typedef std::lock_guard<mutex_type> lock_type;

    class time_iterator : public std::iterator<std::bidirectional_iterator_tag, double>
    {
    public:
        typedef std::map<double, size_type>::const_iterator const_iterator;

        explicit time_iterator(const_iterator it);

        value_type operator*();
        //reference operator[](Index index);

        time_iterator & operator++();
        time_iterator operator++(int);

        time_iterator & operator--();
        time_iterator operator--(int);

        bool operator==(const time_iterator & other) const;
        bool operator!=(const time_iterator & other) const;

    private:
        const_iterator p_it;
    };

    /**
     * @brief Allocate storage for numberOfSavedPoints records
     */
    SimulationResult(unsigned int numberOfSavedPoints);

    /**
     * @brief lock for write acces
     */
    mutex_type & lock();

    /**
     * @brief returns a ref to the value at time time and with index index. Performs an insertion if
     * necessary
     *
     * No iterators are invalidated
     */
    double & value(double time, Index index);
    /**
     * @brief returns a ref to the value at time time and with index index, if any. If there is no value for
     * time time, interpolates between the two nearest time points.
     */
    double interpolate(double time, Index index) const;

    /**
     * @brief allows to iterate over time values of the container
     */
    time_iterator begin() const;
    time_iterator end() const;

    /**
     * @brief writes the result to the database db, using the database ID databaseID
     */
    void write(QSqlDatabase & db, unsigned int databaseID) const;

private:
    static size_type toSizeType(Index index);

private:
    mutex_type p_mutex;

    std::map<double, size_type> p_mapTimeToData;
    std::vector<double> p_data;
};

/**
 * @brief Holds the description of a simulation (C++ side)
 */
class SimulationDescription
{
public:
    SimulationDescription();

    void setValue(QString name, QString value);
    std::pair<bool, QString> value(QString name) const;

    void write(QSqlDatabase & db, unsigned int databaseID) const;

private:
    std::map<QString, QString> p_description;
};

#endif // RESULT_H
