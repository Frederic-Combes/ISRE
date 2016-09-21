#ifndef RESULT_H
#define RESULT_H

#include <vector>
#include <map>
#include <memory>

#include <mutex>

#include <QSqlDatabase>
#include <QString>

class Result
{    
public:
    enum Index {Time = 0, SpinX = 1, SpinY = 2, SpinZ = 3, Contrast = 4, AtomLossFromF1 = 5, AtomLossFromF2 = 6, Phase = 7};

    static unsigned int toInt(Index index);

private :
    static constexpr unsigned int size = 8;

public:
    Result();
    Result(QSqlDatabase & db, unsigned int id);

    std::mutex & lock();

    const std::map<QString,QString> & metadata() const;
    std::map<QString,QString> & metadata();

    double interpolate(double time, Index index) const;

    std::vector<double> & operator[](double time);
    std::vector<double> & at(double time);

    /**
     * @brief Writes (or updates) the results to the database.
     *
     * @warning The database must be open (and preferably the call must be enclosed in a transaction)
     *
     * @param db database
     * @param name name of the simulation
     */
    unsigned int write(QSqlDatabase & db, const QString & name);

private:
    std::mutex p_mutex;

    bool p_saved;
    unsigned int p_id;
    std::map<QString, QString> p_metadata;

    std::map<double, std::vector<double>> p_data;
};

typedef std::unique_ptr<Result> ResultPtr;

#endif // RESULT_H
