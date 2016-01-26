#ifndef THREADTOOLDS_H
#define THREADTOOLDS_H

#include <QString>

#include <fstream>
#include <queue>
#include <map>
#include <vector>
#include <mutex>
#include <algorithm>

#include <QSqlDatabase>

typedef double real;

// Set of parameters to simulate
struct DataPoint
{
    double  totalTime;
    int     timeStepCount;

    double  larmor;
    double  larmorInhomogeneity;
    double  densityInhomogeneity;
    double  exchange;
    double  damping;

    int                 simulationId;
    std::vector<int>    energyBinId;
};

// To allow the use of std:map<DataPoint, ... >
bool operator<(const DataPoint & a, const DataPoint & b);

// Solution of the equation, valid up to the time timestep
struct Solution
{
    int                 timeStep;

    real                *spin_x;                        // Spin S_x(E,t=timeStep)
    real                *spin_y;                        // Spin S_y(E,t=timeStep)
    real                *spin_z;                        // Spin S_z(E,t=timeStep)

    real                totalSpin_x, totalSpin_x_old;   // Total spin S_x(t=timeStep), S_x(t=timeStep-1)
    real                totalSpin_y, totalSpin_y_old;   // Total spin S_y(t=timeStep), S_y(t=timeStep-1)
    real                totalSpin_z, totalSpin_z_old;   // Total spin S_z(t=timeStep), S_z(t=timeStep-1)

    real                *totalSpin_x_save;              // Values of the total spin S_x(t) that will be saved in the file
    real                *totalSpin_y_save;              // Values of the total spin S_y(t) that will be saved in the file
    real                *totalSpin_z_save;              // Values of the total spin S_z(t) that will be saved in the file
    // The values of t for which the total spin is saved are defined by simulation.time.points

    real                **partialSpin_x_save;           // Values of the partial spin \sum_{E in [E1;E2]} S_x(E,t) that will be saved in the file
    real                **partialSpin_y_save;           // Values of the partial spin \sum_{E in [E1;E2]} S_y(E,t) that will be saved in the file
    real                **partialSpin_z_save;           // Values of the partial spin \sum_{E in [E1;E2]} S_z(E,t) that will be saved in the file
    // The values of t for which the partial spins are saved are defined by simulation.time.points
};

class ThreadTools
{
public:
    struct SpinEchoResult
    {
        real                time;

        real                totalSpin_x;
        real                totalSpin_y;
        real                totalSpin_z;

        std::vector<real>   partialSpin_x;
        std::vector<real>   partialSpin_y;
        std::vector<real>   partialSpin_z;
    };

    typedef std::pair<double,double> EnergyBin;

private:
    // Time ordering of SpinEchoResult (when spin echo is used)
    struct {bool operator()(const SpinEchoResult & a, const SpinEchoResult & b) {return a.time < b.time;}} spinEchoResultTimeOrder;

public:
    ThreadTools();
private:
    ThreadTools(const ThreadTools &);
public:
    ~ThreadTools();

    /* Initialization
     * calling init() computes the energy samples, the weight of the energy bins and the time steps at
     * which a spin echo must be made
     */
    void                    init();

    /* Data queue
     * operations on queue() must be made inside a queueMutex()->lock() ... queueMutex()->unlock() block
     */
    inline std::mutex *             queueMutex() {return p_queueMutex;}
    inline std::queue<DataPoint> &  queue() {return p_dataQueue;}

    /* Space dimension */
    void                    setSpaceDimension(int dimension) {p_dimension = dimension;}

    /* Time settings
     * total time is the lenght of the simulation in second, time steps count is the number of time steps
     * used to computed the numerical solution, while time points count is the number of points that
     * will be saved.
     * When spin echo is used, when we compute the result at t = a*time, then the number of time steps
     * used for the numerical solution is a*timeStepCount(), i.e. the lenght (in sec) of the time step
     * is the same, but not the number of time steps
     */
    void                    setTotalTime(double time) {p_totalTime = time;}
    inline double           totalTime() const {return p_totalTime;}
    void                    setTimeStepsCount(int count) {p_timeStepsCount = count;}
    inline int              timeStepsCount() const{return p_timeStepsCount;}
    void                    setTimePointsCount(int count) {p_timePointsCount = count;}
    inline int              timePointsCount() const {return p_timePointsCount;}

    /* Energy settings */
    void                    setEnergyCount(int count) {p_energyCount = count;}
    inline int              energyCount() const{return p_energyCount;}
    inline double           energy(int i) const{ return p_energy[i];}
    void                    addEnergyBin(double e1, double e2) {p_energyBin.push_back(std::make_pair(e1 < e2 ? e1 : e2, e1 < e2 ? e2: e1));}
    inline int              energyBinCount() const {return p_energyBin.size();}
    inline EnergyBin        energyBin(int i) const {return p_energyBin.at(i);}
    inline int              energyBinWeight(int i) const {return p_energyBinWeight.at(i);}

    /* Spin echo settings */
    void                    enableSpinEcho(bool enable) {p_spinEcho = enable;}
    void                    setSpinEchoTime(std::vector<double> time) {p_tmp_spinEchoTime = time;}
    inline bool             spinEcho() const {return p_spinEcho;}
    inline int              spinEchoCount() const {return p_spinEchoTimeStep.size();}
    inline double           spinEchoTimeStep(int i) const {return p_spinEchoTimeStep.at(i);}

    /* Time dependence of the density
     * if the initial density is n_0, then at the time step t, the density is n_0 * densityFactor(t).
     * Therefore all quantities proportional to or depending on the density must be multiplied by the
     * correct power of densityFactor(t)
     */
    void                    setDensityFactor(double * data) {p_densityFactor = data;}
    inline double           densityFactor(int timestep) {return p_densityFactor[timestep];}

    /* File stream with automatic field width setting
     * file operations must be made within a fileMutex()->lock() ... fileMutex()->unlock() block.
     * file() has automatic width setting, while rawfile() does not
     */
    bool                        setOutput(QString filename);
    inline std::mutex *         fileMutex(){return p_fileMutex;}
    inline std::fstream &       file() {p_file.width(20);return p_file;}
    inline std::fstream &       rawFile() {return p_file;}

    /* Database*/
    inline const QSqlDatabase & database() const {return p_database;}
    inline const QString &      simulationName() const {return p_simulationName;}
    inline void                 setSimulationName(const QString & name) {p_simulationName = name;}

    /* Time ordering of the result when spin echo is used
     * when spin echo is used, each thread computes only one value in the time dependence. When all points
     * have been computed for a given data set (of identical larmor, inhomeneity, exchange and damping)
     * the function allResultPointsComputed returns true and result returns the time ordered result (due
     * to multithreading, values are usually not computed in a time ordered fashion)
     *
     * this 4 methods must be called inside a resultMutex()->lock() ... resultMutex()->unlock() block
     */
    inline std::mutex *         resultMutex() {return p_resultMutex;}
    void                        insertResultPoint(const DataPoint & dataPoint, const Solution & solution);
    bool                        allResultPointsComputed(const DataPoint & dataPoint) const;
    std::vector<SpinEchoResult> result(const DataPoint & dataPoint);
    void                        removeResult(const DataPoint & dataPoint);// Thread safe function

private :
    std::mutex                  *p_queueMutex;
    std::mutex                  *p_fileMutex;
    std::mutex                  *p_resultMutex;

    std::queue<DataPoint>       p_dataQueue;

    int                         p_dimension;

    double                      p_totalTime;
    int                         p_timeStepsCount;
    int                         p_timePointsCount;

    int                         p_energyCount;
    double                      *p_energy;
    std::vector<EnergyBin>      p_energyBin;
    std::vector<int>            p_energyBinWeight;

    bool                        p_spinEcho;
    std::vector<int>            p_spinEchoTimeStep;
    std::vector<double>         p_tmp_spinEchoTime;

    double                      *p_densityFactor;

    std::fstream                p_file;
    QSqlDatabase                p_database;
    QString                     p_simulationName;

    std::map<DataPoint,std::vector<SpinEchoResult>>  p_spinEchoResults;
};


#endif // THREADTOOLDS_H
