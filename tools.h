#ifndef TOOLS_H
#define TOOLS_H

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

class SimulationTools
{
public:
    typedef std::pair<double,double> EnergyBin;

    SimulationTools();
    SimulationTools(const SimulationTools &) = delete;
    ~SimulationTools();

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
    inline void              setSpaceDimension(int dimension) {p_dimension = dimension;}

    /* Time settings
     * total time is the lenght of the simulation in second, time steps count is the number of time steps
     * used to computed the numerical solution, while time points count is the number of points that
     * will be saved.
     * When spin echo is used, when we compute the result at t = a*time, then the number of time steps
     * used for the numerical solution is a*timeStepCount(), i.e. the lenght (in sec) of the time step
     * is the same, but not the number of time steps
     */
    inline void             setTotalTime(double time) {p_totalTime = time;}
    inline double           totalTime() const {return p_totalTime;}
    inline void             setTimeStepsCount(int count) {p_timeStepsCount = count;}
    inline int              timeStepsCount() const{return p_timeStepsCount;}
    inline void             setTimePointsCount(int count) {p_timePointsCount = count;}
    inline int              timePointsCount() const {return p_timePointsCount;}

    /* Energy settings */
    inline void             setEnergyCount(int count) {p_energyCount = count;}
    inline int              energyCount() const{return p_energyCount;}
    inline double           energy(int i) const{ return p_energy[i];}
    inline void             addEnergyBin(double e1, double e2) {p_energyBin.push_back(std::make_pair(e1 < e2 ? e1 : e2, e1 < e2 ? e2: e1));}
    inline int              energyBinCount() const {return p_energyBin.size();}
    inline EnergyBin        energyBin(int i) const {return p_energyBin.at(i);}
    inline int              energyBinWeight(int i) const {return p_energyBinWeight.at(i);}

    /* Spin echo settings */
    inline void             enableSpinEcho(bool enable) {p_spinEcho = enable;}
    inline void             setSpinEchoTime(std::vector<double> time) {p_tmp_spinEchoTime = time;}
    inline bool             spinEcho() const {return p_spinEcho;}
    inline int              spinEchoCount() const {return p_spinEchoTimeStep.size();}
    inline double           spinEchoTimeStep(int i) const {return p_spinEchoTimeStep.at(i);}

    /* Time dependence of the density
     * if the initial density is n_0, then at the time step t, the density is n_0 * densityFactor(t).
     * Therefore all quantities proportional to or depending on the density must be multiplied by the
     * correct power of densityFactor(t)
     */
    inline void             setDensityFactor(double * data) {p_densityFactor = data;}
    inline double           densityFactor(int timestep) {return p_densityFactor[timestep];}

    /* Database*/
    inline std::mutex *         databaseMutex(){return p_databaseMutex;}
    inline const QSqlDatabase & database() const {return p_database;}
    inline const QString &      simulationName() const {return p_simulationName;}
    inline void                 setSimulationName(const QString & name) {p_simulationName = name;}

private :
    std::mutex                  *p_queueMutex;
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

    //std::fstream                p_file;
    std::mutex                  *p_databaseMutex;
    QSqlDatabase                p_database;
    QString                     p_simulationName;
};


#endif // THREADTOOLDS_H
