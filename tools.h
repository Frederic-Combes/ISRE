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

// Linear algebra
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#pragma GCC diagnostic pop

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

    /* Atom loss to spectator states */
    double  lossFromF1;
    double  lossFromF2;

    int                 simulationId;
    std::vector<int>    energyBinId;
};

struct Operation
{
    bool            isFixedTime;
    double          time;

    Eigen::Matrix3d operationMatrix;
};

class SimulationTools
{
public:
    typedef std::pair<double,double> EnergyBin;

    SimulationTools();
    SimulationTools(const SimulationTools &) = delete;
    ~SimulationTools();

    /* Initialization
     * calling init() computes the time steps at which a spin echo must be made
     */
    void init();

    /* Data queue
     * operations on queue() must be made inside a queueMutex()->lock() ... queueMutex()->unlock() block
     */
    inline std::mutex * queueMutex() {return p_queueMutex;}
    inline std::queue<DataPoint> & queue() {return p_dataQueue;}

    /* Space dimension */
    void setSpaceDimension(int dimension); // Also updates the energy sampling

    /* Time settings
     * total time is the lenght of the simulation in second, time steps count is the number of time steps
     * used to computed the numerical solution, while time points count is the number of points that
     * will be saved.
     * When spin echo is used, when we compute the result at t = a*time, then the number of time steps
     * used for the numerical solution is a*timeStepCount(), i.e. the lenght (in sec) of the time step
     * is the same, but not the number of time steps
     */
    inline void setTotalTime(double time) {p_totalTime = time;}
    inline double totalTime() const {return p_totalTime;}
    inline void setTimeStepsCount(int count) {p_timeStepsCount = count;}
    inline int timeStepsCount() const{return p_timeStepsCount;}
    inline void setTimePointsCount(int count) {p_timePointsCount = count;}
    inline int timePointsCount() const {return p_timePointsCount;}

    /* Energy settings */
    inline int energyCount() const {return p_energyCount;}
    inline double energy(int i) const{ return p_energy[i];}
    inline void addEnergyBin(double e1, double e2) {p_energyBin.push_back(std::make_pair(e1 < e2 ? e1 : e2, e1 < e2 ? e2: e1));}
    inline int energyBinCount() const {return p_energyBin.size();}
    inline EnergyBin energyBin(int i) const {return p_energyBin.at(i);}
    inline int energyBinWeight(int i) const {return p_energyBinWeight.at(i);}

    /* Energy sampling related functions */
    void sampleEnergy(int spaceDimension, int sampleCount);
    void computeEnergyBinWeights();

    inline bool hasMovingOperation() {return p_hasMovingOperations;}
    void addOperation(const Operation & operation);
    const Operation * nextOperation(double currentTime, double totalTime) const;

    /* Spin echo settings */
    // TODO: change name to 'moving' spin echo or something
    inline void enableSpinEcho(bool enable) {p_spinEcho = enable;}                          // TODO: when enabled is true, spin echoes are moving
    inline void setSpinEchoTime(std::vector<double> time) {p_tmp_spinEchoTime = time;}
    inline bool spinEcho() const {return p_spinEcho;}
    inline int spinEchoCount() const {return p_spinEchoTimeStep.size();}
    inline double spinEchoTimeStep(int i) const {return p_spinEchoTimeStep.at(i);}

    /* Time dependence of the density
     * if the initial density is n_0, then at the time step t, the density is n_0 * densityFactor(t).
     * Therefore all quantities proportional to or depending on the density must be multiplied by the
     * correct power of densityFactor(t) = densityAtTimestep(t)
     */
    inline std::vector<double> & densityTimeDependance() {return p_densityTimeDependence;}
    inline double densityAtTimestep(int timestep) const {return p_densityTimeDependence.at(timestep);}

    /* Energy dependence of the density (for the mean-field density inhomogeneity) */
    inline std::vector<double> & densityEnergyDependance() {return p_densityEnergyDependence;}
    inline double densityAtEnergy(int energy) const {return p_densityEnergyDependence.at(energy);}

    /* Database */
    inline std::mutex * databaseMutex(){return p_databaseMutex;}
    inline const QSqlDatabase & database() const {return p_database;}
    inline const QString & simulationName() const {return p_simulationName;}
    inline void setSimulationName(const QString & name) {p_simulationName = name;}

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

    bool                        p_hasMovingOperations;
    std::vector<Operation>      p_operations;

    bool                        p_spinEcho;             // TODO: remove, use new p_hasMovingOperations
    std::vector<int>            p_spinEchoTimeStep;
    std::vector<double>         p_tmp_spinEchoTime;

    std::vector<double>         p_densityTimeDependence;
    std::vector<double>         p_densityEnergyDependence;

    std::mutex                  *p_databaseMutex;
    QSqlDatabase                p_database;
    QString                     p_simulationName;
};

// Solution of the equation, valid up to the time timestep
class Solution
{
public:
    Solution(const SimulationTools & simulation);
    ~Solution();
private:
    Solution(const Solution &) = delete;
    Solution & operator=(const Solution &) = delete;

public:
    int                 timeStep;

    std::vector<double>   spin_x;                        // Spin S_x(E,t=timeStep)
    std::vector<double>   spin_y;                        // Spin S_y(E,t=timeStep)
    std::vector<double>   spin_z;                        // Spin S_z(E,t=timeStep)
    std::vector<double>   spin_n;                        // Population at energy E in states F=(1,2), mf=0 (norm of the spin)

    double                spectatorPopulationF1;          // Population in states F=2, mf=(-1,1)
    double                spectatorPopulationF2;          // Population in states F=2, mf=(-2,-1,1,2)

    double                totalSpin_x, totalSpin_x_old;   // Total spin S_x(t=timeStep), S_x(t=timeStep-1)
    double                totalSpin_y, totalSpin_y_old;   // Total spin S_y(t=timeStep), S_y(t=timeStep-1)
    double                totalSpin_z, totalSpin_z_old;   // Total spin S_z(t=timeStep), S_z(t=timeStep-1)
    double                totalSpin_n, totalSpin_n_old;   // Sum of the norm of the spin (total population in F=(1,2), mf=0)

    /* The values of t for which the total spin is saved are defined by simulation.time.points */
    double                *totalSpin_x_save;              // Values of the total spin S_x(t) that will be saved
    double                *totalSpin_y_save;              // Values of the total spin S_y(t) that will be saved
    double                *totalSpin_z_save;              // Values of the total spin S_z(t) that will be saved
    double                *spectatorPopulationF1_save;    // Population in states F=2, mf=(-1,1)
    double                *spectatorPopulationF2_save;    // Population in states F=2, mf=(-2,-1,1,2)

    /* The values of t for which the partial spins are saved are defined by simulation.time.points */
    double                **partialSpin_x_save;           // Values of the partial spin \sum_{E in [E1;E2]} S_x(E,t) that will be saved
    double                **partialSpin_y_save;           // Values of the partial spin \sum_{E in [E1;E2]} S_y(E,t) that will be saved
    double                **partialSpin_z_save;           // Values of the partial spin \sum_{E in [E1;E2]} S_z(E,t) that will be saved

private:
    unsigned int p_energyBinCount;
};


#endif // THREADTOOLDS_H
