#ifndef SOLUTION_H
#define SOLUTION_H

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
#include <Eigen/StdVector>
#pragma GCC diagnostic pop

#include "simulation.h"
#include "scriptobjects.h"
#include "result.h"

// Solution of the equation, valid up to the time timestep
// The class has the following members :

// solution.timestep                    // current timestep

// solution.current.spin(e)             // current spin at energy index e
// solution.current.total.spin          // total current spin
// solution.current.total.phase         // current phase
// solution.current.total.atomLosses(1) // totel current atom losses from the F=1 level
// solution.current.total.atomLosses(2) // total current atom losses from the F=2 level
// solution.current.partial(i).spin     // total current spin in the i-th engery bin
// solution.current.partial(i).phase    // current phase in the i-th energy bin

// solution.previous.total.spin         // total spin at the previous timestep
// solution.previous.total.phase        // TODO : not needed ?
// solution.previous.partial(i).spin    // TODO : not needed ?
// solution.previous.partial(i).phase   // TODO : not needed ?

class Solution
{
public:
    template <class T> using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

public:
    Solution(const ScriptObject::Settings & settings);
    ~Solution();

    void clear();

private:
    Solution(const Solution &) = delete;
    Solution & operator=(const Solution &) = delete;

public:
    // TODO : turn timestep to unisgned
    // TODO : is it usefull since data[Index::Timestep] ?
    int timeStep;

    void save(double time, const std::vector<ResultPtr> & results);

    struct Current {
        void clear();
        alignas(16) aligned_vector<Eigen::Vector4d> spin;

        struct Total {
            void clear();
            alignas(16) Eigen::Vector4d spin;
            double phase;

            struct AtomLosses {
                void clear();
                double & operator()(int i) {switch(i) {case 1: return p_fromF1; case 2: return p_fromF2; default: throw;}}
            private:
                double p_fromF1;
                double p_fromF2;
            } atomLosses;
        } total;

        struct Partial {
            void clear();
        private:
            struct Internal {
                void clear();
                alignas(16) Eigen::Vector4d spin;
                double phase;
            };
            std::vector<Internal> p_internal;
        public:
            inline void resize(std::vector<Internal>::size_type newSize) {p_internal.resize(newSize);}
            inline Internal & operator()(int i) {return p_internal[i];}
        } partial;
    } current;

    struct Previous {
        void clear();
        struct Total {
            void clear();
            alignas(16) Eigen::Vector4d spin;
            double phase;
        } total;

        struct Partial {
            void clear();
        private:
            struct Internal {
                void clear();
                alignas(16) Eigen::Vector4d spin;
                double phase;
            };
            std::vector<Internal> p_internal;
        public:
            inline void resize(std::vector<Internal>::size_type newSize) {p_internal.resize(newSize);}
            inline Internal & operator()(int i) {return p_internal[i];}
        } partial;
    } previous;
};


#endif // SOLUTION_H
