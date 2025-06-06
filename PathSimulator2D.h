#ifndef PATHSIMULATOR2D_H
#define PATHSIMULATOR2D_H

#include <vector>        // For using std::vector
#include <random>        // For random number generation
#include "Model2D.h"     // Include the abstract Model2D interface

// Declaration of the PathSimulator2D class
class PathSimulator2D {
public:
    // Constructor to simulate a 2D path with a given model.
    // 'model' is passed as a reference but cloned internally.
    PathSimulator2D(const std::vector<double>& time_points, double S0, double V0, const Model2D& model);

    // Copy constructor: deep copies the model
    PathSimulator2D(const PathSimulator2D& other);

    // Assignment operator: deep copies the model
    PathSimulator2D& operator=(const PathSimulator2D& other);

    // Destructor: cleans up the cloned model
    ~PathSimulator2D();

    // Returns the simulated path of the underlying asset S
    std::vector<double> path() const;

    // Returns the simulated path of the variance V
    std::vector<double> variance_path() const;

    // Returns the time grid used for the simulation
    std::vector<double> time_points() const { return _time_points; }

private:
    std::vector<double> _time_points;   // Time grid for simulation
    double _S0, _V0;                    // Initial values of S and V
    Model2D* _model;                    // Pointer to the model (cloned instance)
    mutable std::mt19937 _rng;          // Random number generator (mutable to allow const functions to use it)
};

#endif // PATHSIMULATOR2D_H
