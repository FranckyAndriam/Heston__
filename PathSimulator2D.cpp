#include "PathSimulator2D.h"

// Constructor: initializes the time grid, initial values, and clones the model
PathSimulator2D::PathSimulator2D(const std::vector<double>& time_points, double S0, double V0, const Model2D& model)
    : _time_points(time_points), _S0(S0), _V0(V0), _model(model.clone()), _rng(std::random_device{}()) {}

// Copy constructor: performs deep copy by cloning the model
PathSimulator2D::PathSimulator2D(const PathSimulator2D& other)
    : _time_points(other._time_points), _S0(other._S0), _V0(other._V0), _model(other._model->clone()), _rng(std::random_device{}()) {}

// Assignment operator: performs deep copy of the time grid and model
PathSimulator2D& PathSimulator2D::operator=(const PathSimulator2D& other) {
    if (this != &other) {
        _time_points = other._time_points;
        _S0 = other._S0;
        _V0 = other._V0;
        delete _model;  // Clean up existing model
        _model = other._model->clone();  // Deep copy of the model
    }
    return *this;
}

// Destructor: deletes the cloned model
PathSimulator2D::~PathSimulator2D() {
    delete _model;
}

// Generates the path of the underlying asset price S over the time grid
std::vector<double> PathSimulator2D::path() const {
    std::vector<double> S_path(_time_points.size());
    std::vector<double> V_path(_time_points.size());
    S_path[0] = _S0;
    V_path[0] = _V0;

    // Loop over time steps to generate the trajectory
    for (size_t i = 1; i < _time_points.size(); ++i) {
        double dt = _time_points[i] - _time_points[i - 1];
        auto res = _model->step(S_path[i - 1], V_path[i - 1], dt, _rng);
        S_path[i] = res.first;
        V_path[i] = res.second;
    }
    return S_path;
}

// Generates the path of the variance V over the time grid
std::vector<double> PathSimulator2D::variance_path() const {
    std::vector<double> S_path(_time_points.size());
    std::vector<double> V_path(_time_points.size());
    S_path[0] = _S0;
    V_path[0] = _V0;

    // Loop over time steps to generate the trajectory
    for (size_t i = 1; i < _time_points.size(); ++i) {
        double dt = _time_points[i] - _time_points[i - 1];
        auto res = _model->step(S_path[i - 1], V_path[i - 1], dt, _rng);
        S_path[i] = res.first;
        V_path[i] = res.second;
    }
    return V_path;
}
