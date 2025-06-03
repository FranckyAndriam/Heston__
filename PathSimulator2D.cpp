#include "PathSimulator2D.h"

PathSimulator2D::PathSimulator2D(const std::vector<double>& time_points, double S0, double V0, const Model2D& model)
    : _time_points(time_points), _S0(S0), _V0(V0), _model(model.clone()), _rng(std::random_device{}()) {}

PathSimulator2D::PathSimulator2D(const PathSimulator2D& other)
    : _time_points(other._time_points), _S0(other._S0), _V0(other._V0), _model(other._model->clone()), _rng(std::random_device{}()) {}

PathSimulator2D& PathSimulator2D::operator=(const PathSimulator2D& other) {
    if (this != &other) {
        _time_points = other._time_points;
        _S0 = other._S0;
        _V0 = other._V0;
        delete _model;
        _model = other._model->clone();
    }
    return *this;
}

PathSimulator2D::~PathSimulator2D() {
    delete _model;
}

std::vector<double> PathSimulator2D::path() const {
    std::vector<double> S_path(_time_points.size());
    std::vector<double> V_path(_time_points.size());
    S_path[0] = _S0;
    V_path[0] = _V0;

    for (size_t i = 1; i < _time_points.size(); ++i) {
        double dt = _time_points[i] - _time_points[i - 1];
        auto res = _model->step(S_path[i - 1], V_path[i - 1], dt, _rng);
        S_path[i] = res.first;
        V_path[i] = res.second;
    }
    return S_path;
}

std::vector<double> PathSimulator2D::variance_path() const {
    std::vector<double> S_path(_time_points.size());
    std::vector<double> V_path(_time_points.size());
    S_path[0] = _S0;
    V_path[0] = _V0;

    for (size_t i = 1; i < _time_points.size(); ++i) {
        double dt = _time_points[i] - _time_points[i - 1];
        auto res = _model->step(S_path[i - 1], V_path[i - 1], dt, _rng);
        S_path[i] = res.first;
        V_path[i] = res.second;
    }
    return V_path;
}
