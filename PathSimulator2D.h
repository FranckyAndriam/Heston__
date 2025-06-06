#ifndef PATHSIMULATOR2D_H
#define PATHSIMULATOR2D_H

#include <vector>
#include <random>
#include "Model2D.h"

// PathSimulator2D.h
class PathSimulator2D {
public:
	// Constructeur pour simuler un chemin 2D avec un modèle donné. Ici model est un pointeur vers une instance de Model2D.
    PathSimulator2D(const std::vector<double>& time_points, double S0, double V0, const Model2D& model);
    PathSimulator2D(const PathSimulator2D& other);
    PathSimulator2D& operator=(const PathSimulator2D& other);

    ~PathSimulator2D();

    std::vector<double> path() const;           // path sous-jacent S
    std::vector<double> variance_path() const;  // path variance V
    std::vector<double> time_points() const { return _time_points; }

private:
    std::vector<double> _time_points;
    double _S0, _V0;
    Model2D* _model;
    mutable std::mt19937 _rng;
};

#endif // PATHSIMULATOR2D_H


