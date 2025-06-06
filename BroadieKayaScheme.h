#ifndef BROADIEKAYASCHEME_H
#define BROADIEKAYASCHEME_H

// Include the Variance scheme interface
#include "Variance.h"

// Include the base class Model2D, representing a 2D model (S and V)
#include "Model2D.h"

// Include standard library headers for random number generation and utility pair
#include <random>
#include <utility>

// BroadieKayaScheme class declaration, derived from Model2D
class BroadieKayaScheme : public Model2D {
public:
    // Constructor: takes Heston model parameters, weights, and a pointer to a Variance scheme
    BroadieKayaScheme(double kappa, double theta, double epsilon, double rho, double r,
        double gamma1, double gamma2, const Variance* variance_scheme);

    // Main method to simulate one step of (S, V)
    std::pair<double, double> step(double St, double Vt, double dt, std::mt19937& rng) const override;

    // Clone method to duplicate this object polymorphically
    BroadieKayaScheme* clone() const override;

private:
    // Heston model parameters
    double kappa, theta, epsilon, rho, r;

    // Weights for variance discretization (gamma1 and gamma2)
    double gamma1, gamma2;

    // Pointer to a Variance scheme (QE or TG).
    const Variance* _variance_scheme;
};

#endif // BROADIEKAYASCHEME_H
