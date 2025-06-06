#ifndef TGSCHEME_H
#define TGSCHEME_H

#include "Variance.h"     // Include the base Variance class
#include <vector>         // For std::vector
#include <random>         // For std::mt19937
#include <mutex>          // For thread safety (if needed)

// Declaration of the TGScheme class, derived from Variance
class TGScheme : public Variance {
public:
    // Constructor: initializes the psi range and grid size
    TGScheme(double psi_min = 0.001, double psi_max = 10.0, int grid_size = 5000);

    // Copy constructor
    TGScheme(const TGScheme& other);

    // Assignment operator
    TGScheme& operator=(const TGScheme& other);

    // Destructor: uses default implementation
    ~TGScheme() override = default;

    // Main method: computes the next variance value using the TG scheme
    double step(const double& Vt,
        const double& dt,
        const double& kappa,
        const double& theta,
        const double& epsilon,
        std::mt19937& rng) const override;

    // Clone method: creates a copy of the object
    TGScheme* clone() const override;

private:
    double psi_min_, psi_max_;    // Bounds of psi for the grid
    int grid_size_;               // Number of points in the grid

    std::vector<double> psi_grid_;  // Precomputed psi grid values
    std::vector<double> r_grid_;    // Corresponding r values for each psi

    // Helper methods for initialization and Newton's method
    void initialize_grid();
    double compute_r(double psi) const;
    static double objective_function(double r, double psi);
    static double objective_derivative(double r, double psi);
    static double phi(double x);
    static double Phi(double x);

    mutable std::mutex cache_mutex_;  // Thread safety (if needed)
};

#endif // TGSCHEME_H
