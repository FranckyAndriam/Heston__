#include "QEScheme.h"
#include <cmath>      // For math functions
#include <random>     // For random number generation

// Constructor for QEScheme: initializes psi_c threshold parameter
QEScheme::QEScheme(const double& psi_c)
    : psi_c(psi_c) {}

// Copy constructor for QEScheme
QEScheme::QEScheme(const QEScheme& other)
    : psi_c(other.psi_c) {}

// Assignment operator for QEScheme
QEScheme& QEScheme::operator=(const QEScheme& other) {
    if (this != &other) {
        psi_c = other.psi_c;
    }
    return *this;
}

// Main method to simulate the next V_{t+dt} using the QE scheme
double QEScheme::step(
    const double& Vt,       // Current variance
    const double& dt,       // Time increment
    const double& kappa,    // Mean reversion speed
    const double& theta,    // Long-term mean
    const double& epsilon,  // Volatility of variance
    std::mt19937& rng       // Random number generator
) const {
    // Compute e^{-kappa * dt}
    const double e_kd = std::exp(-kappa * dt);

    // Compute the conditional mean m
    const double m = theta + (Vt - theta) * e_kd;

    // Compute the conditional variance s²
    const double var1 = Vt * epsilon * epsilon * e_kd / kappa * (1.0 - e_kd);
    const double var2 = theta * epsilon * epsilon / (2.0 * kappa) * (1.0 - e_kd) * (1.0 - e_kd);
    const double s2 = var1 + var2;

    // Calculate psi
    const double psi = s2 / (m * m);

    // Generate a uniform random number for the branching
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double u = uniform(rng);

    // If psi <= psi_c, use the non-central chi-squared approximation
    if (psi <= psi_c) {
        const double inv_psi = 1.0 / psi;
        const double term = 2.0 * inv_psi - 1.0;
        const double b2 = term + std::sqrt(2.0 * inv_psi * term);
        const double b = std::sqrt(b2);
        const double a = m / (1.0 + b2);

        // Draw a standard normal variable
        std::normal_distribution<double> normal(0.0, 1.0);
        const double z = normal(rng);

        // Return the next variance value
        return a * (b + z) * (b + z);
    }
    else {
        // Otherwise, use an exponential approximation
        const double p = (psi - 1.0) / (psi + 1.0);
        const double beta = 2.0 / (m * (psi + 1.0));

        // With probability p, variance collapses to 0
        if (u <= p) {
            return 0.0;
        }

        // Otherwise, return the exponential approximation
        return std::log((1.0 - p) / (1.0 - u)) / beta;
    }
}

// Clone method to duplicate a QEScheme object
QEScheme* QEScheme::clone() const {
    return new QEScheme(*this);
}
