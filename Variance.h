#ifndef VARIANCE_H
#define VARIANCE_H

#include <random>  // For random number generator

// Abstract base class representing a variance process discretization scheme
class Variance {
public:
    // Default constructor
    Variance() = default;

    // Default copy constructor
    Variance(const Variance&) = default;

    // Default assignment operator
    Variance& operator=(const Variance&) = default;

    // Default destructor (virtual to allow proper cleanup of derived classes)
    virtual ~Variance() = default;

    // Pure virtual method: computes the next variance value V_{t+dt}
    virtual double step(
        const double& Vt,       // Current variance
        const double& dt,       // Time step
        const double& kappa,    // Mean reversion speed
        const double& theta,    // Long-term mean
        const double& epsilon,  // Volatility of variance
        std::mt19937& rng       // Random number generator
    ) const = 0;

    // Pure virtual method: creates a clone of the object (for polymorphic copying)
    virtual Variance* clone() const = 0;
};

#endif // VARIANCE_H
