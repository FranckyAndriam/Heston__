#ifndef QESCHEME_H
#define QESCHEME_H

#include "Variance.h"  // Include the abstract Variance class interface

// Declaration of the QEScheme class, derived from Variance
class QEScheme : public Variance {
public:
    // Constructor: initializes the psi_c threshold parameter
    QEScheme(const double& psi_c);

    // Copy constructor
    QEScheme(const QEScheme& other);

    // Assignment operator
    QEScheme& operator=(const QEScheme& other);

    // Destructor: uses the default implementation
    ~QEScheme() override = default;

    // Main method: computes the next variance value using the QE scheme
    double step(
        const double& Vt,       // Current variance
        const double& dt,       // Time step
        const double& kappa,    // Mean reversion speed
        const double& theta,    // Long-term mean
        const double& epsilon,  // Volatility of variance
        std::mt19937& rng       // Random number generator
    ) const override;

    // Clone method: creates a copy of the object
    QEScheme* clone() const override;

private:
    double psi_c;  // Threshold parameter to select between branches in the QE scheme
};

#endif // QESCHEME_H
