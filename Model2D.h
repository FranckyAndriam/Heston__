#ifndef MODEL2D_H
#define MODEL2D_H

#include <random>     // For random number generation
#include <utility>    // For std::pair

// Interface for any 2D model (e.g., Heston, SABR)
class Model2D {
public:
    // Virtual destructor for proper cleanup of derived classes
    virtual ~Model2D() = default;

    // Pure virtual function to simulate the next (S, V) state given the current state (S, V) and time step dt
    virtual std::pair<double, double> step(
        double St,             // Current asset price S_t
        double Vt,             // Current variance V_t
        double dt,             // Time increment
        std::mt19937& rng      // Random number generator
    ) const = 0;

    // Pure virtual function to clone the object polymorphically
    virtual Model2D* clone() const = 0;
};

#endif // MODEL2D_H
