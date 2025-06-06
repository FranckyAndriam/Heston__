#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <cstddef>  // For size_t

// Forward declarations to avoid full includes and circular dependencies
class PathSimulator;
class PathSimulator2D;
class Payoff;

// Declaration of the MonteCarlo class
class MonteCarlo
{
public:
    // Constructor for PathSimulator2D (used in models like Heston)
    MonteCarlo(const size_t num_sims, const PathSimulator2D& path_sim2d, const Payoff& payoff);

    // Main method to compute the option price via Monte Carlo
    double price() const;

private:
    size_t _number_simulations;  // Number of simulations to run

    // Pointers to path simulators and payoff (non-owning pointers)
    const PathSimulator* _path_simulator = nullptr;
    const PathSimulator2D* _path_simulator2d = nullptr;

    const Payoff* _payoff = nullptr;  // Pointer to the payoff function

    bool _use2d = false;  // Flag to indicate whether to use 2D simulation
};

#endif // MONTECARLO_H
