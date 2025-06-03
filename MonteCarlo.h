#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <cstddef>

// Forward declarations to avoid full include
class PathSimulator;
class PathSimulator2D;
class Payoff;

class MonteCarlo
{
public:
    // Pour PathSimulator2D (Heston, etc.)
    MonteCarlo(const size_t num_sims, const PathSimulator2D& path_sim2d, const Payoff& payoff);

    double price() const;

private:
    size_t _number_simulations;
    const PathSimulator* _path_simulator = nullptr;
    const PathSimulator2D* _path_simulator2d = nullptr;
    const Payoff* _payoff = nullptr;
    bool _use2d = false;
};

#endif // MONTECARLO_H

