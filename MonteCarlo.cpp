#include "MonteCarlo.h"
#include "PathSimulator2D.h"
#include "Payoff.h"

// Constructor for the MonteCarlo class
// Initializes the simulation parameters: number of simulations, path simulator, and payoff
MonteCarlo::MonteCarlo(const size_t num_sims, const PathSimulator2D& path_sim2d, const Payoff& payoff)
    : _number_simulations(num_sims), _path_simulator2d(&path_sim2d), _payoff(&payoff), _use2d(true)
{
}

// Main method to compute the option price via Monte Carlo
double MonteCarlo::price() const
{
    double price = 0.0;

    // Loop over the number of simulations
    for (size_t sim_idx = 0; sim_idx < _number_simulations; sim_idx++) {
        if (_use2d) {
            // For each simulation, compute the discounted payoff using the path simulator and payoff
            price += _payoff->discounted_payoff(_path_simulator2d->path(), _path_simulator2d->time_points());
        }
    }

    // Calculate the average price over all simulations
    price /= static_cast<double>(_number_simulations);
    return price;
}
