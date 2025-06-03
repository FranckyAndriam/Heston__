#include "MonteCarlo.h"
#include "PathSimulator2D.h"
#include "Payoff.h"


MonteCarlo::MonteCarlo(const size_t num_sims, const PathSimulator2D& path_sim2d, const Payoff& payoff)
    : _number_simulations(num_sims), _path_simulator2d(&path_sim2d), _payoff(&payoff), _use2d(true)
{
}

double MonteCarlo::price() const
{
    double price = 0.0;
    for (size_t sim_idx = 0; sim_idx < _number_simulations; sim_idx++) {
        if (_use2d) {
            price += _payoff->discounted_payoff(_path_simulator2d->path(), _path_simulator2d->time_points());
        }
    }
    price /= static_cast<double>(_number_simulations);
    return price;
}
