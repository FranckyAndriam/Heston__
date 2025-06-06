#include "Payoff.h"

#include <algorithm>  // For std::max

// Constructor for the base Payoff class
Payoff::Payoff(const double& risk_free_rate)
	: _risk_free_rate(risk_free_rate)
{
}

// Constructor for the EuropeanOptionPayoff class
// Inherits from Payoff and stores the strike and option type (Call or Put)
EuropeanOptionPayoff::EuropeanOptionPayoff(const double& risk_free_rate, const double& strike, const Call_Put& option_type)
	: Payoff(risk_free_rate), _strike(strike), _option_type(option_type)
{
}

// Computes the discounted payoff of the European option given the simulated underlying path
double EuropeanOptionPayoff::discounted_payoff(const std::vector<double>& underlying_path, const std::vector<double>& time_points) const
{
	size_t final_index = time_points.size() - 1;  // Index of the final time point
	double underlying_value = underlying_path[final_index];  // Asset price at maturity
	double maturity = time_points[final_index];  // Time to maturity

	// Calculate the payoff depending on the option type
	// multiplier is 1 for Call, 0 for Put (Note: may need adjustment for Puts)
	double multiplier = (_option_type == Call_Put::Call) ? 1. : 0.;

	// Calculate the non-discounted payoff (e.g., max(S_T - K, 0))
	double undisc_payoff = std::max(multiplier * (underlying_value - _strike), 0.);

	// Discount the payoff using the risk-free rate
	double disc_payoff = exp(-_risk_free_rate * maturity) * undisc_payoff;

	return disc_payoff;
}

// Clone method: returns a new copy of the EuropeanOptionPayoff object
EuropeanOptionPayoff* EuropeanOptionPayoff::clone() const
{
	return new EuropeanOptionPayoff(*this);
}
