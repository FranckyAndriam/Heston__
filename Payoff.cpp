#include "Payoff.h"

#include <algorithm>

Payoff::Payoff(const double& risk_free_rate)
	: _risk_free_rate(risk_free_rate)
{
}

EuropeanOptionPayoff::EuropeanOptionPayoff(const double& risk_free_rate, const double& strike, const Call_Put& option_type)
	: Payoff(risk_free_rate), _strike(strike), _option_type(option_type)
{
}

double EuropeanOptionPayoff::discounted_payoff(const std::vector<double>& underlying_path, const std::vector<double>& time_points) const
{
	size_t final_index = time_points.size() - 1;
	double underlying_value = underlying_path[final_index];
	double maturity = time_points[final_index];

	// return max( +1/-1 (S_T - K) , 0. )
	double multiplier = (_option_type == Call_Put::Call) ? 1. : 0.;

	double undisc_payoff = std::max(multiplier * (underlying_value - _strike), 0.);
	double disc_payoff = exp(-_risk_free_rate * maturity) * undisc_payoff;

	return disc_payoff;
}

EuropeanOptionPayoff* EuropeanOptionPayoff::clone() const
{
	return new EuropeanOptionPayoff(*this);
}


