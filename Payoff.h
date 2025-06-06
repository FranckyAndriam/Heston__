#ifndef PAYOFF_H
#define PAYOFF_H

#include <vector>  // For std::vector

// Enumeration to distinguish between Call and Put options
enum class Call_Put
{
	Call,
	Put
};

// Abstract base class for option payoffs
class Payoff
{
public:
	// Constructor: initializes the risk-free rate
	Payoff(const double& risk_free_rate);

	// Pure virtual method: calculates the discounted payoff given the underlying path and time points
	virtual double discounted_payoff(const std::vector<double>& underlying_path, const std::vector<double>& time_points) const = 0;

	// Pure virtual method: returns a copy of the Payoff object (polymorphic cloning)
	virtual Payoff* clone() const = 0;

protected:
	double _risk_free_rate;  // Risk-free interest rate for discounting
};

// Concrete class for European Call/Put option payoff
class EuropeanOptionPayoff : public Payoff
{
public:
	// Constructor: initializes the risk-free rate, strike, and option type
	EuropeanOptionPayoff(const double& risk_free_rate, const double& strike, const Call_Put& option_type);

	// Implements the discounted payoff calculation for European options
	double discounted_payoff(const std::vector<double>& underlying_path, const std::vector<double>& time_points) const override;

	// Implements the clone method
	EuropeanOptionPayoff* clone() const override;

private:
	double _strike;             // Option strike price
	Call_Put _option_type;      // Call or Put
};

#endif
