#ifndef PAYOFF_H
#define PAYOFF_H

#include<vector>


enum class Call_Put
{
	Call,
	Put
};


class Payoff
{
public:
	Payoff(const double& risk_free_rate);
	virtual double discounted_payoff(const std::vector<double>& underlying_path, const std::vector<double>& time_points) const = 0;

	virtual Payoff* clone() const = 0;

protected:
	double _risk_free_rate;

};

class EuropeanOptionPayoff : public Payoff
{
public:
	EuropeanOptionPayoff(const double& risk_free_rate, const double& strike, const Call_Put& option_type);

	double discounted_payoff(const std::vector<double>& underlying_path, const std::vector<double>& time_points) const override;

	EuropeanOptionPayoff* clone() const override;

private:

	double _strike;
	Call_Put _option_type;
};





#endif

