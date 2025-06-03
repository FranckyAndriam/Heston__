#ifndef VARIANCE_H
#define VARIANCE_H

#include <random>

class Variance {
public:
    Variance() = default;
    Variance(const Variance&) = default;
    Variance& operator=(const Variance&) = default;
    virtual ~Variance() = default;

    virtual double step(
        const double& Vt,
        const double& dt,
        const double& kappa,
        const double& theta,
        const double& epsilon,
        std::mt19937& rng
    ) const = 0;

    virtual Variance* clone() const = 0;
};

#endif // VARIANCE_H


