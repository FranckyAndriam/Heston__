#ifndef TGSCHEME_H
#define TGSCHEME_H

#include "Variance.h"
#include <random>

class TGScheme : public Variance {
public:
    TGScheme(double kappa_, double theta_, double epsilon_, double dt_);
    TGScheme(const TGScheme& other) = default;
    TGScheme& operator=(const TGScheme& other);
    ~TGScheme() override;

    double step(
        const double& Vt,
        const double& dt,
        const double& kappa,
        const double& theta,
        const double& epsilon,
        std::mt19937& rng
    ) const override;

    TGScheme* clone() const override;

private:
    double simulateNextV(double Vt) const;
    double computeMean(double Vt) const;
    double computeVariance(double Vt) const;
    double computePsi(double m, double s2) const;
    double lookup_r(double psi) const;

    mutable std::mt19937 gen;
    mutable std::normal_distribution<double> norm;

    double kappa, theta, epsilon, dt;
};

#endif // TGSCHEME_H

