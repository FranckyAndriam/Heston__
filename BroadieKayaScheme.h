#ifndef BROADIEKAYASCHEME_H
#define BROADIEKAYASCHEME_H

#include "Variance.h"
#include "Model2D.h"
#include <random>
#include <utility>

class BroadieKayaScheme : public Model2D {
public:
    BroadieKayaScheme(double kappa, double theta, double epsilon, double rho, double r,
        double gamma1, double gamma2, const Variance* variance_scheme);

    std::pair<double, double> step(double St, double Vt, double dt, std::mt19937& rng) const override;
    BroadieKayaScheme* clone() const override;

private:
    double kappa, theta, epsilon, rho, r, gamma1, gamma2;
    const Variance* _variance_scheme; // Non-owning pointer
};

#endif // BROADIEKAYASCHEME_H

