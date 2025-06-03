#include "QEScheme.h"
#include <cmath>
#include <random>

QEScheme::QEScheme(const double& psi_c)
    : psi_c(psi_c) {}
QEScheme::QEScheme(const QEScheme& other)
    : psi_c(other.psi_c) {}
QEScheme& QEScheme::operator=(const QEScheme& other) {
    if (this != &other) {
        psi_c = other.psi_c;
    }
    return *this;
}

double QEScheme::step(
    const double& Vt,
    const double& dt,
    const double& kappa,
    const double& theta,
    const double& epsilon,
    std::mt19937& rng
) const {
    const double e_kd = std::exp(-kappa * dt);
    const double m = theta + (Vt - theta) * e_kd;
    const double var1 = Vt * epsilon * epsilon * e_kd / kappa * (1.0 - e_kd);
    const double var2 = theta * epsilon * epsilon / (2.0 * kappa) * (1.0 - e_kd) * (1.0 - e_kd);
    const double s2 = var1 + var2;
    const double psi = s2 / (m * m);

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double u = uniform(rng);

    if (psi <= psi_c) {
        const double inv_psi = 1.0 / psi;
        const double term = 2.0 * inv_psi - 1.0;
        const double b2 = term + std::sqrt(2.0 * inv_psi * term);
        const double b = std::sqrt(b2);
        const double a = m / (1.0 + b2);

        std::normal_distribution<double> normal(0.0, 1.0);
        const double z = normal(rng);
        return a * (b + z) * (b + z);
    }
    else {
        const double p = (psi - 1.0) / (psi + 1.0);
        const double beta = 2.0 / (m * (psi + 1.0));
        if (u <= p) {
            return 0.0;
        }
        return std::log((1.0 - p) / (1.0 - u)) / beta;
    }
}

QEScheme* QEScheme::clone() const {
    return new QEScheme(*this);
}
