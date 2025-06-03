#define _USE_MATH_DEFINES
#include "TGScheme.h"
#include <cmath>
#include <algorithm>

TGScheme::TGScheme(double kappa_, double theta_, double epsilon_, double dt_)
    : kappa(kappa_), theta(theta_), epsilon(epsilon_), dt(dt_), gen(std::random_device{}()), norm(0.0, 1.0) {}

TGScheme& TGScheme::operator=(const TGScheme& other) {
    if (this != &other) {
        kappa = other.kappa;
        theta = other.theta;
        epsilon = other.epsilon;
        dt = other.dt;
        gen = other.gen;
        norm = other.norm;
    }
    return *this;
}

TGScheme::~TGScheme() {}

double TGScheme::step(
    const double& Vt,
    const double& /*dt_in*/,
    const double& /*kappa_in*/,
    const double& /*theta_in*/,
    const double& /*epsilon_in*/,
    std::mt19937& /*rng*/
) const {
    return simulateNextV(Vt);
}

TGScheme* TGScheme::clone() const {
    return new TGScheme(*this);
}

double TGScheme::simulateNextV(double Vt) const {
    double m = computeMean(Vt);
    double s2 = computeVariance(Vt);
    double psi = computePsi(m, s2);
    double r = lookup_r(psi);

    double phi_r = std::exp(-0.5 * r * r) / std::sqrt(2.0 * M_PI);
    double Phi_r = 0.5 * (1.0 + std::erf(r / std::sqrt(2.0)));

    double f_mu = r / (phi_r + r * Phi_r);
    double f_sigma = std::sqrt(psi) * f_mu;

    double mu = f_mu * m;
    double sigma = f_sigma * std::sqrt(s2);

    double Z = norm(gen);
    return std::max(0.0, mu + sigma * Z);
}

double TGScheme::computeMean(double Vt) const {
    return theta + (Vt - theta) * std::exp(-kappa * dt);
}

double TGScheme::computeVariance(double Vt) const {
    double ekt = std::exp(-kappa * dt);
    double term1 = Vt * epsilon * epsilon * ekt * (1.0 - ekt) / kappa;
    double term2 = theta * epsilon * epsilon * std::pow(1.0 - ekt, 2) / (2.0 * kappa);
    return term1 + term2;
}

double TGScheme::computePsi(double m, double s2) const {
    return s2 / (m * m);
}

double TGScheme::lookup_r(double psi) const {
    return std::sqrt(2.0 * std::log(1.0 + psi));
}
