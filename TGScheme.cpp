#define _USE_MATH_DEFINES
#include "TGScheme.h"
#include <cmath>
#include <algorithm>

TGScheme::TGScheme(double psi_min, double psi_max, int grid_size)
    : psi_min_(psi_min), psi_max_(psi_max), grid_size_(grid_size) {
    initialize_grid();
}

TGScheme::TGScheme(const TGScheme& other)
    : psi_min_(other.psi_min_), psi_max_(other.psi_max_), grid_size_(other.grid_size_),
    psi_grid_(other.psi_grid_), r_grid_(other.r_grid_) {}

TGScheme& TGScheme::operator=(const TGScheme& other) {
    if (this != &other) {
        psi_min_ = other.psi_min_;
        psi_max_ = other.psi_max_;
        grid_size_ = other.grid_size_;
        psi_grid_ = other.psi_grid_;
        r_grid_ = other.r_grid_;
    }
    return *this;
}

double TGScheme::phi(double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2 * M_PI);
}

double TGScheme::Phi(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
}

double TGScheme::objective_function(double r, double psi) {
    double ph = phi(r);
    double PH = Phi(r);
    return r * ph + (1 + r * r) * PH - (1 + psi) * std::pow(ph + r * PH, 2);
}

double TGScheme::objective_derivative(double r, double psi) {
    double ph = phi(r);
    double PH = Phi(r);
    double common_term = ph + r * PH;

    double derivative_common_term = -r * ph + PH;

    double derivative = ph + 2 * r * PH -
        2 * (1 + psi) * common_term * derivative_common_term;

    return derivative;
}

double TGScheme::compute_r(double psi) const {
    // Méthode de Newton robuste (selon Andersen)
    double r = (psi < 1.0) ? 1.0 : -4.0; // Initial guess selon Andersen
    const int max_iter = 100;
    const double tol = 1e-10;

    for (int iter = 0; iter < max_iter; ++iter) {
        double f = objective_function(r, psi);
        double df = objective_derivative(r, psi);

        if (std::abs(df) < 1e-8) df = (df < 0 ? -1 : 1) * 1e-8; // Protection contre division zéro

        double delta = f / df;
        r -= delta;

        if (std::abs(delta) < tol) break;

        // Clamping pour stabilité
        r = std::min(std::max(r, -10.0), 10.0);
    }

    return r;
}

void TGScheme::initialize_grid() {
    psi_grid_.resize(grid_size_);
    r_grid_.resize(grid_size_);

    double dpsi = (psi_max_ - psi_min_) / (grid_size_ - 1);

    for (int i = 0; i < grid_size_; ++i) {
        double psi = psi_min_ + i * dpsi;
        psi_grid_[i] = psi;
        r_grid_[i] = compute_r(psi);
    }
}

double TGScheme::step(const double& Vt,
    const double& dt,
    const double& kappa,
    const double& theta,
    const double& epsilon,
    std::mt19937& rng) const {
    double exp_kdt = std::exp(-kappa * dt);
    double m = theta + (Vt - theta) * exp_kdt;

    double s2 = Vt * epsilon * epsilon * exp_kdt * (1.0 - exp_kdt) / kappa +
        theta * epsilon * epsilon * std::pow(1.0 - exp_kdt, 2) / (2.0 * kappa);

    double psi = s2 / (m * m);

    // Interpolation linéaire de r
    double r;
    if (psi <= psi_min_) {
        r = r_grid_.front();
    }
    else if (psi >= psi_max_) {
        r = r_grid_.back();
    }
    else {
        double idx = (psi - psi_min_) / (psi_max_ - psi_min_) * (grid_size_ - 1);
        int idx_lower = static_cast<int>(idx);
        double weight = idx - idx_lower;
        r = r_grid_[idx_lower] * (1.0 - weight) + r_grid_[idx_lower + 1] * weight;
    }

    double ph = phi(r);
    double PH = Phi(r);
    double mu = (m * r) / (ph + r * PH);
    double sigma = m / (ph + r * PH);

    std::normal_distribution<double> normal(0.0, 1.0);
    double Z = normal(rng);

    return std::max(0.0, mu + sigma * Z);
}

TGScheme* TGScheme::clone() const {
    return new TGScheme(*this);
}