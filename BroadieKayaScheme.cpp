#include "BroadieKayaScheme.h"
#include <cmath>

BroadieKayaScheme::BroadieKayaScheme(double kappa, double theta, double epsilon, double rho, double r,
    double gamma1, double gamma2, const Variance* variance_scheme)
    : kappa(kappa), theta(theta), epsilon(epsilon), rho(r), r(r),
    gamma1(gamma1), gamma2(gamma2), _variance_scheme(variance_scheme)
{
    const double sum = gamma1 + gamma2;
    if (std::abs(sum - 1.0) > 1e-12) {
        this->gamma1 /= sum;
        this->gamma2 /= sum;
    }
}

std::pair<double, double> BroadieKayaScheme::step(double St, double Vt, double dt, std::mt19937& rng) const
{
    // On utilise ici le schéma polymorphe (QE ou TG, etc)
    double V_next = _variance_scheme->step(Vt, dt, kappa, theta, epsilon, rng);
    double V_bar = gamma1 * Vt + gamma2 * V_next;

    double mu_X = std::log(St)
        + r * dt
        + rho / epsilon * (V_next - Vt - kappa * theta * dt)
        + (kappa * rho / epsilon - 0.5) * dt * V_bar;

    double sigma_X = std::sqrt((1.0 - rho * rho) * dt * V_bar);

    std::normal_distribution<double> norm(0.0, 1.0);
    double Zx = norm(rng);

    double X_next = mu_X + sigma_X * Zx;
    double S_next = std::exp(X_next);

    return { S_next, V_next };
}

BroadieKayaScheme* BroadieKayaScheme::clone() const {
    // Attention: ici, on ne copie pas le schéma de variance mais juste le pointeur !
    // Pour une copie profonde, il faudrait stocker un unique_ptr<Variance> (et cloner dans clone())
    return new BroadieKayaScheme(kappa, theta, epsilon, rho, r, gamma1, gamma2, _variance_scheme);
}
