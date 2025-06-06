#include "BroadieKayaScheme.h"
#include <cmath>

// Constructor for the BroadieKayaScheme class
// Takes in the Heston model parameters and the weights gamma1/gamma2 for the weighted variance.
BroadieKayaScheme::BroadieKayaScheme(double kappa, double theta, double epsilon, double rho, double r,
    double gamma1, double gamma2, const Variance* variance_scheme)
    : kappa(kappa), theta(theta), epsilon(epsilon), rho(r), r(r),
    gamma1(gamma1), gamma2(gamma2), _variance_scheme(variance_scheme)
{
    // Check the consistency of the weights gamma1 and gamma2
    const double sum = gamma1 + gamma2;
    if (std::abs(sum - 1.0) > 1e-12) {
        // Normalize the weights if their sum differs slightly from 1
        this->gamma1 /= sum;
        this->gamma2 /= sum;
    }
}

// Main method to simulate one step of (S, V) using the Broadie-Kaya scheme
std::pair<double, double> BroadieKayaScheme::step(double St, double Vt, double dt, std::mt19937& rng) const
{
    // Use the polymorphic variance scheme (QE or TG) to simulate V_{t+dt}
    double V_next = _variance_scheme->step(Vt, dt, kappa, theta, epsilon, rng);

    // Compute the weighted variance for drift and diffusion
    double V_bar = gamma1 * Vt + gamma2 * V_next;

    // Compute the drift (mu_X) of the log(S) process
    double mu_X = std::log(St)
        + r * dt // risk-neutral drift term
        + rho / epsilon * (V_next - Vt - kappa * theta * dt) // covariance correction
        + (kappa * rho / epsilon - 0.5) * dt * V_bar; // volatility adjustment

    // Compute the volatility (sigma_X) of the log(S) process
    double sigma_X = std::sqrt((1.0 - rho * rho) * dt * V_bar);

    // Generate a standard normal random variable for the process simulation
    std::normal_distribution<double> norm(0.0, 1.0);
    double Zx = norm(rng);

    // Simulate the next log(S)
    double X_next = mu_X + sigma_X * Zx;

    // Convert back to the asset price S_{t+dt}
    double S_next = std::exp(X_next);

    // Return the pair (S_{t+dt}, V_{t+dt})
    return { S_next, V_next };
}

// Clone method to duplicate a BroadieKayaScheme object
BroadieKayaScheme* BroadieKayaScheme::clone() const {
    return new BroadieKayaScheme(kappa, theta, epsilon, rho, r, gamma1, gamma2, _variance_scheme);
}
