#include <iostream>
#include <vector>
#include <cmath>

// Inclusions des headers de ton projet :
#include "QEScheme.h"
#include "TGScheme.h"
#include "BroadieKayaScheme.h"
#include "PathSimulator2D.h"
#include "Payoff.h"
#include "MonteCarlo.h"

int main() {
    // === Paramètres du modèle Heston ===
    double kappa = 2.0;
    double theta = 0.04;
    double epsilon = 0.5;
    double rho = -0.7;
    double r = 0.05;
    double gamma1 = 0.5;
    double gamma2 = 0.5;
    double psi_c = 1.5;
    double S0 = 110.0;
    double V0 = 0.04;

    // === Paramètres option européenne ===
    double T = 1.0;
    double strike = 100.0;
    Call_Put option_type = Call_Put::Call;
    double risk_free_rate = r;

    // === Discrétisation temporelle ===
    size_t steps = 252;
    std::vector<double> time_points(steps + 1);
    for (size_t i = 0; i <= steps; ++i) {
        time_points[i] = T * i / steps;
    }

    // === Schéma QE ou TG ===
    bool use_QE = false; // Passe à false pour tester TG

    QEScheme qe_scheme(psi_c);
    TGScheme tg_scheme(kappa, theta, epsilon, T / steps);

    const Variance* variance_scheme = use_QE
        ? static_cast<const Variance*>(&qe_scheme)
        : static_cast<const Variance*>(&tg_scheme);

    // === Broadie-Kaya Heston Model2D ===
    BroadieKayaScheme heston_model(
        kappa, theta, epsilon, rho, r,
        gamma1, gamma2, variance_scheme
    );

    // === Simulateur de trajectoires ===
    PathSimulator2D path_sim(time_points, S0, V0, heston_model);

    // === Payoff européen ===
    EuropeanOptionPayoff payoff(risk_free_rate, strike, option_type);

    // === Pricing Monte Carlo ===
    size_t num_sim = 10000;
    MonteCarlo mc(num_sim, path_sim, payoff);

    double price = mc.price();

    std::cout << "European Call price (Heston MC, "
        << (use_QE ? "QE" : "TG") << " scheme): " << price << std::endl;

    return 0;
}
