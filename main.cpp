#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip> // pour précision CSV

// Inclusions headers projet
#include "QEScheme.h"
#include "TGScheme.h"
#include "BroadieKayaScheme.h"
#include "PathSimulator2D.h"
#include "Payoff.h"
#include "MonteCarlo.h"
#include "FourierPricer.h"

int main() {
    // === Paramètres fixes du modèle Heston ===
    double kappa = 2.0;
    double theta = 0.04;
    double epsilon = 0.5;
    double rho = -0.7;
    double r = 0.01;
    double gamma1 = 0.5;
    double gamma2 = 0.5;
    double psi_c = 1.5;
    double risk_free_rate = r;
    Call_Put option_type = Call_Put::Call;

    // === Grilles d'exploration ===
    std::vector<double> strikes = { 90, 100, 110};      // K
    std::vector<double> maturities = { 0.5, 1.0, 2.0 };     // T
    std::vector<double> S0_list = { 90, 100, 110 };      // S0
    std::vector<double> V0_list = { 0.02, 0.04, 0.06 };  // V0

    // === Setup schéma QE (toujours le même objet) ===
    QEScheme qe_scheme(psi_c);
    const Variance* variance_scheme = &qe_scheme;

    // === Fichier CSV au nom explicite ===
    std::ofstream file("comparaison_QE_vs_Fourier_by_K_T_S0_V0.csv");
    file << "S0,V0,Strike,Maturity,MC_Price,Fourier_Price,Absolute_Error,Relative_Error\n";
    file << std::fixed << std::setprecision(6);

    for (double S0 : S0_list) {
        for (double V0 : V0_list) {
            for (double K : strikes) {
                for (double T : maturities) {
                    // --- Redéfinit la grille de temps pour cette maturité
                    size_t steps = static_cast<size_t>(252 * T); // pour dt ~ 1/252
                    std::vector<double> time_points(steps + 1);
                    for (size_t i = 0; i <= steps; ++i) {
                        time_points[i] = T * i / steps;
                    }

                    // --- Modèle Broadie-Kaya
                    BroadieKayaScheme heston_model(
                        kappa, theta, epsilon, rho, r,
                        gamma1, gamma2, variance_scheme
                    );
                    PathSimulator2D path_sim(time_points, S0, V0, heston_model);

                    // --- Payoff européen
                    EuropeanOptionPayoff payoff(risk_free_rate, K, option_type);

                    // --- Pricing Monte Carlo (QE)
                    size_t num_sim = 10000;
                    MonteCarlo mc(num_sim, path_sim, payoff);
                    double mc_price = mc.price();

                    // --- Pricing exact Fourier
                    FourierPricer analytic_pricer(
                        S0, V0, kappa, theta, epsilon, rho, T, K
                    );
                    double analytic_price = analytic_pricer.computePrice(8000, 100.0);

                    // --- Erreur
                    double abs_error = std::abs(mc_price - analytic_price);
                    double rel_error = (analytic_price != 0.0) ? abs_error / analytic_price : 0.0;

                    // --- Écriture CSV
                    file << S0 << "," << V0 << "," << K << "," << T << ","
                        << mc_price << "," << analytic_price << ","
                        << abs_error << "," << rel_error << "\n";

                    // Affichage pour suivi
                    std::cout << "S0=" << S0 << " V0=" << V0
                        << " K=" << K << " T=" << T
                        << "  MC=" << mc_price
                        << "  Fourier=" << analytic_price
                        << "  AbsErr=" << abs_error
                        << "  RelErr=" << rel_error << std::endl;
                }
            }
        }
    }

    file.close();
    std::cout << "\nRésultats enregistrés dans comparaison_QE_vs_Fourier_by_K_T_S0_V0.csv\n";
    return 0;
}
