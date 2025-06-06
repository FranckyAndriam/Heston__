#include <iostream>        // For console output
#include <vector>          // For std::vector
#include <fstream>         // For file output
#include <cmath>           // For mathematical functions
#include <iomanip>         // For setting decimal precision

#include "TGScheme.h"             // Truncated Gaussian variance scheme
#include "BroadieKayaScheme.h"    // Broadie-Kaya scheme for (S,V) simulation
#include "PathSimulator2D.h"      // Simulates (S,V) paths
#include "Payoff.h"               // Option payoff calculation
#include "MonteCarlo.h"           // Monte Carlo pricing
#include "FourierPricer.h"        // Analytical Fourier pricer

int main() {
    // Fixed Heston model parameters
    double kappa = 2.0;
    double theta = 0.04;
    double epsilon = 0.5;
    double rho = -0.7;
    double r = 0.01;
    double gamma1 = 0.5;
    double gamma2 = 0.5;
    double psi_c = 1.5;
    double S0 = 100.0;
    double V0 = 0.04;
    double T = 1.0;
    double risk_free_rate = r;
    Call_Put option_type = Call_Put::Call;
    size_t steps = 100;

    // List of strikes to test
    std::vector<double> strike_list = { 70.0, 100.0, 140.0 };

    // Number of simulations to test
    std::vector<size_t> num_sim_list = { 10, 100, 1000, 10000, 50000, 100000 };

    // Time grid (fixed)
    std::vector<double> time_points(steps + 1);
    for (size_t i = 0; i <= steps; ++i)
        time_points[i] = T * i / steps;

    // Initialize the TG variance scheme
    TGScheme tg_scheme;
    const Variance* variance_scheme = &tg_scheme;

    // CSV file for exporting results
    std::ofstream file("TG_vs_Fourier_by_numSim_per_strike.csv");
    file << "Strike,NumSim,MC_Price,Fourier_Price,Absolute_Error,Relative_Error\n";
    file << std::fixed << std::setprecision(6);

    // Loops over strikes and number of simulations
    for (size_t i = 0; i < strike_list.size(); ++i) {
        double strike = strike_list[i];

        // Model and payoff fixed for the current strike
        BroadieKayaScheme heston_model(kappa, theta, epsilon, rho, r, gamma1, gamma2, variance_scheme);
        PathSimulator2D path_sim(time_points, S0, V0, heston_model);
        EuropeanOptionPayoff payoff(risk_free_rate, strike, option_type);

        // Exact price via Fourier (independent of num_sim)
        FourierPricer fourier(S0, V0, kappa, theta, epsilon, rho, T, strike);
        double fourier_price = fourier.computePrice(8000, 100.0);

        for (size_t j = 0; j < num_sim_list.size(); ++j) {
            size_t num_sim = num_sim_list[j];

            // Monte Carlo pricing with current num_sim
            MonteCarlo mc(num_sim, path_sim, payoff);
            double mc_price = mc.price();

            // Compute absolute and relative error
            double abs_err = std::abs(mc_price - fourier_price);
            double rel_err = abs_err / fourier_price;

            // Write results to the CSV file
            file << strike << "," << static_cast<unsigned long long>(num_sim) << ","
                << mc_price << "," << fourier_price << ","
                << abs_err << "," << rel_err << "\n";

            // Print results to console
            std::cout << "Strike=" << strike
                << " | NumSim=" << static_cast<unsigned long long>(num_sim)
                << " | MC=" << mc_price
                << " | Fourier=" << fourier_price
                << " | AbsErr=" << abs_err
                << " | RelErr=" << rel_err << "\n";
        }
    }

    file.close();
    std::cout << "\nResults saved in TG_vs_Fourier_by_numSim_per_strike.csv\n";
    return 0;
}
