#define _USE_MATH_DEFINES
#include "FourierPricer.h"
#include <cmath>
#include <complex>
#include <iostream>

// Define the imaginary unit i = sqrt(-1)
const std::complex<double> i(0.0, 1.0);

// Constructor: stores the Heston model parameters
FourierPricer::FourierPricer(double X0, double V0, double kappa, double theta,
    double epsilon, double rho, double T, double K)
    : X0(X0), V0(V0), kappa(kappa), theta(theta),
    epsilon(epsilon), rho(rho), T(T), K(K)
{
    // Calculate kappa_hat = kappa - rho * epsilon / 2
    kappa_hat = kappa - rho * epsilon / 2.0;
}

// Computes ξ = sqrt(k² ε² (1 - ρ²) + 2 i k ε ρ κ̂ + κ̂² + ε²/4)
std::complex<double> FourierPricer::xi(std::complex<double> k) {
    return sqrt(pow(k * epsilon, 2.0) * (1.0 - rho * rho)
        + 2.0 * i * k * epsilon * rho * kappa_hat
        + pow(kappa_hat, 2.0)
        + pow(epsilon, 2.0) / 4.0);
}

// Computes d_plus = ξ − (i k ε ρ + κ̂)
std::complex<double> FourierPricer::d_plus(std::complex<double> k) {
    return xi(k) - (i * k * rho * epsilon + kappa_hat);
}

// Computes d_minus = ξ + (i k ε ρ + κ̂)
std::complex<double> FourierPricer::d_minus(std::complex<double> k) {
    return xi(k) + (i * k * rho * epsilon + kappa_hat);
}

// Computes h1 = − κθ/ε² × [ d_plus * T + 2 ln((d_minus + d_plus * e^{-ξT}) / (2 ξ)) ]
std::complex<double> FourierPricer::h1(std::complex<double> k) {
    auto dp = d_plus(k);
    auto dm = d_minus(k);
    auto expt = exp(-xi(k) * T);
    return -kappa * theta / (epsilon * epsilon) *
        (dp * T + 2.0 * log((dm + dp * expt) / (2.0 * xi(k))));
}

// Computes h2 = (1 − e^{-ξT}) / (d_minus + d_plus * e^{-ξT})
std::complex<double> FourierPricer::h2(std::complex<double> k) {
    auto dp = d_plus(k);
    auto dm = d_minus(k);
    auto expt = exp(-xi(k) * T);
    return (1.0 - expt) / (dm + dp * expt);
}

// Computes the Fourier integrand:
// Integrand = exp[ (1/2 - i k) ln(X0/K) + h1(k) − (k² + 1/4) h2(k) V0 ] / (k² + 1/4)
std::complex<double> FourierPricer::integrand(double k) {
    std::complex<double> kc = k;
    auto num = exp((0.5 - i * kc) * log(X0 / K)
        + h1(kc)
        - (kc * kc + 0.25) * h2(kc) * V0);
    return num / (kc * kc + 0.25);
}

// Performs numerical integration (trapezoidal rule) of the Fourier integral
// over [-integrationBound, integrationBound] with N points
double FourierPricer::computePrice(int N, double integrationBound) {
    double a = -integrationBound;
    double b = integrationBound;
    double h = (b - a) / static_cast<double>(N);
    double integral = 0.0;

    // Sum the interior terms of the trapezoidal rule
    for (int j = 1; j < N; ++j) {
        double k = a + j * h;
        integral += integrand(k).real();  // Use only the real part
    }

    // Add the endpoints with a coefficient of 1/2 (trapezoidal rule)
    integral += 0.5 * (integrand(a).real() + integrand(b).real());
    integral *= h;

    // Final formula: X0 minus (K / 2π) times the integral
    return X0 - K / (2.0 * M_PI) * integral;
}
