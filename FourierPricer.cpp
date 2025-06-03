#define _USE_MATH_DEFINES
#include "FourierPricer.h"
#include <cmath>
#include <complex>
#include <iostream>

// Définition de l'unité imaginaire i = sqrt(-1)
const std::complex<double> i(0.0, 1.0);

// Constructeur : on stocke les paramètres du modèle de Heston
FourierPricer::FourierPricer(double X0, double V0, double kappa, double theta,
    double epsilon, double rho, double T, double K)
    : X0(X0), V0(V0), kappa(kappa), theta(theta),
    epsilon(epsilon), rho(rho), T(T), K(K)
{
    // Calcul de kappa chapeau (κ̂) = κ - ρ * ε / 2
    kappa_hat = kappa - rho * epsilon / 2.0;
}

// ξ = sqrt(k² ε² (1 - ρ²) + 2 i k ε ρ κ̂ + κ̂² + ε²/4)
std::complex<double> FourierPricer::xi(std::complex<double> k) {
    return sqrt(pow(k * epsilon, 2.0) * (1.0 - rho * rho)
        + 2.0 * i * k * epsilon * rho * kappa_hat
        + pow(kappa_hat, 2.0)
        + pow(epsilon, 2.0) / 4.0);
}

// ∂₊ = ξ − (i k ε ρ + κ̂)
std::complex<double> FourierPricer::d_plus(std::complex<double> k) {
    return xi(k) - (i * k * rho * epsilon + kappa_hat);
}

// ∂₋ = ξ + (i k ε ρ + κ̂)
std::complex<double> FourierPricer::d_minus(std::complex<double> k) {
    return xi(k) + (i * k * rho * epsilon + kappa_hat);
}

// h₁ = − κθ/ε² × [ ∂₊ T + 2 ln((∂₋ + ∂₊ e^{-ξT}) / (2 ξ)) ]
std::complex<double> FourierPricer::h1(std::complex<double> k) {
    auto dp = d_plus(k);
    auto dm = d_minus(k);
    auto expt = exp(-xi(k) * T);
    return -kappa * theta / (epsilon * epsilon) *
        (dp * T + 2.0 * log((dm + dp * expt) / (2.0 * xi(k))));
}

// h₂ = (1 − e^{-ξT}) / (∂₋ + ∂₊ e^{-ξT})
std::complex<double> FourierPricer::h2(std::complex<double> k) {
    auto dp = d_plus(k);
    auto dm = d_minus(k);
    auto expt = exp(-xi(k) * T);
    return (1.0 - expt) / (dm + dp * expt);
}

// Integrand = exp[ (1/2 - i k) ln(X0/K) + h₁(k) − (k² + 1/4) h₂(k) V₀ ] / (k² + 1/4)
std::complex<double> FourierPricer::integrand(double k) {
    std::complex<double> kc = k;
    auto num = exp((0.5 - i * kc) * log(X0 / K)
        + h1(kc)
        - (kc * kc + 0.25) * h2(kc) * V0);
    return num / (kc * kc + 0.25);
}

// Intégration numérique (trapèzes) de l'intégrale de Fourier
// sur [-integrationBound, integrationBound] avec N points
double FourierPricer::computePrice(int N, double integrationBound) {
    double a = -integrationBound;
    double b = integrationBound;
    double h = (b - a) / static_cast<double>(N);
    double integral = 0.0;

    // Somme des termes intérieurs
    for (int j = 1; j < N; ++j) {
        double k = a + j * h;
        integral += integrand(k).real();  // On prend la partie réelle uniquement
    }

    // Ajout des extrémités avec coefficient 1/2 (formule des trapèzes)
    integral += 0.5 * (integrand(a).real() + integrand(b).real());
    integral *= h;

    // Formule finale : X(0) - (K / 2π) × intégrale
    return X0 - K / (2.0 * M_PI) * integral;
}
