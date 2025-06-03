#ifndef FOURIER_PRICER_H
#define FOURIER_PRICER_H

#include <complex>

// Cette classe implémente la formule analytique (Fourier inverse) du prix d'une option Call
// dans le modèle de Heston, tel que présenté dans la Proposition 3 de l'article de Leif Andersen.

class FourierPricer {
public:
    // Constructeur : initialise tous les paramètres du modèle de Heston
    FourierPricer(double X0, double V0, double kappa, double theta,
        double epsilon, double rho, double T, double K);

    // Calcule le prix exact de l'option call par quadrature sur l'intégrale de Fourier
    double computePrice(int N = 10000, double integrationBound = 100.0);

private:
    // Paramètres du modèle
    double X0, V0, kappa, theta, epsilon, rho, T, K;
    double kappa_hat; // κ̂ = κ - ρε/2

    // Fonctions internes pour calculer les quantités complexes de la formule
    std::complex<double> h1(std::complex<double> k);
    std::complex<double> h2(std::complex<double> k);
    std::complex<double> xi(std::complex<double> k);
    std::complex<double> d_plus(std::complex<double> k);
    std::complex<double> d_minus(std::complex<double> k);

    // Fonction intégrande à évaluer dans l'intégrale de Fourier
    std::complex<double> integrand(double k);
};

#endif // FOURIER_PRICER_H

