#ifndef FOURIER_PRICER_H
#define FOURIER_PRICER_H

#include <complex>

// This class implements the analytical formula (Fourier inversion) for pricing a Call option
// in the Heston model, as presented in Proposition 3 of Leif Andersen's paper.

class FourierPricer {
public:
    // Constructor: initializes all the parameters of the Heston model
    FourierPricer(double X0, double V0, double kappa, double theta,
        double epsilon, double rho, double T, double K);

    // Computes the exact Call option price via quadrature over the Fourier integral
    double computePrice(int N = 10000, double integrationBound = 100.0);

private:
    // Model parameters
    double X0, V0, kappa, theta, epsilon, rho, T, K;
    double kappa_hat; // κ̂ = κ - ρ * ε / 2

    // Internal functions to compute the complex quantities in the formula
    std::complex<double> h1(std::complex<double> k);
    std::complex<double> h2(std::complex<double> k);
    std::complex<double> xi(std::complex<double> k);
    std::complex<double> d_plus(std::complex<double> k);
    std::complex<double> d_minus(std::complex<double> k);

    // Integrand function to evaluate in the Fourier integral
    std::complex<double> integrand(double k);
};

#endif // FOURIER_PRICER_H
