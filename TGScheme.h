#ifndef TGSCHEME_H
#define TGSCHEME_H

#include "Variance.h"
#include <vector>
#include <random>
#include <mutex>

class TGScheme : public Variance {
public:
    TGScheme(double psi_min = 0.001, double psi_max = 10.0, int grid_size = 5000);
    TGScheme(const TGScheme& other);
    TGScheme& operator=(const TGScheme& other);
    ~TGScheme() override = default;

    double step(const double& Vt,
        const double& dt,
        const double& kappa,
        const double& theta,
        const double& epsilon,
        std::mt19937& rng) const override;

    TGScheme* clone() const override;

private:
    double psi_min_, psi_max_;
    int grid_size_;

    std::vector<double> psi_grid_;
    std::vector<double> r_grid_;

    void initialize_grid();
    double compute_r(double psi) const;
    static double objective_function(double r, double psi);
    static double objective_derivative(double r, double psi);
    static double phi(double x);
    static double Phi(double x);

    mutable std::mutex cache_mutex_;
};

#endif // TGSCHEME_H