#ifndef QESCHEME_H
#define QESCHEME_H

#include "Variance.h"

class QEScheme : public Variance {
public:
    QEScheme(const double& psi_c);
    QEScheme(const QEScheme& other);
    QEScheme& operator=(const QEScheme& other);
    ~QEScheme() override = default;

    double step(
        const double& Vt,
        const double& dt,
        const double& kappa,
        const double& theta,
        const double& epsilon,
        std::mt19937& rng
    ) const override;

    QEScheme* clone() const override;

private:
    double psi_c;
};

#endif // QESCHEME_H

