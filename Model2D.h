#ifndef MODEL2D_H
#define MODEL2D_H

#include <random>
#include <utility>

// Interface pour tout modèle 2D (Heston, SABR, etc)
class Model2D {
public:
    virtual ~Model2D() = default;

    // Simule le prochain (S, V) connaissant (S, V, dt)
    virtual std::pair<double, double> step(
        double St,
        double Vt,
        double dt,
        std::mt19937& rng
    ) const = 0;

    virtual Model2D* clone() const = 0;
};

#endif // MODEL2D_H

