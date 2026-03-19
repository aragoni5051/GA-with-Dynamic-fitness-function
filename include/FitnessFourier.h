#pragma once
#include "IFitness.h"
#include "RastriginCore.h"
#include <vector>

class FitnessFourier : public IFitness, private RastriginCore {
public:
    FitnessFourier(int max_generations, int n0 = 1, int nmax = 10);

    double eval(const Chromosome& c, int g, std::mt19937& rng) override;

private:
    int G;
    int n0;
    int nmax;
    std::vector<double> amp; // amp[k] for harmonic k

    int order_at(int g) const;
    double objective_fourier(int x, int y, int n) const; // nonnegative
};