#pragma once
#include "IFitness.h"
#include "RastriginCore.h"

class FitnessTaylor : public IFitness, private RastriginCore {
public:
    FitnessTaylor(int max_generations, int n0 = 0, int nmax = 10);

    double eval(const Chromosome& c, int g, std::mt19937& rng) override;

private:
    int G;
    int n0;
    int nmax;

    int order_at(int g) const;

    static double wrap_to_pi(double t);
    static double cos_taylor_wrapped(double t, int n);

    double objective_taylor(int x, int y, int n) const; // like Rastrigin but with cos_taylor
};