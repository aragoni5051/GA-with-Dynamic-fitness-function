#pragma once
#include "IFitness.h"
#include "RastriginCore.h"

class FitnessFlattened : public IFitness, private RastriginCore {
public:
    explicit FitnessFlattened(int max_generations);

    double eval(const Chromosome& c, int g, std::mt19937& rng) override;

private:
    int G;
    double fbar; // mean of base fitness over domain (constant)

    static double alpha(int g, int G);
    double compute_fbar() const;
};