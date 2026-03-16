#pragma once
#include "IFitness.h"

class FitnessRastrigin : public IFitness {
public:
    double eval(const Chromosome& c, int gen, std::mt19937& rng) override;
};