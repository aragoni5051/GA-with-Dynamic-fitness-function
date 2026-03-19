#pragma once
#include "IFitness.h"
#include "RastriginCore.h"

class FitnessFixed : public IFitness, private RastriginCore {
public:
    FitnessFixed() = default;
    double eval(const Chromosome& c, int g, std::mt19937& rng) override;
};