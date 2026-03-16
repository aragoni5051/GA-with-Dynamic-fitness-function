#pragma once
#include "IFitness.h"

// Returns random fitness in [0, 1)
class DummyFitness : public IFitness {
public:
    double eval(const Chromosome& c, int g, std::mt19937& rng) override;
};