#pragma once
#include "IMutation.h"

// Randomly flips bits with given probability
class DummyMutation : public IMutation {
public:
    void mutate(Chromosome& c, double mutation_rate, std::mt19937& rng) override;
};