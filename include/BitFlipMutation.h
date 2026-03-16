#pragma once
#include "IMutation.h"

class BitFlipMutation : public IMutation {
public:
    void mutate(Chromosome& c, double mutation_rate, std::mt19937& rng) override;
};