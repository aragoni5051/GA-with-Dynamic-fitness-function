#pragma once
#include "Chromosome.h"
#include <random>

struct IMutation {
    virtual ~IMutation() = default;

    virtual void mutate(Chromosome& c, double mutation_rate, std::mt19937& rng) = 0;
};