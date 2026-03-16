#pragma once
#include "Chromosome.h"
#include <random>

struct IFitness {
    virtual ~IFitness() = default;

    // Return fitness value for a chromosome at generation g
    virtual double eval(const Chromosome& c, int g, std::mt19937& rng) = 0;
};