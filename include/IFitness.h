#pragma once
#include "Chromosome.h"
#include <random>

struct IFitness {
    virtual ~IFitness() = default;

    // Return training fitness used by GA at generation g
    virtual double eval(const Chromosome& c, int g, std::mt19937& rng) = 0;

    // OPTIONAL: called once per generation so dynamic fitness can update its model
    // Default: do nothing
    virtual void on_generation_end(const Chromosome&, int, std::mt19937&) {}
};