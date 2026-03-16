#pragma once
#include "Chromosome.h"
#include <random>
#include <utility>

struct ICrossover {
    virtual ~ICrossover() = default;

    // Return two children
    virtual std::pair<Chromosome, Chromosome>
    cross(const Chromosome& p1, const Chromosome& p2, std::mt19937& rng) = 0;
};