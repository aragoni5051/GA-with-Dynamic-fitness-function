#pragma once
#include "Chromosome.h"
#include <random>
#include <vector>

struct ISelection {
    virtual ~ISelection() = default;

    // Return an index into population
    virtual int pick_parent(const std::vector<Chromosome>& population,
                            const std::vector<double>& fitness,
                            std::mt19937& rng) = 0;
};