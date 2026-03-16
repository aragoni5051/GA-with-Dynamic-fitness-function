#pragma once
#include "ICrossover.h"

// Produces random children (same length as parents)
class DummyCrossover : public ICrossover {
public:
    std::pair<Chromosome, Chromosome>
    cross(const Chromosome& p1, const Chromosome& p2, std::mt19937& rng) override;
};