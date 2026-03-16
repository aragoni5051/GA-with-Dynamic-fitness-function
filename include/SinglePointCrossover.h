#pragma once
#include "ICrossover.h"

class SinglePointCrossover : public ICrossover {
public:
    std::pair<Chromosome, Chromosome>
    cross(const Chromosome& p1, const Chromosome& p2, std::mt19937& rng) override;
};