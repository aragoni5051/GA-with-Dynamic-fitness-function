#pragma once
#include "ISelection.h"

// Picks a random parent index uniformly
class DummySelection : public ISelection {
public:
    int pick_parent(const std::vector<Chromosome>& population,
                    const std::vector<double>& fitness,
                    std::mt19937& rng) override;
};