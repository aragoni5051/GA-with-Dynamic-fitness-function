#pragma once
#include "ISelection.h"

class TournamentSelection : public ISelection {
    int k; // tournament size
public:
    explicit TournamentSelection(int k = 3); // k=3 is a good default

    int pick_parent(const std::vector<Chromosome>& population,
                    const std::vector<double>& fitness,
                    std::mt19937& rng) override;
};