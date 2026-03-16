#include "DummySelection.h"
#include <random>

int DummySelection::pick_parent(const std::vector<Chromosome>& population,
                                const std::vector<double>&,
                                std::mt19937& rng) {
    if (population.empty()) return -1;
    std::uniform_int_distribution<int> dist(0, (int)population.size() - 1);
    return dist(rng);
}