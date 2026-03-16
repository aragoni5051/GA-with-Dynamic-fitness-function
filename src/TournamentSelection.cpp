#include "TournamentSelection.h"
#include <random>
#include <stdexcept>

TournamentSelection::TournamentSelection(int k_) : k(k_) {
    if (k < 2) throw std::invalid_argument("TournamentSelection: k must be >= 2");
}

int TournamentSelection::pick_parent(const std::vector<Chromosome>& population,
                                     const std::vector<double>& fitness,
                                     std::mt19937& rng) {
    int n = (int)population.size();
    if (n == 0) return -1;
    if ((int)fitness.size() != n) return -1;

    std::uniform_int_distribution<int> pick(0, n - 1);

    int best = pick(rng);
    double best_fit = fitness[best];

    for (int t = 1; t < k; ++t) {
        int cand = pick(rng);
        if (fitness[cand] > best_fit) {
            best = cand;
            best_fit = fitness[cand];
        }
    }
    return best;
}