#include "DummyMutation.h"
#include <random>

void DummyMutation::mutate(Chromosome& c, double mutation_rate, std::mt19937& rng) {
    std::bernoulli_distribution flip(mutation_rate);
    for (size_t i = 0; i < c.bits.size(); ++i) {
        if (flip(rng)) c.bits[i] = !c.bits[i];
    }
}