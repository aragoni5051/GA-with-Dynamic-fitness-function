#include "BitFlipMutation.h"
#include <random>
#include <stdexcept>

void BitFlipMutation::mutate(Chromosome& c, double mutation_rate, std::mt19937& rng) {
    if (mutation_rate < 0.0 || mutation_rate > 1.0)
        throw std::invalid_argument("mutation_rate must be in [0,1]");

    std::bernoulli_distribution flip(mutation_rate);

    for (size_t i = 0; i < c.bits.size(); ++i) {
        if (flip(rng)) c.bits[i] = !c.bits[i];
    }
}