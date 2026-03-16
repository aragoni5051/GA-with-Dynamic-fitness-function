#include "DummyCrossover.h"
#include <random>

std::pair<Chromosome, Chromosome>
DummyCrossover::cross(const Chromosome& p1, const Chromosome& p2, std::mt19937& rng) {
    int n = (int)p1.bits.size();
    Chromosome c1(n), c2(n);

    std::bernoulli_distribution bit(0.5);

    for (int i = 0; i < n; ++i) {
        c1.bits[i] = bit(rng);
        c2.bits[i] = bit(rng);
    }
    return {c1, c2};
}