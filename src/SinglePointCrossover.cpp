#include "SinglePointCrossover.h"
#include <random>
#include <stdexcept>

std::pair<Chromosome, Chromosome>
SinglePointCrossover::cross(const Chromosome& p1, const Chromosome& p2, std::mt19937& rng) {
    int n = (int)p1.bits.size();
    if (n != (int)p2.bits.size()) throw std::invalid_argument("crossover: length mismatch");
    if (n != 16) throw std::invalid_argument("crossover: expected 16 bits (8 for x, 8 for y)");

    Chromosome c1(n), c2(n);

    std::uniform_int_distribution<int> cut_x(1, 7);   // cut inside [0..7]
    std::uniform_int_distribution<int> cut_y(9, 15);  // cut inside [8..15]

    int px = cut_x(rng);
    int py = cut_y(rng);

    // X half
    for (int i = 0; i < 8; ++i) {
        if (i < px) { c1.bits[i] = p1.bits[i]; c2.bits[i] = p2.bits[i]; }
        else        { c1.bits[i] = p2.bits[i]; c2.bits[i] = p1.bits[i]; }
    }

    // Y half
    for (int i = 8; i < 16; ++i) {
        if (i < py) { c1.bits[i] = p1.bits[i]; c2.bits[i] = p2.bits[i]; }
        else        { c1.bits[i] = p2.bits[i]; c2.bits[i] = p1.bits[i]; }
    }

    return {c1, c2};
}