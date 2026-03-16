#include "DummyFitness.h"
#include <random>

double DummyFitness::eval(const Chromosome&, int, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}