#include "FitnessFlattened.h"
#include <stdexcept>
#include <algorithm>

FitnessFlattened::FitnessFlattened(int max_generations)
    : G(max_generations), fbar(0.0) {
    if (G <= 0) throw std::invalid_argument("FitnessFlattened: max_generations must be > 0");
    fbar = compute_fbar();
}

double FitnessFlattened::alpha(int g, int G) {
    if (G <= 1) return 1.0;
    g = std::clamp(g, 0, G - 1);
    return (double)g / (double)(G - 1); // linear 0->1
}

double FitnessFlattened::compute_fbar() const {
    long long count = 0;
    long double sum = 0.0;

    for (int x = -128; x <= 127; ++x) {
        for (int y = -128; y <= 127; ++y) {
            sum += (long double)base_fitness(x, y);
            ++count;
        }
    }
    return (double)(sum / (long double)count);
}

double FitnessFlattened::eval(const Chromosome& c, int g, std::mt19937& rng) {
    (void)rng;
    double a = alpha(g, G);
    double fb = base_fitness(c);
    return (1.0 - a) * fbar + a * fb;
}