#include "Progress.h"
#include "GA.h"

#include "DummyFitness.h"
#include "DummySelection.h"
#include "DummyCrossover.h"
#include "DummyMutation.h"

#include <random>
#include <cstdio>

int main() {
    const int EPOCHS = 1000;
    unsigned base_seed = 12345;

    GAConfig cfg;
    cfg.population_size = 100;
    cfg.gene_length = 32;
    cfg.max_generations = 100;
    cfg.mutation_rate = 0.02;

    double best_overall = 0.0;

    progress_bar(0, EPOCHS, "Epochs");

    for (int e = 1; e <= EPOCHS; ++e) {
        std::mt19937 rng(base_seed + e);

        GA ga(cfg, rng,
              std::make_unique<DummyFitness>(),
              std::make_unique<DummySelection>(),
              std::make_unique<DummyCrossover>(),
              std::make_unique<DummyMutation>());

        ga.evaluate();
        if (ga.best_fitness() > best_overall) best_overall = ga.best_fitness();

        progress_bar(e, EPOCHS, "Epochs");
    }

    std::printf("Best fitness seen (dummy): %.6f\n", best_overall);
    return 0;
}