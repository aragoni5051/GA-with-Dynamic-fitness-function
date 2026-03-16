#include "GA.h"
#include <algorithm>
#include <random>

GA::GA(const GAConfig& cfg_,
       std::mt19937& rng_,
       std::unique_ptr<IFitness> fit_,
       std::unique_ptr<ISelection> sel_,
       std::unique_ptr<ICrossover> cross_,
       std::unique_ptr<IMutation> mut_)
    : cfg(cfg_),
      rng(rng_),
      fit(std::move(fit_)),
      sel(std::move(sel_)),
      cross(std::move(cross_)),
      mut(std::move(mut_)) {

    population.resize(cfg.population_size);
    fitness_vec.resize(cfg.population_size, 0.0);
}

void GA::init_population() {
    std::bernoulli_distribution bit(0.5);
    for (int i = 0; i < cfg.population_size; ++i) {
        population[i] = Chromosome(cfg.gene_length);
        for (int j = 0; j < cfg.gene_length; ++j) {
            population[i].bits[j] = bit(rng);
        }
    }
}

void GA::evaluate_population(int gen) {
    for (int i = 0; i < cfg.population_size; ++i) {
        fitness_vec[i] = fit->eval(population[i], gen, rng);
    }
}

void GA::evaluate() {
    init_population();

    for (int gen = 0; gen < cfg.max_generations; ++gen) {
        evaluate_population(gen);

        std::vector<Chromosome> next;
        next.reserve(cfg.population_size);

        while ((int)next.size() < cfg.population_size) {
            int i = sel->pick_parent(population, fitness_vec, rng);
            int j = sel->pick_parent(population, fitness_vec, rng);
            if (i < 0 || j < 0) break;

            auto [c1, c2] = cross->cross(population[i], population[j], rng);

            mut->mutate(c1, cfg.mutation_rate, rng);
            mut->mutate(c2, cfg.mutation_rate, rng);

            next.push_back(std::move(c1));
            if ((int)next.size() < cfg.population_size) next.push_back(std::move(c2));
        }

        population = std::move(next);
    }

    // final evaluation so best_fitness() is meaningful after run
    evaluate_population(cfg.max_generations);
}

double GA::best_fitness() const {
    return *std::max_element(fitness_vec.begin(), fitness_vec.end());
}