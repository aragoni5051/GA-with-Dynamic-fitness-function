#include "GA.h"
#include <algorithm>
#include <random>
#include <stdexcept>

// Termination condition: population is "converged" if every bit position has
// >= dom dominance of either 0 or 1 (e.g., dom = 0.95).
static bool converged_dominance(const std::vector<Chromosome>& pop,
                                int gene_length,
                                double dom = 0.95) {
    int N = (int)pop.size();
    if (N == 0) return false;

    for (int j = 0; j < gene_length; ++j) {
        int ones = 0;
        for (int i = 0; i < N; ++i) {
            ones += pop[i].bits[j] ? 1 : 0;
        }
        double p = (double)ones / (double)N;
        double dominance = (p > 1.0 - p) ? p : (1.0 - p);
        if (dominance < dom) return false; // not converged at this bit
    }
    return true;
}

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

    population.assign(cfg.population_size, Chromosome(cfg.gene_length));
    fitness_vec.assign(cfg.population_size, 0.0);
}

void GA::init_population_gray() {
    // We randomize bits. DecodeUtils interprets them as Gray when decoding.
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
    init_population_gray();

    int last_gen = 0;

    for (int gen = 0; gen < cfg.max_generations; ++gen) {
        last_gen = gen;

        evaluate_population(gen);

        // find best under current training fitness (fitness_vec)
        int best_i = 0;
        for (int i = 1; i < (int)fitness_vec.size(); ++i)
            if (fitness_vec[i] > fitness_vec[best_i]) best_i = i;

        // let fitness update (for local Taylor / Fourier)
        fit->on_generation_end(population[best_i], gen, rng);

        // Terminate if population converged (>=95% dominance on every gene/bit)
        if (converged_dominance(population, cfg.gene_length, 0.95)) {
            break;
        }

        std::vector<Chromosome> next;
        next.reserve(cfg.population_size);

        while ((int)next.size() < cfg.population_size) {
            int i = sel->pick_parent(population, fitness_vec, rng);
            int j = sel->pick_parent(population, fitness_vec, rng);

            if (i < 0 || j < 0) {
                throw std::runtime_error("Selection failed: returned invalid index");
            }

            auto [c1, c2] = cross->cross(population[i], population[j], rng);

            mut->mutate(c1, cfg.mutation_rate, rng);
            mut->mutate(c2, cfg.mutation_rate, rng);

            next.push_back(std::move(c1));
            if ((int)next.size() < cfg.population_size) next.push_back(std::move(c2));
        }

        population = std::move(next);
    }

    // Final eval using the LAST generation we actually reached
    evaluate_population(last_gen);
}

double GA::best_fitness() const {
    return *std::max_element(fitness_vec.begin(), fitness_vec.end());
}

const Chromosome& GA::best_chromosome() const {
    int best_i = 0;
    for (int i = 1; i < (int)fitness_vec.size(); ++i) {
        if (fitness_vec[i] > fitness_vec[best_i]) best_i = i;
    }
    return population[best_i];
}