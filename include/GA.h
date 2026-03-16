#pragma once
#include <memory>
#include <random>
#include <vector>

#include "Chromosome.h"
#include "IFitness.h"
#include "ISelection.h"
#include "ICrossover.h"
#include "IMutation.h"

struct GAConfig {
    int population_size = 50;
    int gene_length = 32;
    int max_generations = 200;
    double mutation_rate = 0.01; // probability per bit
};

class GA {
public:
    GA(const GAConfig& cfg,
       std::mt19937& rng,
       std::unique_ptr<IFitness> fit,
       std::unique_ptr<ISelection> sel,
       std::unique_ptr<ICrossover> cross,
       std::unique_ptr<IMutation> mut);

    void evaluate();            // runs the full GA
    double best_fitness() const;

private:
    GAConfig cfg;
    std::mt19937& rng;

    std::vector<Chromosome> population;
    std::vector<double> fitness_vec;

    std::unique_ptr<IFitness> fit;
    std::unique_ptr<ISelection> sel;
    std::unique_ptr<ICrossover> cross;
    std::unique_ptr<IMutation> mut;

    void init_population();
    void evaluate_population(int gen);
};