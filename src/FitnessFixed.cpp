#include "FitnessFixed.h"

double FitnessFixed::eval(const Chromosome& c, int, std::mt19937&) {
    return base_fitness(c);
}