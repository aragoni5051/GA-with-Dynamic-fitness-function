#pragma once
#include <vector>

struct Chromosome {
    std::vector<bool> bits;

    Chromosome() = default;
    explicit Chromosome(int gene_length) : bits(gene_length, false) {}
};