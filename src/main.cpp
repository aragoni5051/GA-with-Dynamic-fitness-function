#include "Progress.h"
#include "GA.h"

#include "FitnessRastrigin.h"
#include "TournamentSelection.h"
#include "SinglePointCrossover.h"
#include "BitFlipMutation.h"

#include "DecodeUtils.h"   // <-- decode best chromosome to (x,y)

#include <random>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility>

struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        size_t h1 = std::hash<int>{}(p.first);
        size_t h2 = std::hash<int>{}(p.second);
        // standard-ish hash combine
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

int main() {
    const int EPOCHS = 1000;
    unsigned base_seed = 12345;

    GAConfig cfg;
    cfg.population_size = 100;
    cfg.gene_length = 16;          // 8 bits for x, 8 for y (Gray)
    cfg.max_generations = 100;
    cfg.mutation_rate = 0.02;

    int hitcnt = 0;
    const std::pair<int,int> GLOBAL_OPT = {0, 0};

    std::unordered_map<std::pair<int,int>, int, PairHash> opt_counts;
    opt_counts.reserve(EPOCHS);

    progress_bar(0, EPOCHS, "Epochs");

    for (int e = 1; e <= EPOCHS; ++e) {
        std::mt19937 rng(base_seed + e);

        GA ga(cfg, rng,
              std::make_unique<FitnessRastrigin>(),
              std::make_unique<TournamentSelection>(3),
              std::make_unique<SinglePointCrossover>(),
              std::make_unique<BitFlipMutation>());

        ga.evaluate();

        // Decode the best chromosome at the end of the epoch
        const Chromosome& bestC = ga.best_chromosome();
        int x, y;
        decode_xy_gray_8_8(bestC.bits, x, y);

        std::pair<int,int> xy = {x, y};
        opt_counts[xy]++;

        if (xy == GLOBAL_OPT) hitcnt++;

        progress_bar(e, EPOCHS, "Epochs");
    }

    // Convert map to vector and sort by count descending
    std::vector<std::pair<std::pair<int,int>, int>> items;
    items.reserve(opt_counts.size());
    for (auto& kv : opt_counts) items.push_back(kv);

    std::sort(items.begin(), items.end(),
              [](auto& a, auto& b){ return a.second > b.second; });

    double hit_ratio = hitcnt / (double)EPOCHS;
    std::printf("Hit ratio (x,y)=(0,0): %.4f (%d/%d)\n\n", hit_ratio, hitcnt, EPOCHS);

    std::printf("Top 5 converged optima (by final best (x,y)):\n");
    int topK = (int)std::min<size_t>(5, items.size());
    for (int i = 0; i < topK; ++i) {
        auto xy = items[i].first;
        int count = items[i].second;
        double ratio = count / (double)EPOCHS;
        std::printf("%d) (x,y)=(%d,%d)  ratio=%.4f  (%d/%d)\n",
                    i + 1, xy.first, xy.second, ratio, count, EPOCHS);
    }

    return 0;
}