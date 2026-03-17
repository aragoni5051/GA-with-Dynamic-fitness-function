#include "Progress.h"
#include "GA.h"

#include "FitnessRastrigin.h"
#include "TournamentSelection.h"
#include "SinglePointCrossover.h"
#include "BitFlipMutation.h"

#include "DecodeUtils.h"

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
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

int main() {
    const int RUNS = 10000;
    unsigned base_seed = 12345;

    GAConfig cfg;
    cfg.population_size = 100;
    cfg.gene_length = 16;          // 8 bits for x, 8 for y (Gray)
    cfg.max_generations = 1000;
    cfg.mutation_rate = 0.02;

    int hitcnt = 0;
    const std::pair<int,int> GLOBAL_OPT = {0, 0};

    std::unordered_map<std::pair<int,int>, int, PairHash> opt_counts;
    opt_counts.reserve(RUNS);

    progress_bar(0, RUNS, "Runs");

    for (int r = 1; r <= RUNS; ++r) {
        std::mt19937 rng(base_seed + r);

        GA ga(cfg, rng,
              std::make_unique<FitnessRastrigin>(),
              std::make_unique<TournamentSelection>(3),
              std::make_unique<SinglePointCrossover>(),
              std::make_unique<BitFlipMutation>());

        ga.evaluate();

        const Chromosome& bestC = ga.best_chromosome();
        int x, y;
        decode_xy_gray_8_8(bestC.bits, x, y);

        std::pair<int,int> xy = {x, y};
        opt_counts[xy]++;

        if (xy == GLOBAL_OPT) hitcnt++;

        progress_bar(r, RUNS, "Runs");
    }

    std::vector<std::pair<std::pair<int,int>, int>> items;
    items.reserve(opt_counts.size());
    for (auto& kv : opt_counts) items.push_back(kv);

    std::sort(items.begin(), items.end(),
              [](auto& a, auto& b){ return a.second > b.second; });

    double hit_ratio = hitcnt / (double)RUNS;
    std::printf("Hit ratio (x,y)=(0,0): %.4f (%d/%d)\n\n", hit_ratio, hitcnt, RUNS);

    std::printf("Top 5 converged optima (by final best (x,y)):\n");
    int topK = (int)std::min<size_t>(5, items.size());
    for (int i = 0; i < topK; ++i) {
        auto xy = items[i].first;
        int count = items[i].second;
        double ratio = count / (double)RUNS;
        std::printf("%d) (x,y)=(%d,%d)  ratio=%.4f  (%d/%d)\n",
                    i + 1, xy.first, xy.second, ratio, count, RUNS);
    }

    return 0;
}