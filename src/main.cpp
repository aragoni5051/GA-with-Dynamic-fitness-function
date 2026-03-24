#include "Progress.h"
#include "GA.h"

#include "TournamentSelection.h"
#include "SinglePointCrossover.h"
#include "BitFlipMutation.h"

#include "DecodeUtils.h"

#include "Benchmarks2D.h"
#include "FitnessShapers2D.h"

#include <chrono>
#include <random>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <memory>
#include <string>
#include <cstddef>

struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        std::size_t h1 = std::hash<int>{}(p.first);
        std::size_t h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

template <typename FitnessFactory>
static void run_experiment(const char* name,
                           int RUNS,
                           unsigned base_seed,
                           const GAConfig& cfg,
                           const BenchmarkConfig2D& bcfg,
                           FitnessFactory make_fitness) {
    using clock = std::chrono::steady_clock;

    int hitcnt = 0;
    std::unordered_map<std::pair<int,int>, int, PairHash> opt_counts;
    opt_counts.reserve((std::size_t)RUNS);

    progress_bar(0, RUNS, name);

    auto t0 = clock::now();

    for (int r = 1; r <= RUNS; ++r) {
        std::mt19937 rng(base_seed + r);

        GA ga(cfg, rng,
              make_fitness(),
              std::make_unique<TournamentSelection>(3),
              std::make_unique<SinglePointCrossover>(),
              std::make_unique<BitFlipMutation>());

        ga.evaluate();

        const Chromosome& bestC = ga.best_chromosome();
        int x, y;
        decode_xy_gray_8_8(bestC.bits, x, y);

        opt_counts[{x, y}]++;
        if (x == bcfg.x0 && y == bcfg.y0) hitcnt++;

        progress_bar(r, RUNS, name);
    }

    auto t1 = clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double avg_ms = total_ms / (double)RUNS;

    std::vector<std::pair<std::pair<int,int>, int>> items;
    items.reserve(opt_counts.size());
    for (auto& kv : opt_counts) items.push_back(kv);

    std::sort(items.begin(), items.end(),
              [](auto& a, auto& b){ return a.second > b.second; });

    double hit_ratio = hitcnt / (double)RUNS;

    std::printf("\n=== %s ===\n", name);
    std::printf("Avg time per run: %.3f ms\n", avg_ms);
    std::printf("Hit ratio (x,y)=(%d,%d): %.4f (%d/%d)\n\n",
                bcfg.x0, bcfg.y0, hit_ratio, hitcnt, RUNS);

    std::printf("Top 5 converged optima (by final best (x,y)):\n");
    int topK = (int)std::min<std::size_t>(5, items.size());
    for (int i = 0; i < topK; ++i) {
        auto xy = items[i].first;
        int count = items[i].second;
        double ratio = count / (double)RUNS;
        std::printf("%d) (x,y)=(%d,%d)  ratio=%.4f  (%d/%d)\n",
                    i + 1, xy.first, xy.second, ratio, count, RUNS);
    }

    std::printf("\n");
}

int main() {
    const int RUNS = 1000;
    unsigned base_seed = 12345;

    GAConfig cfg;
    cfg.population_size = 100;
    cfg.gene_length = 16;
    cfg.max_generations = 1000;
    cfg.mutation_rate = 0.02;

    Benchmark2D bench = Benchmark2D::Rastrigin;

    BenchmarkConfig2D bcfg;
    bcfg.x0 = 51;
    bcfg.y0 = 34;
    bcfg.rastrigin_scale = 25.0;

    run_experiment("Runs (Fixed)", RUNS, base_seed, cfg, bcfg, [&]() {
        ShaperConfig2D scfg;
        scfg.type = Shaper2D::Fixed;
        scfg.max_generations = cfg.max_generations;
        return std::make_unique<FitnessShaper2D>(bench, bcfg, scfg);
    });

    /*run_experiment("Runs (Flattened)", RUNS, base_seed, cfg, bcfg, [&]() {
        ShaperConfig2D scfg;
        scfg.type = Shaper2D::Flattened;
        scfg.max_generations = cfg.max_generations;
        return std::make_unique<FitnessShaper2D>(bench, bcfg, scfg);
    });*/

    // Taylor disabled for now
    /*
    run_experiment("Runs (Taylor)", RUNS, base_seed, cfg, bcfg, [&]() {
        ShaperConfig2D scfg;
        scfg.type = Shaper2D::Taylor;
        scfg.max_generations = cfg.max_generations;
        scfg.taylor_refit_every = 10;
        scfg.taylor_samples = 300;
        scfg.taylor_radius = 12;
        scfg.taylor_surrogate_weight = 0.3;
        return std::make_unique<FitnessShaper2D>(bench, bcfg, scfg);
    });
    */

    run_experiment("Runs (Fourier)", RUNS, base_seed, cfg, bcfg, [&]() {
        ShaperConfig2D scfg;
        scfg.type = Shaper2D::Fourier;
        
        scfg.max_generations = cfg.max_generations;
        scfg.fourier_n0 = 1;
        scfg.fourier_nmax = 127;
        return std::make_unique<FitnessShaper2D>(bench, bcfg, scfg);
    });

    return 0;
}