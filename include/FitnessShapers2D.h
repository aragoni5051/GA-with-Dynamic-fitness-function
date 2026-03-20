#pragma once
#include "IFitness.h"
#include "DecodeUtils.h"
#include "Benchmarks2D.h"

#include <array>
#include <vector>

enum class Shaper2D {
    Fixed,
    Flattened,
    Taylor,   // (local quadratic surrogate around best-of-generation)
    Fourier   // (global Fourier surrogate over domain, truncation order grows)
};

struct ShaperConfig2D {
    Shaper2D type = Shaper2D::Fixed;

    int max_generations = 1000;

    // Flattened: alpha schedule uses max_generations
    // Taylor:
    int taylor_refit_every = 10;
    int taylor_samples = 300;
    int taylor_radius = 12;
    double taylor_surrogate_weight = 0.3; // how much surrogate influences early training

    // Fourier:
    int fourier_n0 = 1;
    int fourier_nmax = 10;
    double fourier_surrogate_weight = 0.3; // how much surrogate influences early training
};

class FitnessShaper2D : public IFitness {
public:
    FitnessShaper2D(Benchmark2D bench, BenchmarkConfig2D bench_cfg, ShaperConfig2D shaper_cfg);

    double eval(const Chromosome& c, int g, std::mt19937& rng) override;
    void on_generation_end(const Chromosome& bestC, int g, std::mt19937& rng) override;

private:
    Benchmark2D bench;
    BenchmarkConfig2D bcfg;
    ShaperConfig2D scfg;

    // Flattened
    bool has_fbar = false;
    double fbar = 0.0;

    // Taylor (quadratic model around (x0,y0))
    int tx0 = 0, ty0 = 0;
    bool has_taylor = false;
    std::array<double, 6> ta{}; // a0 + a1 u + a2 v + a3 u^2 + a4 uv + a5 v^2

    // Fourier (global coefficients computed once)
    bool has_fourier = false;
    double fc0 = 0.0;
    std::vector<double> fax, fbx, fay, fby; // 1..nmax

    // Helpers
    static double alpha_linear(int g, int G);
    static void decode_xy(const Chromosome& c, int& x, int& y);

    double base_xy(int x, int y) const;

    // Flattened
    void ensure_fbar();

    // Taylor
    double taylor_surrogate(int x, int y) const;
    void fit_taylor(std::mt19937& rng);
    static bool solve6(double M[6][6], double b[6], double out[6]);

    // Fourier
    int order_at(int g) const;
    void ensure_fourier();
    double fourier_surrogate(int x, int y, int n) const;
};