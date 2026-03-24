#pragma once
#include "IFitness.h"
#include "DecodeUtils.h"
#include "Benchmarks2D.h"

#include <array>
#include <vector>
#include <random>

enum class Shaper2D {
    Fixed,
    Flattened,
    Taylor,
    Fourier
};

struct ShaperConfig2D {
    Shaper2D type = Shaper2D::Fixed;
    int max_generations = 1000;

    // Taylor settings
    int taylor_refit_every = 10;
    int taylor_samples = 300;
    int taylor_radius = 12;
    int taylor_n0 = 1;
    int taylor_nmax = 2;

    // Fourier settings
    int fourier_n0 = 1;
    int fourier_nmax = 127;
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

    // ---------- Flattened ----------
    bool has_fbar = false;
    double fbar = 0.0;

    // ---------- Taylor ----------
    int tx0 = 0;
    int ty0 = 0;
    bool has_taylor = false;

    // Current Taylor coefficients.
    // For now we keep quadratic storage; cpp can later generalize this if needed.
    // [a0, a1, a2, a3, a4, a5] corresponds to:
    // a0 + a1*u + a2*v + a3*u^2 + a4*u*v + a5*v^2
    std::array<double, 6> taylor_coeffs{};

    // ---------- Fourier ----------
    bool has_fourier = false;
    double fc0 = 0.0;

    // Separable Fourier coefficients:
    // Fx(x) = Σ (fax[k] cos(...) + fbx[k] sin(...))
    // Fy(y) = Σ (fay[k] cos(...) + fby[k] sin(...))
    std::vector<double> fax;
    std::vector<double> fbx;
    std::vector<double> fay;
    std::vector<double> fby;

    // ---------- Common helpers ----------
    static void decode_xy(const Chromosome& c, int& x, int& y);
    static double alpha_linear(int g, int G);
    double base_xy(int x, int y) const;

    // ---------- Flattened helpers ----------
    void ensure_fbar();

    // ---------- Taylor helpers ----------
    int order_taylor(int g) const;
    void fit_taylor(std::mt19937& rng);
    double taylor_value_xy(int x, int y, int n) const;
    static bool solve6(double M[6][6], double b[6], double out[6]);

    // ---------- Fourier helpers ----------
    int order_fourier(int g) const;
    void ensure_fourier();
    double fourier_value_xy(int x, int y, int n) const;
};