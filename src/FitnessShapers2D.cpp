#include "FitnessShapers2D.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

FitnessShaper2D::FitnessShaper2D(Benchmark2D bench_,
                                 BenchmarkConfig2D bench_cfg_,
                                 ShaperConfig2D shaper_cfg_)
    : bench(bench_), bcfg(bench_cfg_), scfg(shaper_cfg_) {

    if (scfg.max_generations <= 0) {
        throw std::invalid_argument("FitnessShaper2D: max_generations must be > 0");
    }

    if (scfg.type == Shaper2D::Flattened) {
        ensure_fbar();
    }

    if (scfg.type == Shaper2D::Fourier) {
        ensure_fourier();
    }
}

void FitnessShaper2D::decode_xy(const Chromosome& c, int& x, int& y) {
    decode_xy_gray_8_8(c.bits, x, y);
}

double FitnessShaper2D::alpha_linear(int g, int G) {
    if (G <= 1) return 1.0;
    g = std::clamp(g, 0, G - 1);
    return (double)g / (double)(G - 1);
}

double FitnessShaper2D::base_xy(int x, int y) const {
    return eval_benchmark(bench, x, y, bcfg);
}

// ---------- Flattened ----------
void FitnessShaper2D::ensure_fbar() {
    if (has_fbar) return;

    long double sum = 0.0;
    long long cnt = 0;

    for (int x = -128; x <= 127; ++x) {
        for (int y = -128; y <= 127; ++y) {
            sum += (long double)base_xy(x, y);
            ++cnt;
        }
    }

    fbar = (double)(sum / (long double)cnt);
    has_fbar = true;
}

// ---------- Taylor ----------
int FitnessShaper2D::order_taylor(int g) const {
    const int G = scfg.max_generations;
    const int n0 = scfg.taylor_n0;
    const int nmax = scfg.taylor_nmax;

    if (G <= 1) return nmax;

    g = std::clamp(g, 0, G - 1);
    double t = (double)g / (double)(G - 1);
    int n = (int)std::lround(n0 + t * (nmax - n0));

    // Current storage only supports up to quadratic
    return std::clamp(n, 0, 2);
}

void FitnessShaper2D::fit_taylor(std::mt19937& rng) {
    double ATA[6][6] = {};
    double ATb[6] = {};

    std::uniform_int_distribution<int> du(-scfg.taylor_radius, scfg.taylor_radius);
    std::uniform_int_distribution<int> dv(-scfg.taylor_radius, scfg.taylor_radius);

    auto clamp_xy = [](int& x, int& y) {
        x = std::clamp(x, -128, 127);
        y = std::clamp(y, -128, 127);
    };

    for (int i = 0; i < scfg.taylor_samples; ++i) {
        int x = tx0 + du(rng);
        int y = ty0 + dv(rng);
        clamp_xy(x, y);

        double u = (double)(x - tx0);
        double v = (double)(y - ty0);

        // quadratic basis
        double phi[6] = {
            1.0,
            u,
            v,
            u * u,
            u * v,
            v * v
        };

        double fx = base_xy(x, y);

        for (int r = 0; r < 6; ++r) {
            ATb[r] += phi[r] * fx;
            for (int c = 0; c < 6; ++c) {
                ATA[r][c] += phi[r] * phi[c];
            }
        }
    }

    double out[6];
    if (solve6(ATA, ATb, out)) {
        for (int i = 0; i < 6; ++i) {
            taylor_coeffs[i] = out[i];
        }
        has_taylor = true;
    } else {
        has_taylor = false;
    }
}

double FitnessShaper2D::taylor_value_xy(int x, int y, int n) const {
    double u = (double)(x - tx0);
    double v = (double)(y - ty0);

    // n = 0: constant
    if (n <= 0) {
        return taylor_coeffs[0];
    }

    // n = 1: constant + linear
    if (n == 1) {
        return taylor_coeffs[0]
             + taylor_coeffs[1] * u
             + taylor_coeffs[2] * v;
    }

    // n >= 2: full quadratic currently supported by storage
    return taylor_coeffs[0]
         + taylor_coeffs[1] * u
         + taylor_coeffs[2] * v
         + taylor_coeffs[3] * u * u
         + taylor_coeffs[4] * u * v
         + taylor_coeffs[5] * v * v;
}

bool FitnessShaper2D::solve6(double M[6][6], double b[6], double out[6]) {
    const int N = 6;

    for (int i = 0; i < N; ++i) {
        int piv = i;
        double best = std::fabs(M[i][i]);

        for (int r = i + 1; r < N; ++r) {
            double v = std::fabs(M[r][i]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }

        if (best < 1e-12) return false;

        if (piv != i) {
            for (int c = 0; c < N; ++c) {
                std::swap(M[i][c], M[piv][c]);
            }
            std::swap(b[i], b[piv]);
        }

        double diag = M[i][i];
        for (int c = i; c < N; ++c) M[i][c] /= diag;
        b[i] /= diag;

        for (int r = 0; r < N; ++r) {
            if (r == i) continue;
            double factor = M[r][i];
            if (std::fabs(factor) < 1e-18) continue;

            for (int c = i; c < N; ++c) {
                M[r][c] -= factor * M[i][c];
            }
            b[r] -= factor * b[i];
        }
    }

    for (int i = 0; i < N; ++i) {
        out[i] = b[i];
    }
    return true;
}

// ---------- Fourier ----------
int FitnessShaper2D::order_fourier(int g) const {
    const int G = scfg.max_generations;
    const int n0 = scfg.fourier_n0;
    const int nmax = scfg.fourier_nmax;

    if (G <= 1) return nmax;

    g = std::clamp(g, 0, G - 1);
    double t = (double)g / (double)(G - 1);
    int n = (int)std::lround(n0 + t * (nmax - n0));
    return std::clamp(n, n0, nmax);
}

void FitnessShaper2D::ensure_fourier() {
    if (has_fourier) return;

    const int nmax = scfg.fourier_nmax;
    fax.assign(nmax + 1, 0.0);
    fbx.assign(nmax + 1, 0.0);
    fay.assign(nmax + 1, 0.0);
    fby.assign(nmax + 1, 0.0);

    constexpr double L = 256.0;
    constexpr double two_pi = 6.2831853071795864769;

    // Separable marginals
    std::vector<double> fx_mean(256, 0.0);
    std::vector<double> fy_mean(256, 0.0);

    for (int xi = 0; xi < 256; ++xi) {
        int x = xi - 128;
        long double sum = 0.0;
        for (int yi = 0; yi < 256; ++yi) {
            int y = yi - 128;
            sum += (long double)base_xy(x, y);
        }
        fx_mean[xi] = (double)(sum / 256.0L);
    }

    for (int yi = 0; yi < 256; ++yi) {
        int y = yi - 128;
        long double sum = 0.0;
        for (int xi = 0; xi < 256; ++xi) {
            int x = xi - 128;
            sum += (long double)base_xy(x, y);
        }
        fy_mean[yi] = (double)(sum / 256.0L);
    }

    long double overall = 0.0;
    for (int xi = 0; xi < 256; ++xi) {
        overall += (long double)fx_mean[xi];
    }
    fc0 = (double)(overall / 256.0L);

    auto idx_to_coord = [](int i) { return i - 128; };

    for (int k = 1; k <= nmax; ++k) {
        long double ax = 0.0, bx = 0.0, ay = 0.0, by = 0.0;

        for (int xi = 0; xi < 256; ++xi) {
            double x = (double)idx_to_coord(xi);
            double w = two_pi * k * x / L;
            ax += (long double)fx_mean[xi] * std::cos(w);
            bx += (long double)fx_mean[xi] * std::sin(w);
        }

        for (int yi = 0; yi < 256; ++yi) {
            double y = (double)idx_to_coord(yi);
            double w = two_pi * k * y / L;
            ay += (long double)fy_mean[yi] * std::cos(w);
            by += (long double)fy_mean[yi] * std::sin(w);
        }

        double scale = 2.0 / 256.0;
        fax[k] = (double)(ax * scale);
        fbx[k] = (double)(bx * scale);
        fay[k] = (double)(ay * scale);
        fby[k] = (double)(by * scale);
    }

    has_fourier = true;
}

double FitnessShaper2D::fourier_value_xy(int x, int y, int n) const {
    constexpr double L = 256.0;
    constexpr double two_pi = 6.2831853071795864769;

    const double wx_base = two_pi * (double)x / L;
    const double wy_base = two_pi * (double)y / L;

    double val = fc0;

    for (int k = 1; k <= n; ++k) {
        double wx = k * wx_base;
        double wy = k * wy_base;

        val += fax[k] * std::cos(wx) + fbx[k] * std::sin(wx);
        val += fay[k] * std::cos(wy) + fby[k] * std::sin(wy);
    }

    return val;
}

// ---------- Main eval ----------
double FitnessShaper2D::eval(const Chromosome& c, int g, std::mt19937& rng) {
    (void)rng;

    int x, y;
    decode_xy(c, x, y);

    switch (scfg.type) {
    case Shaper2D::Fixed:
        return base_xy(x, y);

    case Shaper2D::Flattened: {
        ensure_fbar();
        double fb = base_xy(x, y);
        double a = alpha_linear(g, scfg.max_generations);
        return (1.0 - a) * fbar + a * fb;
    }

    case Shaper2D::Taylor: {
        if (!has_taylor) {
            // until fitted, fall back to base
            return base_xy(x, y);
        }
        int n = order_taylor(g);
        return taylor_value_xy(x, y, n);
    }

    case Shaper2D::Fourier: {
        ensure_fourier();
        int n = order_fourier(g);
        return fourier_value_xy(x, y, n);
    }

    default:
        return base_xy(x, y);
    }
}

void FitnessShaper2D::on_generation_end(const Chromosome& bestC, int g, std::mt19937& rng) {
    if (scfg.type != Shaper2D::Taylor) return;
    if (g % scfg.taylor_refit_every != 0) return;

    decode_xy(bestC, tx0, ty0);
    fit_taylor(rng);
}