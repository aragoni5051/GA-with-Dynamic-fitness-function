#pragma once
#include <cmath>

enum class Benchmark2D {
    Rastrigin,
    Sphere,
    Rosenbrock,
    Schwefel,
    Alpine
};

struct BenchmarkConfig2D {
    // Define the known optimum location for shifted benchmarks
    int x0 = 12;
    int y0 = 34;

    // Rastrigin scaling: 128/25 ≈ 5.12
    double rastrigin_scale = 25.0;
};

// Return fitness to MAXIMIZE.
// Convention: for all benchmarks here, the global optimum is at (cfg.x0, cfg.y0).
inline double eval_benchmark(Benchmark2D type, int x, int y, const BenchmarkConfig2D& cfg) {
    // Shift once, use everywhere
    const double sx = (double)(x - cfg.x0);
    const double sy = (double)(y - cfg.y0);

    switch (type) {
    case Benchmark2D::Sphere: {
        // Maximize negative sphere: optimum at (0,0) in shifted coords
        return -(sx * sx + sy * sy);
    }

    case Benchmark2D::Rosenbrock: {
        // Standard Rosenbrock has minimum at (1,1).
        // Map shifted coords so minimum becomes (0,0) -> (cfg.x0,cfg.y0):
        // X = sx + 1, Y = sy + 1
        constexpr double a = 1.0;
        constexpr double b = 100.0;

        const double X = sx + 1.0;
        const double Y = sy + 1.0;

        const double term1 = (a - X) * (a - X);
        const double term2 = b * (Y - X * X) * (Y - X * X);

        return -(term1 + term2);
    }

    case Benchmark2D::Schwefel: {
        // Schwefel 2.26 has global minimum at 420.9687 per dimension.
        // Shift so the minimum occurs at (sx,sy)=(0,0) -> (cfg.x0,cfg.y0).
        constexpr double C = 418.9829;
        constexpr double OPT = 420.9687;

        const double u = sx + OPT;
        const double v = sy + OPT;

        auto term = [](double t) {
            return t * std::sin(std::sqrt(std::fabs(t)));
        };

        const double S = 2.0 * C - (term(u) + term(v));
        return -S; // maximize
    }

    case Benchmark2D::Alpine: {
        // Alpine N.1: minimize |t sin(t) + 0.1 t| (separable)
        auto alpine1 = [](double t) {
            return std::fabs(t * std::sin(t) + 0.1 * t);
        };

        const double A = alpine1(sx) + alpine1(sy);
        return -A; // maximize
    }

    case Benchmark2D::Rastrigin:
    default: {
        // Rastrigin: optimum at (0,0) in shifted coords
        constexpr double A = 10.0;
        constexpr double two_pi = 6.2831853071795864769;

        const double xs = sx / cfg.rastrigin_scale;
        const double ys = sy / cfg.rastrigin_scale;

        const double R = 2.0 * A
                       + (xs * xs + ys * ys)
                       - A * std::cos(two_pi * xs)
                       - A * std::cos(two_pi * ys);

        return -R; // maximize
    }
    }
}