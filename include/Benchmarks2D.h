#pragma once
#include <cmath>

enum class Benchmark2D {
    Rastrigin,
    Sphere
    // Add more later: Ackley, Rosenbrock, Griewank, ...
};

struct BenchmarkConfig2D {
    // “Target” / shift (optional). For unknown optima set these to 0.
    int x0 = 21;
    int y0 = 34;

    // Rastrigin scale: 128/25 ≈ 5.12
    double rastrigin_scale = 25.0;
};

// Return fitness to MAXIMIZE
inline double eval_benchmark(Benchmark2D type, int x, int y, const BenchmarkConfig2D& cfg) {
    int dx = x - cfg.x0;
    int dy = y - cfg.y0;

    switch (type) {
    case Benchmark2D::Sphere: {
        // Maximize negative sphere
        return -(double)(dx * dx + dy * dy);
    }
    case Benchmark2D::Rastrigin:
    default: {
        constexpr double A = 10.0;
        constexpr double two_pi = 6.2831853071795864769;
        double xs = dx / cfg.rastrigin_scale;
        double ys = dy / cfg.rastrigin_scale;

        double R = 2.0 * A
                 + (xs * xs + ys * ys)
                 - A * std::cos(two_pi * xs)
                 - A * std::cos(two_pi * ys);

        return -R; // maximize
    }
    }
}