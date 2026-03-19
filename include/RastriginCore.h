#pragma once
#include "Chromosome.h"
#include "DecodeUtils.h"
#include <cmath>

struct RastriginCore {
protected:
    // constants shared by all variants
    static constexpr double A = 10.0;
    static constexpr double s = 25.0; // scale: x in [-128..127] -> xs in ~[-5.12..5.08]
    static constexpr double two_pi = 6.2831853071795864769;

    static inline void decode_xy(const Chromosome& c, int& x, int& y) {
        decode_xy_gray_8_8(c.bits, x, y);
    }

    // True objective (minimize): R(x,y)
    static inline double objective_true(int x, int y) {
        double xs = x / s;
        double ys = y / s;
        return 2.0 * A
             + (xs * xs + ys * ys)
             - A * std::cos(two_pi * xs)
             - A * std::cos(two_pi * ys);
    }

    // Base fitness for GA (maximize): f_base = -R
    static inline double base_fitness(int x, int y) {
        return -objective_true(x, y);
    }

    // Convenience from chromosome
    static inline double base_fitness(const Chromosome& c) {
        int x, y;
        decode_xy(c, x, y);
        return base_fitness(x, y);
    }
};