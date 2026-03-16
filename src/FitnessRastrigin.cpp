#include "FitnessRastrigin.h"
#include "DecodeUtils.h"
#include <cmath>

static double rastrigin2d_scaled(int x, int y) {
    constexpr double A = 10.0;
    constexpr double s = 25.0; // 128/25 = 5.12 (standard rastrigin domain)
    constexpr double two_pi = 6.2831853071795864769;

    double xs = x / s;
    double ys = y / s;

    return 2.0 * A
         + (xs*xs + ys*ys)
         - A * std::cos(two_pi * xs)
         - A * std::cos(two_pi * ys);
}

double FitnessRastrigin::eval(const Chromosome& c, int /*gen*/, std::mt19937& /*rng*/) {
    int x, y;
    decode_xy_gray_8_8(c.bits, x, y);

    // Rastrigin is a minimization function; GA usually maximizes -> negate it.
    return -rastrigin2d_scaled(x, y); // global max at (0,0)
}