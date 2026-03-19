#include "FitnessTaylor.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

FitnessTaylor::FitnessTaylor(int max_generations, int n0_, int nmax_)
    : G(max_generations), n0(n0_), nmax(nmax_) {
    if (G <= 0) throw std::invalid_argument("FitnessTaylor: max_generations must be > 0");
    if (n0 < 0) throw std::invalid_argument("FitnessTaylor: n0 must be >= 0");
    if (nmax < n0) throw std::invalid_argument("FitnessTaylor: nmax must be >= n0");
}

int FitnessTaylor::order_at(int g) const {
    // Smooth schedule: reaches nmax by the end of the run (G=1000)
    if (G <= 1) return nmax;
    g = std::clamp(g, 0, G - 1);
    double t = (double)g / (double)(G - 1);
    int n = (int)std::lround(n0 + t * (nmax - n0));
    return std::clamp(n, n0, nmax);
}

double FitnessTaylor::wrap_to_pi(double t) {
    constexpr double pi = 3.14159265358979323846;
    t = std::fmod(t, two_pi);
    if (t <= -pi) t += two_pi;
    else if (t > pi) t -= two_pi;
    return t;
}

double FitnessTaylor::cos_taylor_wrapped(double t, int n) {
    t = wrap_to_pi(t);

    double term = 1.0;
    double sum  = 1.0;
    double t2   = t * t;

    for (int k = 0; k < n; ++k) {
        double denom = (2.0 * k + 1.0) * (2.0 * k + 2.0);
        term *= (-t2) / denom;
        sum  += term;
    }
    return sum;
}

double FitnessTaylor::objective_taylor(int x, int y, int n) const {
    double xs = x / s;
    double ys = y / s;

    double cx = cos_taylor_wrapped(two_pi * xs, n);
    double cy = cos_taylor_wrapped(two_pi * ys, n);

    return 2.0 * A + (xs * xs + ys * ys) - A * cx - A * cy;
}

double FitnessTaylor::eval(const Chromosome& c, int g, std::mt19937& rng) {
    (void)rng;
    int x, y;
    decode_xy(c, x, y);

    int n = order_at(g);
    return -objective_taylor(x, y, n); // maximize
}