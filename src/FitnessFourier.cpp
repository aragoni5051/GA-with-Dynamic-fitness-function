#include "FitnessFourier.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

FitnessFourier::FitnessFourier(int max_generations, int n0_, int nmax_)
    : G(max_generations), n0(n0_), nmax(nmax_) {
    if (G <= 0) throw std::invalid_argument("FitnessFourier: max_generations must be > 0");
    if (n0 < 1) throw std::invalid_argument("FitnessFourier: n0 must be >= 1");
    if (nmax < n0) throw std::invalid_argument("FitnessFourier: nmax must be >= n0");

    amp.assign(nmax + 1, 0.0);
    for (int k = 1; k <= nmax; ++k) amp[k] = 10.0 / (double)k; // tunable
}

int FitnessFourier::order_at(int g) const {
    if (G <= 1) return nmax;
    g = std::clamp(g, 0, G - 1);
    double t = (double)g / (double)(G - 1);
    int n = (int)std::lround(n0 + t * (nmax - n0));
    return std::clamp(n, n0, nmax);
}

double FitnessFourier::objective_fourier(int x, int y, int n) const {
    double xs = x / s;
    double ys = y / s;

    // base bowl (unique min at 0)
    double val = xs*xs + ys*ys;

    // add harmonics: (1 - cos(...)) is >=0 and equals 0 at origin
    for (int k = 1; k <= n; ++k) {
        val += amp[k] * ( (1.0 - std::cos(two_pi * k * xs)) + (1.0 - std::cos(two_pi * k * ys)) );
    }
    return val;
}

double FitnessFourier::eval(const Chromosome& c, int g, std::mt19937& rng) {
    (void)rng;
    int x, y;
    decode_xy(c, x, y);

    int n = order_at(g);
    return -objective_fourier(x, y, n); // maximize
}