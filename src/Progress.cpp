#include "Progress.h"
#include <cstdio>

// include your GA header here
#include "GA.h"   // adjust include path if needed

static void print_progress(int done, int total) {
    const int width = 30;
    double frac = (total == 0) ? 1.0 : (double)done / (double)total;
    if (frac < 0) frac = 0;
    if (frac > 1) frac = 1;

    int filled = (int)(frac * width);

    std::printf("\r[");
    for (int i = 0; i < width; ++i) std::printf(i < filled ? "=" : " ");
    std::printf("] %3d%% (%d/%d)", (int)(frac * 100.0 + 0.5), done, total);
    std::fflush(stdout);

    if (done == total) std::printf("\n");
}

void run_epochs(int epochs, int ans) {
    int hitcnt = 0;

    for (int e = 1; e <= epochs; ++e) {
        // Keep GA construction simple; you can swap these strings later.
        GA ga("fixed", "TournamentSelection", "RandomMutation", "SPC");

        ga.evaluate(); // runs full GA until termination

        if (ga.getAnswer() == ans) hitcnt++;

        print_progress(e, epochs);
    }

    std::printf("Success rate: %.4f (%d/%d)\n",
                hitcnt / (double)epochs, hitcnt, epochs);
}