#include "Progress.h"
#include <cstdio>

void progress_bar(int done, int total, std::string_view label) {
    const int width = 30;
    if (total <= 0) total = 1;
    if (done < 0) done = 0;
    if (done > total) done = total;

    double frac = (double)done / (double)total;
    int filled = (int)(frac * width);

    std::printf("\r");
    if (!label.empty()) std::printf("%.*s ", (int)label.size(), label.data());

    std::printf("[");
    for (int i = 0; i < width; ++i) std::printf(i < filled ? "=" : " ");
    std::printf("] %3d%% (%d/%d)", (int)(frac * 100.0 + 0.5), done, total);
    std::fflush(stdout);

    if (done == total) std::printf("\n");
}