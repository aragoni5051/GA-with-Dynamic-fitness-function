#pragma once

// Runs multiple GA epochs, prints a progress bar, and prints final success rate.
// - epochs: number of independent runs
// - ans: target answer to count as "success"
void run_epochs(int epochs, int ans);