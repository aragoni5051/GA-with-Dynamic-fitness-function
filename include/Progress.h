#pragma once
#include <string_view>

// Updates a single-line progress bar (overwrites the same line).
// Call repeatedly with done in [0,total]. Prints newline automatically when done==total.
void progress_bar(int done, int total, std::string_view label = {});