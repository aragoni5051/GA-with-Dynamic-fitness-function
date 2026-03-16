#pragma once
#include <vector>

int gray8_to_uint(const std::vector<bool>& bits, int offset);
int gray8_to_signed_128_127(const std::vector<bool>& bits, int offset);
void decode_xy_gray_8_8(const std::vector<bool>& bits, int& x, int& y);