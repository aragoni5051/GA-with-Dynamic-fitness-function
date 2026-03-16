#include "DecodeUtils.h"

int gray8_to_uint(const std::vector<bool>& bits, int offset) {
    int gray = 0;
    for (int i = 0; i < 8; ++i) {
        gray = (gray << 1) | (bits[offset + i] ? 1 : 0);
    }
    int bin = 0;
    for (; gray; gray >>= 1) bin ^= gray;   // Gray -> binary
    return bin; // 0..255
}

int gray8_to_signed_128_127(const std::vector<bool>& bits, int offset) {
    return gray8_to_uint(bits, offset) - 128; // -128..127
}

void decode_xy_gray_8_8(const std::vector<bool>& bits, int& x, int& y) {
    x = gray8_to_signed_128_127(bits, 0);
    y = gray8_to_signed_128_127(bits, 8);
}