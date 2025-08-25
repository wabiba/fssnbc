#ifndef MUL_H
#define MUL_H

#include <array>
#include <cstdint>

constexpr uint32_t N = 1 << 16;

std::pair<std::array<uint32_t, 3>, std::array<uint32_t, 3>> GenX(uint32_t rin1, uint32_t rin2, uint32_t rout);

uint32_t EvalX(int b, const std::array<uint32_t, 3>& kb, uint32_t x1, uint32_t x2);

#endif
