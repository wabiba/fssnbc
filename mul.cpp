#include "mul.h"
#include <iostream>
#include <random>

std::pair<std::array<uint32_t, 3>, std::array<uint32_t, 3>> GenX(uint32_t rin1, uint32_t rin2, uint32_t rout) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, N - 1);

    uint32_t P0 = dist(gen);
    uint32_t P1 = (rin1 - P0) % N;

    uint32_t Q0 = dist(gen);
    uint32_t Q1 = (rin2 - Q0) % N;

    uint32_t R0 = dist(gen);
    uint32_t R1 = (rin1 * rin2 + rout - R0) % N;

    std::array<uint32_t, 3> K0 = {P0, Q0, R0};
    std::array<uint32_t, 3> K1 = {P1, Q1, R1};

    return {K0, K1};
}

uint32_t EvalX(int b, const std::array<uint32_t, 3>& kb, uint32_t x1, uint32_t x2) {
    uint32_t Pb = kb[0];
    uint32_t Qb = kb[1];
    uint32_t Rb = kb[2];

    return (static_cast<uint32_t>(b) * (x1 * x2) - x1 * Qb - x2 * Pb + Rb) % N;
}
