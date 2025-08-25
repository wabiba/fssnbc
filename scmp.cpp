#include "scmp.h"
#include <random>

std::pair<SCMPKey, SCMPKey> GenSCMP(uint32_t n, uint32_t r1, uint32_t r2, uint32_t rout) {
    uint32_t maskN = (1u << n) - 1;
    uint32_t diff  = (r1 - r2) & maskN;          // (r1 - r2) mod 2^n
    uint32_t y     = ((1u << n) - diff) & maskN; // y ∈ U_{2^n}
    uint32_t y_n1  = (y >> (n - 1)) & 1;         // MSB
    uint32_t alpha = y & ((1u << (n - 1)) - 1);  // low (n-1) bits
    uint32_t beta1 = 1 ^ y_n1;
    uint32_t beta2 = y_n1;
    auto ddcf_keys = GenDDCF(n, alpha, beta1, beta2);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, (1 << n) - 1);

    uint32_t r_0 = dist(gen);
    uint32_t r_1 = (rout + (1 << n) - r_0) % (1 << n);

    SCMPKey key0 = {ddcf_keys.first, r_0};
    SCMPKey key1 = {ddcf_keys.second, r_1};
    
    return {key0, key1};
}

uint32_t EvalSCMP(uint32_t b, const SCMPKey& key, uint32_t x, uint32_t y) {
    uint32_t z = (x + 16 - y) % 16;
    uint32_t z_n1 = (z >> (4 - 1)) & 1;
    uint32_t z_n_minus_1 = (1u << (4 - 1)) - 1 - (z & ((1u << (4 - 1)) - 1));
    DDCFKey ddcf_copy = key.Keyddcf;
    uint32_t mb = EvalDDCF(b, ddcf_copy, z_n_minus_1);
    uint32_t res = (b - (b * z_n1 + mb - 2 * z_n1 * mb) + key.Rb) % 16;
    return res;
}
