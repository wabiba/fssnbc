#ifndef SCMP_H
#define SCMP_H

#include <iostream>
#include <vector>
#include <utility>
#include <random>
#include "ddcf.h"

struct SCMPKey {
    DDCFKey Keyddcf;
    uint32_t Rb;
};

std::pair<SCMPKey, SCMPKey> GenSCMP(uint32_t n, uint32_t r1, uint32_t r2, uint32_t rout);

uint32_t EvalSCMP(uint32_t b, const SCMPKey& key, uint32_t x, uint32_t y);

#endif
