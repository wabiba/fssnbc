#ifndef DDCF_H
#define DDCF_H

#include <iostream>
#include <vector>
#include <utility>
#include "fss-client.h"
#include "fss-server.h"

struct DDCFKey {
    ServerKeyLt KeyPart;
    uint32_t S;     
    Fss* fpk; 
};

std::pair<DDCFKey, DDCFKey> GenDDCF(uint32_t n, uint32_t alpha, uint32_t beta1, uint32_t beta2);

uint32_t EvalDDCF(uint32_t b, DDCFKey& key, uint32_t x);

#endif
