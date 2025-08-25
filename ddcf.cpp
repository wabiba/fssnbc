#include "ddcf.h"
#include <random>
#include <chrono>

std::pair<DDCFKey, DDCFKey> GenDDCF(uint32_t n, uint32_t alpha, uint32_t beta1, uint32_t beta2) {
    uint32_t beta = (beta1 - beta2) % 16;

    //Fss fClient;
    //initializeClient(&fClient, n, 2);
    auto* fClient = new Fss();
    initializeClient(fClient, n, 2);

    ServerKeyLt k0, k1;
    memset(&k0, 0, sizeof(ServerKeyLt));
    memset(&k1, 0, sizeof(ServerKeyLt));
    generateTreeLt(fClient, &k0, &k1, alpha, beta);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, (1 << n) - 1);

    uint32_t S0 = dist(gen);
    uint32_t S1 = (-beta2 + (1 << n) + S0) % (1 << n);

    DDCFKey key0 = {k0, S0, fClient};
    DDCFKey key1 = {k1, S1, fClient};

    return {key0, key1};
}


uint32_t EvalDDCF(uint32_t b, DDCFKey& key, uint32_t x) {
    if (!key.fpk) {
        throw std::runtime_error("key.fpk is null");
    }
    if (!key.fpk->aes_keys) {
        throw std::runtime_error("key.fpk->aes_keys is null");
    }

    Fss fServer;
    initializeServer(&fServer, key.fpk);
    
    if (key.KeyPart.cw == nullptr) {
        std::cerr << "Error: key.KeyPart.cw is nullptr!" << std::endl;
    }
    if (fServer.aes_keys == nullptr) {
        std::cerr << "Error: fServer.aes_keys is nullptr!" << std::endl;
    }

    uint32_t yb = evaluateLt(&fServer, &key.KeyPart, x) % 16;

    return (b == 0) ? (yb + key.S) % 16 : (16 - yb - key.S) % 16;
}