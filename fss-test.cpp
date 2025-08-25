#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include "fss-common.h"
#include "fss-server.h"
#include "fss-client.h"
#include "ddcf.h"
#include "scmp.h"
#include "mul.h"
#include <bits/stdc++.h>

using std::vector;
using std::array;
using std::pair;
using std::uint32_t;

struct RingCfg {
    uint32_t n = 16;
    uint32_t N = 1u << 16;
};
inline uint32_t modN(uint32_t x, const RingCfg& R){ return x % R.N; }

using RingVec = vector<uint32_t>;
using RingMat = vector<RingVec>;

static inline RingVec vadd(const RingVec& a, const RingVec& b, const RingCfg& R){
    RingVec c(a.size()); for(size_t i=0;i<a.size();++i) c[i]=modN(a[i]+b[i],R); return c;
}
static inline RingVec vsub(const RingVec& a, const RingVec& b, const RingCfg& R){
    RingVec c(a.size()); for(size_t i=0;i<a.size();++i) c[i]=modN(a[i]+R.N-b[i],R); return c;
}
static inline RingMat madd(const RingMat& A, const RingMat& B, const RingCfg& R){
    RingMat C=A; for(size_t r=0;r<A.size();++r) for(size_t c=0;c<A[r].size();++c) C[r][c]=modN(A[r][c]+B[r][c],R); return C;
}
static inline RingMat msub(const RingMat& A, const RingMat& B, const RingCfg& R){
    RingMat C=A; for(size_t r=0;r<A.size();++r) for(size_t c=0;c<A[r].size();++c) C[r][c]=modN(A[r][c]+R.N-B[r][c],R); return C;
}
static inline RingVec colDot(const RingVec& vec, const RingMat& M, const RingCfg& R){
    size_t rows = M.size(), cols = M[0].size();
    RingVec y(cols, 0);
    for(size_t j=0;j<cols;++j){
        uint32_t acc = 0;
        for(size_t i=0;i<rows;++i) acc = modN(acc + vec[i]*M[i][j], R);
        y[j] = acc;
    }
    return y;
}
static inline RingVec axpy(uint32_t a, const RingVec& y, const RingVec& z, const RingCfg& R){
    RingVec c(y.size()); for(size_t i=0;i<y.size();++i) c[i]=modN(a*y[i]+z[i],R); return c;
}

static inline uint32_t urand(const RingCfg& R){
    static std::random_device rd; static std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, R.N-1);
    return dist(gen);
}
static inline RingVec vrand(size_t m, const RingCfg& R){
    RingVec v(m); for(size_t i=0;i<m;++i) v[i]=urand(R); return v;
}
static inline RingMat mrand(size_t rows, size_t cols, const RingCfg& R){
    RingMat M(rows, RingVec(cols));
    for(size_t i=0;i<rows;++i) for(size_t j=0;j<cols;++j) M[i][j]=urand(R);
    return M;
}

struct CmpGateMasks { uint32_t r1, r2, rout; };
struct CmpGateKeys  { SCMPKey k0, k1; CmpGateMasks masks; };
struct MulGateMasks { uint32_t rin1, rin2, rout; };
struct MulGateKeys  { array<uint32_t,3> k0, k1; MulGateMasks masks; };

struct StepKeysPerI {
    CmpGateKeys  cmp;
    MulGateKeys  mul2, mul3, mul4, mul5;
};

struct KBundleForParty {
    vector<StepKeysPerI> stepKeys;
    RingVec t_b, d_b, q_b;
    RingVec o_b;
    RingMat H_b;
    uint32_t q_tilde_m_b;
};

struct TAOfflineAll {
    KBundleForParty K[2];
    RingVec r, q, r_tilde, q_tilde;
    vector<array<uint32_t,5>> L;
};

static TAOfflineAll TA_Offline(uint32_t m, uint32_t f, const RingCfg& R){
    TAOfflineAll out;
    out.L.resize(m);
    RingVec o0 = vrand(2*f, R);
    out.K[0].o_b = o0;
    out.K[1].o_b = vsub(RingVec(2*f,0), o0, R);
    RingVec d0 = vrand(m, R); out.K[0].d_b = d0; out.K[1].d_b = vsub(RingVec(m,0), d0, R);
    RingVec q0 = vrand(m, R); out.K[0].q_b = q0; out.K[1].q_b = vsub(RingVec(m,0), q0, R);
    out.r = vrand(m, R);
    out.q = vadd(out.K[0].q_b, out.K[1].q_b, R);
    RingMat H0 = mrand(2*f, m, R);
    out.K[0].H_b = H0;
    out.K[1].H_b = msub(RingMat(2*f, RingVec(m,0)), H0, R);
    RingVec o = vadd(out.K[0].o_b, out.K[1].o_b, R);
    RingMat H = madd(out.K[0].H_b, out.K[1].H_b, R);
    RingVec oH = colDot(o, H, R);
    RingVec t  = vadd(vadd(oH, vadd(out.K[0].d_b, out.K[1].d_b, R), R), out.r, R);
    RingVec t0 = vrand(m, R);
    RingVec t1 = vsub(t, t0, R);
    out.K[0].t_b = t0; out.K[1].t_b = t1;
    for(uint32_t i=0;i<m;++i) for(int j=0;j<5;++j) out.L[i][j]=urand(R);
    out.r_tilde.resize(m); out.q_tilde.resize(m);
    for(uint32_t i=0;i<m;++i){
        out.r_tilde[i] = modN(out.L[i][1] + R.N - out.L[i][2] + out.r[i], R);
        out.q_tilde[i] = modN(out.L[i][3] + R.N - out.L[i][4] + out.q[i], R);
    }
    out.K[0].stepKeys.resize(m);
    out.K[1].stepKeys.resize(m);
    for(uint32_t i=2;i<=m;++i){
        StepKeysPerI S;
        S.cmp.masks = {urand(R), urand(R), urand(R)};
        auto ck = GenSCMP(R.n, S.cmp.masks.r1, S.cmp.masks.r2, S.cmp.masks.rout);
        S.cmp.k0 = ck.first; S.cmp.k1 = ck.second;
        auto makeMul = [&](MulGateKeys& M){
            M.masks = {urand(R), urand(R), urand(R)};
            auto ks = GenX(M.masks.rin1, M.masks.rin2, M.masks.rout);
            M.k0 = ks.first; M.k1 = ks.second;
        };
        makeMul(S.mul2); makeMul(S.mul3); makeMul(S.mul4); makeMul(S.mul5);
        out.K[0].stepKeys[i-1] = S;
        out.K[1].stepKeys[i-1] = S;
    }
    out.K[0].q_tilde_m_b = urand(R);
    out.K[1].q_tilde_m_b = modN(out.q_tilde[m-1] + R.N - out.K[0].q_tilde_m_b, R);
    return out;
}

struct PartyInput {
    RingVec x_b;
    RingVec e_b;
    RingMat A_b;
    RingVec c_b;
};
struct PartyOutput { uint32_t R_b; };
struct OnlineState {
    RingVec x_hat;
    RingMat A_hat;
    RingVec e_hat;
    RingVec c_hat;
    RingVec p_hat;
    RingVec c_bar;
};

static RingVec compute_p_hat_side(uint32_t b,
                                  const RingVec& x_hat,
                                  const RingMat& A_hat,
                                  const RingVec& o_b,
                                  const RingVec& e_hat,
                                  const RingVec& t_b,
                                  const RingCfg& R)
{
    RingVec xA = colDot(x_hat, A_hat, R);
    RingVec oA = colDot(o_b,   A_hat, R);
    RingVec term1 = axpy(b, xA, RingVec(xA.size(),0), R);
    RingVec sum = vadd(term1, oA, R);
    sum = vadd(sum, e_hat, R);
    sum = vadd(sum, t_b, R);
    return sum;
}

static std::pair<PartyOutput,PartyOutput>
Online_Run(const RingCfg& R, uint32_t m, uint32_t f,
           const KBundleForParty& K0, const KBundleForParty& K1,
           const PartyInput& I0, const PartyInput& I1,
           OnlineState* dbg = nullptr)
{
    (void)f;
    assert(I0.x_b.size() == I1.x_b.size());
    assert(I0.A_b.size() == I1.A_b.size());
    assert(I0.A_b.empty() || I0.A_b[0].size() == I1.A_b[0].size());
    PartyOutput O0{0}, O1{0};
    OnlineState S;
    S.x_hat = vadd( vsub(I0.x_b, K0.o_b, R), vsub(I1.x_b, K1.o_b, R), R );
    S.e_hat = vadd( vsub(I0.e_b, K0.d_b, R), vsub(I1.e_b, K1.d_b, R), R );
    S.A_hat = madd( msub(I0.A_b, K0.H_b, R), msub(I1.A_b, K1.H_b, R), R );
    S.c_hat = vadd( vadd(I0.c_b, K0.q_b, R), vadd(I1.c_b, K1.q_b, R), R );
    S.c_bar = S.c_hat;
    RingVec p0 = compute_p_hat_side(0, S.x_hat, S.A_hat, K0.o_b, S.e_hat, K0.t_b, R);
    RingVec p1 = compute_p_hat_side(1, S.x_hat, S.A_hat, K1.o_b, S.e_hat, K1.t_b, R);
    S.p_hat = vadd(p0, p1, R);
    for(uint32_t i=1;i<m;++i){
        const auto& S0 = K0.stepKeys[i];
        const auto& S1 = K1.stepKeys[i];
        uint32_t ptil_prev = S.p_hat[i-1];
        uint32_t phat_i    = S.p_hat[i];
        uint32_t cbar_prev = S.c_bar[i-1];
        uint32_t cbar_i    = S.c_bar[i];
        auto a0 = modN(ptil_prev + S0.cmp.masks.r1, R);
        auto b0 = modN(phat_i    + S0.cmp.masks.r2, R);
        auto U10 = EvalSCMP(0, S0.cmp.k0, a0, b0);
        auto a1 = modN(ptil_prev + S1.cmp.masks.r1, R);
        auto b1 = modN(phat_i    + S1.cmp.masks.r2, R);
        auto U11 = EvalSCMP(1, S1.cmp.k1, a1, b1);
        uint32_t U1 = modN(U10 + U11, R);
        auto doMul = [&](const MulGateKeys& Kmul, uint32_t u, uint32_t v){
            uint32_t u_m = modN(u + Kmul.masks.rin1, R);
            uint32_t v_m = modN(v + Kmul.masks.rin2, R);
            uint32_t y0 = EvalX(0, Kmul.k0, u_m, v_m);
            uint32_t y1 = EvalX(1, Kmul.k1, u_m, v_m);
            return modN(y0 + y1, R);
        };
        uint32_t U2 = doMul(S0.mul2, U1, ptil_prev);
        uint32_t U3 = doMul(S0.mul3, U1, phat_i);
        uint32_t U4 = doMul(S0.mul4, U1, cbar_prev);
        uint32_t U5 = doMul(S0.mul5, U1, cbar_i);
        uint32_t pihat0_b = modN(U2 + R.N - U3 + 0*phat_i + 0*cbar_i, R);
        uint32_t pihat1_b = modN(U2 + R.N - U3 + 1*phat_i + 1*cbar_i, R);
        S.p_hat[i] = modN(pihat0_b + pihat1_b, R);
        (void)U4; (void)U5;
    }
    uint32_t cbar_m = S.c_bar[m-1];
    uint32_t half_c = (cbar_m >> 1);
    O0.R_b = modN(half_c + R.N - K0.q_tilde_m_b, R);
    O1.R_b = modN(half_c + R.N - K1.q_tilde_m_b, R);
    if(dbg) *dbg = S;
    return std::make_pair(O0, O1);
}

int main(){
    RingCfg R;
    uint32_t m = 10;
    uint32_t f = 10;
    auto TA = TA_Offline(m, f, R);
    PartyInput I0, I1;
    I0.x_b = vrand(2*f, R);     I1.x_b = vrand(2*f, R);
    I0.e_b = vrand(m, R);       I1.e_b = vrand(m, R);
    I0.A_b = mrand(2*f, m, R);  I1.A_b = mrand(2*f, m, R);
    I0.c_b = vrand(m, R);       I1.c_b = vrand(m, R);
    OnlineState dbg;
    auto start = std::chrono::high_resolution_clock::now();
    Online_Run(R, m, f, TA.K[0], TA.K[1], I0, I1, &dbg);
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Running Time: "
              << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count()
              << " us\n";
    return 0;
}
