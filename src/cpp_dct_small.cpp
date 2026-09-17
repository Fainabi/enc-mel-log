#include "openfhe.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
using namespace lbcrypto;
using C = std::complex<double>;
using CT = Ciphertext<DCRTPoly>;
static constexpr size_t N = 512, S = 32768;

static CT mv(const CryptoContext<DCRTPoly>& cc, CT x,
             const std::vector<Plaintext>& p) {
  constexpr size_t n1 = 16;
  std::vector<CT> baby(n1); baby[0] = x;
  auto dig = cc->EvalFastRotationPrecompute(x);
  for (size_t i = 1; i < n1; ++i)
    baby[i] = cc->EvalFastRotation(x, i, 2 * cc->GetRingDimension(), dig);
  CT acc; bool first = true;
  for (size_t k = 0; k < N / n1; ++k) {
    CT inner; bool have = false;
    for (size_t i = 0; i < n1; ++i) {
      size_t d = k * n1 + i;
      CT t = cc->EvalMult(baby[i], p[d]);
      if (!have) { inner = t; have = true; }
      else cc->EvalAddInPlace(inner, t);
    }
    if (k) inner = cc->EvalRotate(inner, static_cast<int>(k * n1));
    if (first) { acc = inner; first = false; }
    else cc->EvalAddInPlace(acc, inner);
  }
  cc->RescaleInPlace(acc);
  return acc;
}
static CT mv_direct(const CryptoContext<DCRTPoly>& cc, CT x,
                    const std::vector<Plaintext>& p,
                    const std::vector<int>& ds) {
  CT acc; bool first = true;
  for (int d : ds) {
    CT rot = d ? cc->EvalRotate(x, d) : x;
    CT term = cc->EvalMult(rot, p[d + 79]);
    if (first) { acc = term; first = false; }
    else cc->EvalAddInPlace(acc, term);
  }
  cc->RescaleInPlace(acc);
  return acc;
}

int main() {
  CCParams<CryptoContextCKKSRNS> p;
  p.SetMultiplicativeDepth(20); p.SetScalingModSize(59);
  p.SetRingDim(1 << 16); p.SetSecurityLevel(HEStd_128_classic);
  p.SetScalingTechnique(FIXEDMANUAL); p.SetCKKSDataType(COMPLEX);
  auto cc = GenCryptoContext(p);
  for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE}) cc->Enable(f);
  auto kp = cc->KeyGen(); cc->EvalMultKeyGen(kp.secretKey);
  std::vector<int32_t> keys;
  for (int i = 1; i <= 79; ++i) { keys.push_back(i); keys.push_back(-i); }
  cc->EvalRotateKeyGen(kp.secretKey, keys);

  std::vector<C> x(S);
  for (size_t n = 0; n < 80; ++n)
    x[n] = C(0.094 + 0.18 * std::exp(-0.12 * n),
             0.094 + 0.16 * std::exp(-0.12 * n));
  auto ct0 = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x));

  std::vector<std::vector<C>> D(N, std::vector<C>(N));
  for (size_t k = 0; k < 13; ++k)
    for (size_t n = 0; n < 80; ++n) D[k][n] = std::cos(M_PI * k * (n + .5) / 80);
  std::vector<Plaintext> dp_direct(159);
  std::vector<int> signed_ds;
  for (int d = -79; d <= 79; ++d) {
    std::vector<C> q(S);
    for (size_t r = 0; r < N; ++r) {
      int n = static_cast<int>(r) + d;
      if (n >= 0 && n < static_cast<int>(N)) q[r] = D[r][n];
    }
    dp_direct[d + 79] = cc->MakeCKKSPackedPlaintext(q);
    signed_ds.push_back(d);
  }
  std::vector<C> expected(N);
  for (size_t k = 0; k < 13; ++k)
    for (size_t n = 0; n < 80; ++n) expected[k] += D[k][n] * x[n];

  for (int depth_case = 0; depth_case < 2; ++depth_case) {
    CT ct = ct0;
    if (depth_case) {
      std::vector<C> one(S, C(1.0, 0.0));
      auto ptone = cc->MakeCKKSPackedPlaintext(one);
      for (int i = 0; i < 10; ++i) {
        ct = cc->EvalMult(ct, ptone);
        cc->RescaleInPlace(ct);
      }
    }
    auto out = mv_direct(cc, ct, dp_direct, signed_ds);
    Plaintext dec; cc->Decrypt(out, kp.secretKey, &dec); dec->SetLength(13);
    auto v = dec->GetCKKSPackedValue();
    double se = 0, me = 0;
    for (size_t k = 0; k < 13; ++k) {
      double e = std::abs(v[k] - expected[k]); se += e * e; me = std::max(me, e);
    }
    printf("small_dct case=%s input_level=%zu output_level=%zu rmse=%.6e max=%.6e s0=(%.6g,%.6g) exp=(%.6g,%.6g)\n",
           depth_case ? "deep" : "fresh", ct->GetLevel(), out->GetLevel(),
           std::sqrt(se / 13), me, v[0].real(), v[0].imag(),
           expected[0].real(), expected[0].imag());
    for (size_t k = 0; k < 13; ++k)
      printf("  k%zu got=(%.7g,%.7g) exp=(%.7g,%.7g)\n", k,
             v[k].real(), v[k].imag(), expected[k].real(), expected[k].imag());
  }
}
