#include "openfhe.h"
#include "openfhe/core/math/chebyshev.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
using namespace lbcrypto;
using C = std::complex<double>;
using CT = Ciphertext<DCRTPoly>;
static constexpr size_t N = 512, S = 1 << 15;

int main() {
  CCParams<CryptoContextCKKSRNS> p;
  p.SetMultiplicativeDepth(20);
  p.SetScalingModSize(59);
  p.SetRingDim(1 << 16);
  p.SetSecurityLevel(HEStd_128_classic);
  p.SetScalingTechnique(FIXEDMANUAL);
  p.SetCKKSDataType(REAL);
  auto cc = GenCryptoContext(p);
  for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE})
    cc->Enable(f);
  auto kp = cc->KeyGen();
  cc->EvalMultKeyGen(kp.secretKey);
  std::vector<int32_t> keys;
  for (int i = 1; i < (int)N; i++)
    keys.push_back(i);
  cc->EvalRotateKeyGen(kp.secretKey, keys);
  std::vector<double> x(S);
  for (size_t i = 0; i < 80; i++)
    x[i] = .08 + .0002 * i;
  auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x));
  auto fn = [](double z) { return std::pow(std::max(z, 0.0), .25); };
  auto z = cc->EvalChebyshevFunction(fn, ct, 0, 1, 63);
  std::vector<Plaintext> d(N);
  for (size_t j = 0; j < N; j++) {
    std::vector<double> q(S);
    for (size_t r = 0; r < N; r++)
      q[r] = std::cos(M_PI * r * ((r + j) % N + .5) / 80.0) *
             ((r < 13 && (r + j) % N < 80) ? 1.0 : 0.0);
    d[j] = cc->MakeCKKSPackedPlaintext(q);
  }
  CT out;
  bool first = true;
  for (size_t j = 0; j < N; j++) {
    bool nz = false;
    for (size_t r = 0; r < 13; r++)
      nz = nz ||
           std::abs(std::cos(M_PI * r * (((r + j) % N) + .5) / 80.0)) > 1e-15 &&
               (r + j) % N < 80;
    if (!nz)
      continue;
    auto rot = j ? cc->EvalRotate(z, (int)j) : z;
    auto t = cc->EvalMult(rot, d[j]);
    if (first) {
      out = t;
      first = false;
    } else
      cc->EvalAddInPlace(out, t);
  }
  cc->RescaleInPlace(out);
  Plaintext dec;
  cc->Decrypt(out, kp.secretKey, &dec);
  dec->SetLength(13);
  auto v = dec->GetCKKSPackedValue();
  std::vector<double> ev(13);
  for (size_t r = 0; r < 13; r++)
    for (size_t n = 0; n < 80; n++)
      ev[r] += std::cos(M_PI * r * (n + .5) / 80.0) *
               std::pow(std::max(x[n], 0.0), .25);
  double se = 0, me = 0;
  for (size_t r = 0; r < 13; r++) {
    double e = std::abs(v[r].real() - ev[r]);
    se += e * e;
    me = std::max(me, e);
  }
  printf("real_dct level=%zu rmse=%.6e maxerr=%.6e s0=%.8g exp=%.8g s1=%.8g "
         "exp=%.8g\n",
         out->GetLevel(), std::sqrt(se / 13), me, v[0].real(), ev[0],
         v[1].real(), ev[1]);
}
