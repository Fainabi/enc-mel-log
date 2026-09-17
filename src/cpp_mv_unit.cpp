#include "openfhe.h"
#include <complex>
#include <cstdio>
#include <vector>
using namespace lbcrypto;
using C = std::complex<double>;
using CT = Ciphertext<DCRTPoly>;
static constexpr size_t N = 512;
static CT mv(const CryptoContext<DCRTPoly> &cc, CT x,
             const std::vector<Plaintext> &p) {
  size_t n1 = 16;
  std::vector<CT> b(n1);
  b[0] = x;
  auto d = cc->EvalFastRotationPrecompute(x);
  for (size_t i = 1; i < n1; i++)
    b[i] = cc->EvalFastRotation(x, i, 2 * cc->GetRingDimension(), d);
  CT a;
  bool f = 1;
  for (size_t k = 0; k < N / n1; k++) {
    CT in;
    bool h = 0;
    for (size_t i = 0; i < n1; i++) {
      size_t z = k * n1 + i;
      auto t = cc->EvalMult(b[i], p[z]);
      if (!h) {
        in = t;
        h = 1;
      } else
        cc->EvalAddInPlace(in, t);
    }
    if (k)
      in = cc->EvalRotate(in, k * n1);
    if (f) {
      a = in;
      f = 0;
    } else
      cc->EvalAddInPlace(a, in);
  }
  cc->RescaleInPlace(a);
  return a;
}
int main() {
  CCParams<CryptoContextCKKSRNS> p;
  p.SetMultiplicativeDepth(4);
  p.SetScalingModSize(59);
  p.SetRingDim(1 << 16);
  p.SetSecurityLevel(HEStd_128_classic);
  p.SetScalingTechnique(FIXEDMANUAL);
  p.SetCKKSDataType(COMPLEX);
  auto cc = GenCryptoContext(p);
  for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE})
    cc->Enable(f);
  auto k = cc->KeyGen();
  std::vector<int32_t> rk;
  for (int i = 1; i < (int)N; i++)
    rk.push_back(i);
  cc->EvalRotateKeyGen(k.secretKey, rk);
  size_t S = cc->GetRingDimension() / 2;
  std::vector<C> x(S);
  x[0] = .2;
  x[1] = .3;
  x[17] = .7;
  auto ct = cc->Encrypt(k.publicKey, cc->MakeCKKSPackedPlaintext(x));
  std::vector<Plaintext> ps(N);
  for (size_t z = 0; z < N; z++) {
    std::vector<C> q(S);
    if (z == 0)
      q[1] = 1;
    if (z == 16)
      q[17] = 1;
    ps[z] = cc->MakeCKKSPackedPlaintext(q);
  }
  auto o = mv(cc, ct, ps);
  auto ro = cc->EvalRotate(ct, 1);
  Plaintext d;
  cc->Decrypt(o, k.secretKey, &d);
  d->SetLength(4);
  auto v = d->GetCKKSPackedValue();
  printf("out=(%.8g,%.8g),(%.8g,%.8g) level=%zu\n", v[0].real(), v[0].imag(),
         v[1].real(), v[1].imag(), o->GetLevel());
  Plaintext rd;
  cc->Decrypt(ro, k.secretKey, &rd);
  rd->SetLength(4);
  auto rv = rd->GetCKKSPackedValue();
  printf("rotate=(%.8g,%.8g),(%.8g,%.8g)\n", rv[0].real(), rv[0].imag(),
         rv[1].real(), rv[1].imag());
  std::vector<std::vector<C>> D(N, std::vector<C>(N));
  for (size_t k0 = 0; k0 < 13; k0++)
    for (size_t n = 0; n < 80; n++)
      D[k0][n] = std::cos(M_PI * k0 * (n + .5) / 80);
  std::vector<Plaintext> dp(N);
  for (size_t z = 0; z < N; z++) {
    std::vector<C> q(S);
    for (size_t r = 0; r < N; r++)
      q[r] = D[r][(r + z) % N];
    dp[z] = cc->MakeCKKSPackedPlaintext(q);
  }
  auto od = mv(cc, ct, dp);
  Plaintext dd;
  cc->Decrypt(od, k.secretKey, &dd);
  dd->SetLength(13);
  auto dv = dd->GetCKKSPackedValue();
  printf("dct unit=(%.8g,%.8g),(%.8g,%.8g)\n", dv[0].real(), dv[0].imag(),
         dv[1].real(), dv[1].imag());
}
