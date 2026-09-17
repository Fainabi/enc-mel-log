#include "openfhe.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
using namespace lbcrypto;
using C = std::complex<double>;
using CT = Ciphertext<DCRTPoly>;
static void show(const CryptoContext<DCRTPoly> &cc, const CT &ct,
                 const PrivateKey<DCRTPoly> &sk, const char *n) {
  Plaintext p;
  cc->Decrypt(ct, sk, &p);
  p->SetLength(2);
  auto v = p->GetCKKSPackedValue();
  printf("%s=(%.6g,%.6g),(%.6g,%.6g) level=%zu\n", n, v[0].real(), v[0].imag(),
         v[1].real(), v[1].imag(), ct->GetLevel());
}
static void inv2(CT &ct) {
  auto &v = ct->GetElements();
  auto ep = v[0].GetParams();
  std::vector<NativeInteger> a;
  for (auto &p : ep->GetParams())
    a.push_back(NativeInteger(2).ModInverse(p->GetModulus()));
  for (auto &c : v) {
    c.SetFormat(Format::COEFFICIENT);
    c = c.Times(a);
    c.SetFormat(Format::EVALUATION);
  }
}
static void mul_i(CT &ct) {
  auto &v = ct->GetElements();
  auto ep = v[0].GetParams();
  DCRTPoly md(ep, Format::COEFFICIENT, true);
  for (size_t k = 0; k < ep->GetParams().size(); k++) {
    NativePoly m(ep->GetParams()[k], Format::COEFFICIENT, true);
    m[ep->GetCyclotomicOrder() >> 2] = NativeInteger(1);
    md.SetElementAtIndex(k, std::move(m));
  }
  md.SetFormat(Format::EVALUATION);
  for (auto &c : v)
    c *= md;
}
int main() {
  CCParams<CryptoContextCKKSRNS> p;
  p.SetMultiplicativeDepth(20);
  p.SetScalingModSize(59);
  p.SetRingDim(1 << 16);
  p.SetSecurityLevel(HEStd_128_classic);
  p.SetScalingTechnique(FIXEDMANUAL);
  p.SetCKKSDataType(COMPLEX);
  auto cc = GenCryptoContext(p);
  for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE})
    cc->Enable(f);
  auto kp = cc->KeyGen();
  cc->EvalMultKeyGen(kp.secretKey);
  uint32_t cj = 2 * cc->GetRingDimension() - 1;
  auto ck = cc->EvalAutomorphismKeyGen(kp.secretKey, {cj});
  size_t S = cc->GetRingDimension() / 2;
  std::vector<C> x(S);
  x[0] = C(.1, .2);
  x[1] = C(.3, .4);
  auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x));
  show(cc, ct, kp.secretKey, "in");
  auto ti = ct;
  mul_i(ti);
  show(cc, ti, kp.secretKey, "mul_i");
  auto co = [&](CT z) { return cc->EvalAutomorphism(z, cj, *ck); };
  auto y = co(ct);
  auto re = cc->EvalAdd(ct, y);
  inv2(re);
  auto im = cc->EvalSub(ct, y);
  inv2(im);
  mul_i(im);
  for (auto &c : im->GetElements())
    c = c.Times(-1);
  show(cc, re, kp.secretKey, "re");
  show(cc, im, kp.secretKey, "im");
  auto r2 = cc->EvalSquare(re);
  auto i2 = cc->EvalSquare(im);
  cc->RescaleInPlace(r2);
  cc->RescaleInPlace(i2);
  show(cc, r2, kp.secretKey, "r2");
  show(cc, i2, kp.secretKey, "i2");
  auto i2i = cc->EvalMult(i2, C(0, 1));
  cc->RescaleInPlace(i2i);
  show(cc, i2i, kp.secretKey, "i2i");
  auto q = cc->EvalAdd(r2, i2i);
  show(cc, q, kp.secretKey, "q");
  return 0;
}
