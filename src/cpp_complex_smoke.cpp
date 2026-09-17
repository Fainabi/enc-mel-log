#include "openfhe.h"
#include <complex>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>

using namespace lbcrypto;

int main() {
    using C = std::complex<double>;
    CCParams<CryptoContextCKKSRNS> p;
    p.SetMultiplicativeDepth(5);
    p.SetScalingModSize(50);
    p.SetSecurityLevel(HEStd_128_classic);
    p.SetRingDim(1 << 15);
    p.SetScalingTechnique(FIXEDMANUAL);
    p.SetCKKSDataType(COMPLEX);
    auto cc = GenCryptoContext(p);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE); cc->Enable(ADVANCEDSHE);
    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);

    const size_t slots = cc->GetRingDimension() / 2;
    std::vector<C> x(slots), m(slots), expected(slots);
    for (size_t i = 0; i < slots; ++i) {
        x[i] = C(0.001 * std::sin(0.01 * i), 0.001 * std::cos(0.013 * i));
        m[i] = C(0.7 + 0.1 * std::sin(0.007 * i), 0.2 * std::cos(0.011 * i));
        expected[i] = x[i] * m[i];
    }
    auto pt = cc->MakeCKKSPackedPlaintext(m);
    std::cout << "plaintext m[0]=" << pt->GetCKKSPackedValue()[0] << "\n";
    auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x));
    Plaintext in_dec; cc->Decrypt(ct, kp.secretKey, &in_dec); in_dec->SetLength(slots);
    auto in_got = in_dec->GetCKKSPackedValue();
    std::printf("input got=(%.6g,%.6g) expected=(%.6g,%.6g)\n", in_got[0].real(), in_got[0].imag(), x[0].real(), x[0].imag());
    auto out = cc->EvalMult(ct, pt);
    cc->RescaleInPlace(out);
    Plaintext dec;
    cc->Decrypt(out, kp.secretKey, &dec);
    dec->SetLength(slots);
    auto got = dec->GetCKKSPackedValue();
    double maxerr = 0.0;
    for (size_t i = 0; i < slots; ++i) maxerr = std::max(maxerr, std::abs(got[i] - expected[i]));
    std::printf("sample got=(%.6g,%.6g) expected=(%.6g,%.6g)\n", got[0].real(), got[0].imag(), expected[0].real(), expected[0].imag());
    std::printf("slots=%zu max_complex_mult_error=%.6e\n", slots, maxerr);
    return maxerr < 1e-5 ? 0 : 1;
}
