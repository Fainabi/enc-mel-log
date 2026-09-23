#include "openfhe.h"
#include "math/chebyshev.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace lbcrypto;
using C = std::complex<double>;
using CT = Ciphertext<DCRTPoly>;
#ifdef E52_SMALL
static constexpr size_t N = 256, WIN = 200, BLOCK = 1024, FRAMES = 16;
#else
static constexpr size_t N = 512, WIN = 400, BLOCK = 1024, FRAMES = 64;
#endif
static constexpr size_t BLOCKS = FRAMES / 2;
static constexpr size_t MV_N1 = 16;
static constexpr size_t LANES = FRAMES;
static constexpr size_t ROT_STRIDE = LANES;
// The default synthetic input produces normalized mel values well below one.
// Real-data calibration can override this public normalization at runtime.
static const double MEL_NORM = [] {
  const char* s = std::getenv("E52_MEL_NORM");
  return s ? std::strtod(s, nullptr) : 1.0e-2;
}();
static constexpr double CHEB_A = 0.0;
static constexpr double CHEB_B = 1.0;
static const double POWER_ALPHA = [] {
  const char* s = std::getenv("E52_POWER_ALPHA");
  return s ? std::strtod(s, nullptr) : 0.75;
}();
static const double POWER_OFFSET = [] {
  const char* s = std::getenv("E52_POWER_OFFSET");
  return s ? std::strtod(s, nullptr) : 0.0;
}();
static const uint32_t CHEB_DEGREE = [] {
  const char* s = std::getenv("E52_CHEB_DEGREE");
  return s ? static_cast<uint32_t>(std::strtoul(s, nullptr, 10)) : 63u;
}();
static const bool ZERO_CENTER = std::getenv("E52_ZERO_CENTER") != nullptr;
static constexpr double INPUT_GAIN = 1.0;
static std::vector<std::vector<C>> proj() {
  std::vector<std::vector<C>> M(N, std::vector<C>(N));
  for (size_t j = 0; j < WIN; j++) {
    double w = .5 - .5 * cos(2 * M_PI * j / WIN);
    M[0][j] = w;
    M[N - 1][j] = cos(M_PI * j) * w;
    for (size_t k = 1; k < N / 2; k++) {
      M[2 * k - 1][j] = cos(2 * M_PI * k * j / N) * w;
      M[2 * k][j] = -sin(2 * M_PI * k * j / N) * w;
    }
  }
  return M;
}
static std::vector<std::vector<C>> mm(const std::vector<std::vector<C>>& A,
                                      const std::vector<std::vector<C>>& B) {
  std::vector<std::vector<C>> R(N, std::vector<C>(N));
  for (size_t i = 0; i < N; ++i)
    for (size_t k = 0; k < N; ++k) if (A[i][k] != C(0))
      for (size_t j = 0; j < N; ++j) R[i][j] += A[i][k] * B[k][j];
  return R;
}
static std::vector<std::vector<std::vector<C>>> fft3() {
  std::vector<std::vector<std::vector<C>>> stages;
  for (size_t m = N; m >= 2; m >>= 1) {
    std::vector<std::vector<C>> S(N, std::vector<C>(N));
    size_t h = m / 2;
    for (size_t blk = 0; blk < N; blk += m)
      for (size_t j = 0; j < h; ++j) {
        size_t a = blk + j, b = a + h;
        C w = std::exp(C(0, -2.0 * M_PI * j / m));
        S[a][a] = S[a][b] = 1.0;
        S[b][a] = w; S[b][b] = -w;
      }
    stages.push_back(std::move(S));
  }
  std::vector<std::vector<std::vector<C>>> out;
  std::vector<std::vector<C>> cur(N, std::vector<C>(N));
  for (size_t i = 0; i < N; ++i) cur[i][i] = 1.0;
  for (size_t g = 0; g < stages.size(); g += 3) {
    std::vector<std::vector<C>> M(N, std::vector<C>(N));
    for (size_t i = 0; i < N; ++i) M[i][i] = 1.0;
    for (size_t s = 0; s < 3 && g + s < stages.size(); ++s)
      M = mm(stages[g + s], M);
    if (g == 0)
      for (size_t j = 0; j < WIN; ++j) {
        double w = .5 - .5 * cos(2 * M_PI * j / WIN);
        for (size_t i = 0; i < N; ++i) M[i][j] *= w;
      }
    out.push_back(std::move(M));
  }
  std::vector<std::vector<C>> br(N, std::vector<C>(N));
  size_t bits = 0;
  for (size_t m = N; m > 1; m >>= 1) ++bits;
  for (size_t k = 0; k < N; ++k) {
    size_t q = 0;
    for (size_t b = 0; b < bits; ++b) q = (q << 1) | ((k >> b) & 1);
    br[k][q] = 1.0;
  }
  if (std::getenv("E52_REAL_INPUT"))
    return out;
  std::vector<std::vector<C>> pack(N, std::vector<C>(N));
  pack[0][0] = 1.0;
  for (size_t k = 1; k < N / 2; ++k) {
    pack[2*k-1][k] = .5; pack[2*k-1][N-k] = .5;
    pack[2*k][k] = C(0,-.5); pack[2*k][N-k] = C(0,.5);
  }
  pack[N-1][N/2] = 1.0;
  out.push_back(mm(pack, br));
  return out;
}
static std::vector<std::vector<C>> mel() {
  std::vector<std::vector<C>> M(N, std::vector<C>(N));
  const bool real_input = std::getenv("E52_REAL_INPUT") != nullptr;
  auto bit_reverse = [](size_t k) {
    size_t q = 0;
    size_t bits = 0;
    for (size_t m = N; m > 1; m >>= 1) ++bits;
    for (size_t b = 0; b < bits; ++b)
      q = (q << 1) | ((k >> b) & 1);
    return q;
  };
  for (int r = 0; r < 80; r++) {
    int lo = 1 + r * 3, mid = lo + 3, hi = mid + 3;
    for (int k = lo; k < hi && k < static_cast<int>(N); k++)
      M[r][real_input ? bit_reverse(static_cast<size_t>(k)) : k] =
          (k <= mid) ? double(k - lo) / (mid - lo)
                     : double(hi - k) / (hi - mid);
  }
  return M;
}
static std::vector<std::vector<C>> dct() {
  std::vector<std::vector<C>> M(N, std::vector<C>(N));
  for (int k = 0; k < 13; k++)
    for (int n = 0; n < 80; n++)
      M[k][n] = cos(M_PI * k * (n + .5) / 80);
  return M;
}
static void probe_rmse(const CryptoContext<DCRTPoly> &cc, const CT &ct,
                       const PrivateKey<DCRTPoly> &sk,
                       const std::vector<C> &expected, const char *name) {
  Plaintext p;
  cc->Decrypt(ct, sk, &p);
  p->SetLength(cc->GetRingDimension() / 2);
  auto v = p->GetCKKSPackedValue();
  double se = 0.0, me = 0.0;
  for (size_t i = 0; i < expected.size(); ++i) {
    double e = std::abs(v[i * LANES] - expected[i]);
    se += e * e;
    me = std::max(me, e);
  }
  printf("stage %-12s level=%zu rmse=%.6e maxerr=%.6e\n", name,
         ct->GetLevel(), std::sqrt(se / expected.size()), me);
}
static std::vector<C> decrypt_values(const CryptoContext<DCRTPoly> &cc,
                                      const CT &ct,
                                      const PrivateKey<DCRTPoly> &sk) {
  Plaintext p;
  cc->Decrypt(ct, sk, &p);
  p->SetLength(cc->GetRingDimension() / 2);
  auto all = p->GetCKKSPackedValue();
  std::vector<C> lane(N);
  for (size_t i = 0; i < N; ++i)
    lane[i] = all[i * LANES];
  return lane;
}
static void print_stats(const char *name, const std::vector<C> &v) {
  double minr = 1e300, maxr = -1e300, mean = 0, ss = 0;
  size_t neg = 0, over = 0;
  for (const auto &z : v) {
    double a = z.real();
    minr = std::min(minr, a); maxr = std::max(maxr, a);
    mean += a; neg += (a < 0); over += (a > 1);
  }
  mean /= v.size();
  for (const auto &z : v) { double d = z.real() - mean; ss += d * d; }
  printf("input %-8s min=%.6e max=%.6e mean=%.6e std=%.6e neg=%zu over1=%zu\n",
         name, minr, maxr, mean, std::sqrt(ss / v.size()), neg, over);
}
static CT mv(const CryptoContext<DCRTPoly> &cc, CT x,
             const std::vector<Plaintext> &pts) {
  const size_t n1 = MV_N1;
  std::vector<bool> used_i(n1, false);
  for (size_t d = 0; d < N; ++d)
    if (pts[d]) used_i[d % n1] = true;
  std::vector<CT> baby(n1);
  baby[0] = x;
  auto dig = cc->EvalFastRotationPrecompute(x);
  for (size_t i = 1; i < n1; i++)
    if (used_i[i])
      baby[i] = cc->EvalFastRotation(x, i * ROT_STRIDE,
                                     2 * cc->GetRingDimension(), dig);
  CT acc;
  bool first = true;
  for (size_t k = 0; k < (N + n1 - 1) / n1; k++) {
      CT inner;
      bool have = false;
      for (size_t i = 0; i < n1; i++) {
        size_t d = k * n1 + i;
        if (d >= N || !pts[d]) continue;
        CT t = cc->EvalMult(baby[i], pts[d]);
      if (!have) {
        inner = t;
        have = true;
      } else
        cc->EvalAddInPlace(inner, t);
    }
    if (!have)
      continue;
    if (k)
      inner = cc->EvalRotate(inner, (int)(k * n1 * ROT_STRIDE));
    if (first) {
      acc = inner;
      first = false;
    } else
      cc->EvalAddInPlace(acc, inner);
  }
  cc->RescaleInPlace(acc);
  return acc;
}
static void mul_i_inplace(CT &ct) {
  auto &cv = ct->GetElements();
  auto ep = cv[0].GetParams();
  DCRTPoly md(ep, Format::COEFFICIENT, true);
  for (size_t k = 0; k < ep->GetParams().size(); k++) {
    NativePoly mon(ep->GetParams()[k], Format::COEFFICIENT, true);
    mon[ep->GetCyclotomicOrder() >> 2] = NativeInteger(1);
    md.SetElementAtIndex(k, std::move(mon));
  }
  md.SetFormat(Format::EVALUATION);
  for (auto &c : cv)
    c *= md;
}
int main() {
  CCParams<CryptoContextCKKSRNS> p;
  p.SetMultiplicativeDepth(20);
  p.SetScalingModSize(59);
#ifdef E52_SMALL
  p.SetRingDim(1 << 15);
#else
  p.SetRingDim(1 << 16);
#endif
#ifdef E52_SMALL
  p.SetSecurityLevel(HEStd_NotSet);
#else
  p.SetSecurityLevel(HEStd_128_classic);
#endif
  p.SetScalingTechnique(FIXEDMANUAL);
  p.SetCKKSDataType(COMPLEX);
  auto cc = GenCryptoContext(p);
  for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE})
    cc->Enable(f);
  auto kp = cc->KeyGen();
  cc->EvalMultKeyGen(kp.secretKey);
  std::vector<int32_t> rk;
  auto add_rk = [&](int32_t r) {
    if (r != 0 && std::find(rk.begin(), rk.end(), r) == rk.end())
      rk.push_back(r);
  };
  // BSGS uses baby steps 1..MV_N1-1 and giant steps that are multiples of
  // MV_N1.  DCT additionally uses its signed direct-rotation offsets.
  for (int32_t i = 1; i < static_cast<int32_t>(MV_N1); ++i)
    add_rk(i * static_cast<int32_t>(ROT_STRIDE));
  for (int32_t i = static_cast<int32_t>(MV_N1); i < static_cast<int32_t>(N);
       i += static_cast<int32_t>(MV_N1))
    add_rk(i * static_cast<int32_t>(ROT_STRIDE));
  add_rk(-80 * static_cast<int32_t>(ROT_STRIDE));
  cc->EvalRotateKeyGen(kp.secretKey, rk);
  uint32_t cj = 2 * cc->GetRingDimension() - 1;
  auto ck = cc->EvalAutomorphismKeyGen(kp.secretKey, {cj});
  size_t S = cc->GetRingDimension() / 2;
  const bool bench = std::getenv("E52_BENCH") != nullptr;
  const char* input_bin = std::getenv("E52_INPUT_BIN");
  std::vector<C> x(S);
  if (input_bin) {
    std::ifstream in(input_bin, std::ios::binary);
    in.seekg(0, std::ios::end);
    const std::streamoff bytes = in.tellg();
    in.seekg(0, std::ios::beg);
    const size_t real_count = WIN * LANES;
    const size_t complex_count = 2 * real_count;
    if (bytes != static_cast<std::streamoff>(real_count * sizeof(float)) &&
        bytes != static_cast<std::streamoff>(complex_count * sizeof(float))) {
      std::fprintf(stderr, "failed to read E52_INPUT_BIN=%s (%zu or %zu float32 values expected)\n",
                   input_bin, real_count, complex_count);
      return 2;
    }
    const bool complex_input =
        bytes == static_cast<std::streamoff>(complex_count * sizeof(float));
    std::vector<float> raw(complex_input ? complex_count : real_count);
    in.read(reinterpret_cast<char*>(raw.data()),
            static_cast<std::streamsize>(raw.size() * sizeof(float)));
    if (!in) return 2;
    for (size_t i = 0; i < WIN; ++i)
      for (size_t lane = 0; lane < LANES; ++lane)
        x[i * LANES + lane] = complex_input
            ? C(INPUT_GAIN * raw[2 * (i * LANES + lane)],
                INPUT_GAIN * raw[2 * (i * LANES + lane) + 1])
            : C(INPUT_GAIN * raw[i * LANES + lane], 0.0);
    std::printf("input_source=%s path=%s samples=%zu lanes=%zu\n",
                complex_input ? "complex_binary" : "real_binary",
                input_bin, raw.size(), LANES);
  } else {
    for (size_t b = 0; b < BLOCKS; ++b)
      for (size_t h = 0; h < 2; ++h)
        for (size_t i = 0; i < WIN; i++)
          x[i * LANES + 2 * b + h] =
              std::getenv("E52_REAL_INPUT")
                  ? C(INPUT_GAIN * .001 * sin(.01 * i + .07 * (4 * b + 2 * h)), 0.0)
                  : C(INPUT_GAIN * .001 * sin(.01 * i + .07 * (4 * b + 2 * h)),
                      INPUT_GAIN * .001 * cos(.013 * i + .11 * (4 * b + 2 * h + 1)));
  }
  auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x));
  auto F = fft3();
  auto P = proj(), B = mel(), D = dct();
  P = F[0];
  for (size_t s = 1; s < F.size(); ++s) P = mm(F[s], P);
  if (const char* norm_out = std::getenv("E52_WRITE_MEL_NORMS")) {
    std::ofstream no(norm_out);
    if (!no) return 2;
    for (size_t r = 0; r < 80; ++r) {
      double mx = 0.0;
      for (size_t lane = 0; lane < LANES; ++lane) {
        std::vector<C> fy(N);
        for (size_t k = 0; k < N; ++k)
          for (size_t j = 0; j < N; ++j)
            fy[k] += P[k][j] * x[j * LANES + lane];
        double er = 0.0, ei = 0.0;
        for (size_t k = 0; k < N; ++k) {
          er += B[r][k].real() * fy[k].real() * fy[k].real();
          ei += B[r][k].real() * fy[k].imag() * fy[k].imag();
        }
        mx = std::max(mx, std::max(er, ei));
      }
      no << std::setprecision(17) << (1.25 * mx) << '\n';
    }
    return 0;
  }
  std::vector<double> mel_norm(N, MEL_NORM);
  if (const char* norm_file = std::getenv("E52_MEL_NORM_FILE")) {
    std::ifstream nf(norm_file);
    for (size_t r = 0; r < 80; ++r) {
      if (!(nf >> mel_norm[r]) || !(mel_norm[r] > 0.0)) {
        std::fprintf(stderr, "failed to read 80 positive mel norms from %s\n",
                     norm_file);
        return 2;
      }
    }
  }
  // Normalize the mel energy to a public unit interval before the nonlinear
  // stage.  The inverse power scale is folded into the following DCT.
  for (size_t row = 0; row < B.size(); ++row)
    for (auto &z : B[row])
      z *= 0.125 / mel_norm[row];
  for (size_t row = 0; row < D.size(); ++row)
    for (size_t col = 0; col < D[row].size(); ++col)
      D[row][col] *= std::pow(mel_norm[col], POWER_ALPHA);
  auto enc_bsgs = [&](auto &M) {
    std::vector<Plaintext> v(N);
    for (size_t d = 0; d < N; d++) {
      bool active = false;
      size_t g = (d / MV_N1) * MV_N1;
      for (size_t r = 0; r < N; ++r) {
        size_t rr = (r + N - g) % N;
        if (std::abs(M[rr][(rr + d) % N]) > 1e-12) { active = true; break; }
      }
      if (!active) continue;
      std::vector<C> q(S);
      // The giant rotation is applied after the baby-step sum; pre-shift the
      // diagonal by its group offset so the output row remains r.
      for (size_t b = 0; b < BLOCKS; ++b)
        for (size_t h = 0; h < 2; ++h)
          for (size_t r = 0; r < N; r++) {
            size_t rr = (r + N - g) % N;
            q[r * LANES + 2 * b + h] = M[rr][(rr + d) % N];
          }
      v[d] = cc->MakeCKKSPackedPlaintext(q);
    }
    return v;
  };
  std::vector<std::vector<Plaintext>> fps;
  for (auto &M : F) fps.push_back(enc_bsgs(M));
  auto mp = enc_bsgs(B);
  // Mel produces 80 meaningful channels per lane.  Reserve the following
  // 80 slots in the same lane for the imaginary component so both real
  // components can share one Chebyshev evaluation.
  std::vector<C> mel_mask_v(S);
  for (size_t b = 0; b < BLOCKS; ++b)
    for (size_t h = 0; h < 2; ++h) {
      const size_t lane = 2 * b + h;
      for (size_t r = 0; r < 80; ++r) {
        mel_mask_v[r * LANES + lane] = C(1.0, 0.0);
      }
    }
  auto mel_mask = cc->MakeCKKSPackedPlaintext(mel_mask_v);
  auto dct_pts = enc_bsgs(D);
  auto fn = [](double z) {
    // POWER_OFFSET, when enabled, has already been added to the ciphertext.
    return std::pow(std::max(z, 0.0), POWER_ALPHA);
  };
  std::vector<C> ey(N);
  for (size_t r = 0; r < N; r++)
    for (size_t j = 0; j < N; j++)
      ey[r] += P[r][j] * x[j * LANES];
  std::vector<C> e_re(N), e_im(N), e_re2(N), e_im2(N), e_q(N);
  for (size_t r = 0; r < N; r++) {
    e_re[r] = C(2.0 * ey[r].real(), 0.0);
    e_im[r] = C(2.0 * ey[r].imag(), 0.0);
    e_re2[r] = C(4.0 * ey[r].real() * ey[r].real(), 0.0);
    e_im2[r] = C(4.0 * ey[r].imag() * ey[r].imag(), 0.0);
    e_q[r] = C(4.0 * ey[r].real() * ey[r].real(),
               4.0 * ey[r].imag() * ey[r].imag());
  }
  std::vector<C> e_m(N), e_mr(N), e_mi(N), e_cr(N), e_ci(N), e_c(N), e_out(N);
  for (size_t r = 0; r < N; r++) {
    for (size_t j = 0; j < N; j++) e_m[r] += B[r][j] * e_q[j];
    e_mr[r] = C(2.0 * e_m[r].real(), 0.0);
    e_mi[r] = C(2.0 * e_m[r].imag(), 0.0);
    e_c[r] = C(e_cr[r].real(), e_ci[r].real());
    for (size_t j = 0; j < N; j++) e_out[r] += D[r][j] * e_c[j];
  }
  std::vector<double> e_mr_real(N), e_mi_real(N);
  for (size_t r = 0; r < N; r++) {
    e_mr_real[r] = e_mr[r].real();
    e_mi_real[r] = e_mi[r].real();
  }
  auto exact_root = [](double z) {
    const double shifted = std::pow(std::max(z + POWER_OFFSET, 0.0),
                                    POWER_ALPHA);
    return ZERO_CENTER
        ? shifted - std::pow(std::max(POWER_OFFSET, 0.0), POWER_ALPHA)
        : shifted;
  };
  std::vector<double> e_cr_real(N), e_ci_real(N);
  for (size_t r = 0; r < N; ++r) {
    e_cr_real[r] = exact_root(e_mr_real[r]);
    e_ci_real[r] = exact_root(e_mi_real[r]);
  }
  if (std::getenv("E52_CHEB_CHECK")) {
    std::vector<double> probe = {0.0, 1.0e-4, 1.0e-3, 1.0e-2,
                                 1.0e-1, 5.0e-1, 1.0};
    for (uint32_t degree : {63u, 127u, 255u, 511u, 1023u, 2047u}) {
      auto approx = EvalChebyshevFunctionPtxt(fn, probe, CHEB_A, CHEB_B, degree);
      double max_err = 0.0;
      for (size_t i = 0; i < probe.size(); ++i) {
        const double err = approx[i] - exact_root(probe[i]);
        max_err = std::max(max_err, std::abs(err));
        if (degree == 63)
          printf("cheb_probe x=%.9g approx=%.9g exact=%.9g err=%.9g\n",
                 probe[i], approx[i], exact_root(probe[i]), err);
      }
      printf("cheb_degree=%u probe_maxerr=%.9g\n", degree, max_err);
    }
    for (double alpha : {0.25, 0.5, 0.75, 1.0}) {
      auto alpha_fn = [alpha](double z) {
        return std::pow(std::max(z, 0.0), alpha);
      };
      std::vector<double> dense(1001);
      for (size_t i = 0; i < dense.size(); ++i)
        dense[i] = static_cast<double>(i) / 1000.0;
      auto approx = EvalChebyshevFunctionPtxt(alpha_fn, dense, CHEB_A, CHEB_B, 63);
      double max_err = 0.0;
      for (size_t i = 0; i < dense.size(); ++i)
        max_err = std::max(max_err, std::abs(approx[i] - alpha_fn(dense[i])));
      printf("alpha=%.2f degree=63 dense_maxerr=%.9g\n", alpha, max_err);
    }
    printf("mel_input_stats re_max=%.9g im_max=%.9g re_root_max=%.9g im_root_max=%.9g\n",
           *std::max_element(e_mr_real.begin(), e_mr_real.end()),
           *std::max_element(e_mi_real.begin(), e_mi_real.end()),
           *std::max_element(e_cr_real.begin(), e_cr_real.end()),
           *std::max_element(e_ci_real.begin(), e_ci_real.end()));
  }
  for (size_t r = 0; r < N; r++) {
    e_cr[r] = C(e_cr_real[r], 0.0);
    e_ci[r] = C(e_ci_real[r], 0.0);
    e_c[r] = C(e_cr_real[r], e_ci_real[r]);
  }
  std::fill(e_out.begin(), e_out.end(), C(0.0, 0.0));
  for (size_t r = 0; r < N; r++)
    for (size_t j = 0; j < N; j++) e_out[r] += D[r][j] * e_c[j];
  auto conj = [&](CT z) { return cc->EvalAutomorphism(z, cj, *ck); };
  auto t0 = std::chrono::steady_clock::now();
  auto ts = t0;
  auto y = ct;
  const bool trace_proj = std::getenv("E52_PROJ_TRACE") != nullptr;
  std::vector<double> proj_stage_s;
  for (size_t i = 0; i < fps.size(); ++i) {
    auto stage_t0 = std::chrono::steady_clock::now();
    y = mv(cc, y, fps[i]);
    if (trace_proj)
      proj_stage_s.push_back(std::chrono::duration<double>(
          std::chrono::steady_clock::now() - stage_t0).count());
  }
  double t_proj = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
  if (trace_proj) {
    printf("projection_stage_s");
    for (size_t i = 0; i < proj_stage_s.size(); ++i)
      printf(" stage%zu=%.6f", i, proj_stage_s[i]);
    printf(" total=%.6f\n", t_proj);
  }
  if (!bench) probe_rmse(cc, y, kp.secretKey, ey, "projection");
  if (std::getenv("E52_GENERAL_CHECK")) {
    Plaintext py;
    cc->Decrypt(y, kp.secretKey, &py);
    py->SetLength(S);
    auto vy = py->GetCKKSPackedValue();
    double se = 0.0, me = 0.0;
    size_t count = 0;
    for (size_t lane = 0; lane < LANES; ++lane) {
      for (size_t r = 0; r < N; ++r) {
        C exp = 0.0;
        for (size_t j = 0; j < N; ++j)
          exp += P[r][j] * x[j * LANES + lane];
        const double e = std::abs(vy[r * LANES + lane] - exp);
        se += e * e;
        me = std::max(me, e);
        ++count;
      }
    }
    printf("general_projection_all_lanes rmse=%.6e maxerr=%.6e lanes=%zu\n",
           std::sqrt(se / count), me, LANES);
  }
  auto yc = conj(y);
  auto re = cc->EvalAdd(y, yc);
  auto im = cc->EvalSub(y, yc);
  mul_i_inplace(im);
  for (auto &c : im->GetElements())
    c = c.Times(-1);
  if (!bench) probe_rmse(cc, re, kp.secretKey, e_re, "split_re");
  if (!bench) probe_rmse(cc, im, kp.secretKey, e_im, "split_im");
  auto re2 = cc->EvalSquare(re);
  auto im2 = cc->EvalSquare(im);
  cc->RescaleInPlace(re2);
  cc->RescaleInPlace(im2);
  if (!bench) probe_rmse(cc, re2, kp.secretKey, e_re2, "square_re");
  if (!bench) probe_rmse(cc, im2, kp.secretKey, e_im2, "square_im");
  CT im2i = std::make_shared<CiphertextImpl<DCRTPoly>>(*im2);
  mul_i_inplace(im2i);
  auto q = cc->EvalAdd(re2, im2i);
  if (!bench) probe_rmse(cc, q, kp.secretKey, e_q, "recombine_q");
  auto m = mv(cc, q, mp);
  double t_mel = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count() - t_proj;
  if (!bench) probe_rmse(cc, m, kp.secretKey, e_m, "mel");
  // Fold the output mask into the mel stage.  The mel matrix has only 80
  // nonzero rows, but an explicit mask is still needed to suppress CKKS
  // noise in unused slots before the later lane-local rotation.  Applying it
  // once to m is equivalent to masking both split branches.
  m = cc->EvalMult(m, mel_mask);
  auto mc = conj(m);
  auto mr = cc->EvalAdd(m, mc);
  auto mi = cc->EvalSub(m, mc);
  mul_i_inplace(mi);
  for (auto &c : mi->GetElements())
    c = c.Times(-1);
  if (!bench) probe_rmse(cc, mr, kp.secretKey, e_mr, "split2_re");
  if (!bench) probe_rmse(cc, mi, kp.secretKey, e_mi, "split2_im");

  // Pack the two already-masked real ciphertexts into disjoint slot regions.
  // The single rescale keeps the two branches aligned.
  auto mr_packed = mr;
  auto mi_packed = mi;
  // OpenFHE's positive rotation maps output slot r to input slot r+offset;
  // use -80 to move source slot r into destination slot r+80.
  auto mi_shifted = cc->EvalRotate(mi_packed, -80 * static_cast<int>(ROT_STRIDE));
  auto packed_mel = cc->EvalAdd(mr_packed, mi_shifted);
  if (POWER_OFFSET != 0.0) {
    // EvalChebyshev is applied to every slot, including masked/unused ones.
    // Shift all slots so CKKS noise around zero cannot be evaluated outside
    // the public approximation interval.
    std::vector<C> offset_v(S, C(POWER_OFFSET, 0.0));
    auto offset_pt = cc->MakeCKKSPackedPlaintext(offset_v);
    packed_mel = cc->EvalAdd(packed_mel, offset_pt);
  }
  cc->RescaleInPlace(packed_mel);
  if (!bench) {
    std::vector<C> packed_expected(N);
    for (size_t r = 0; r < 80; ++r) {
      packed_expected[r] = e_mr[r];
      packed_expected[80 + r] = e_mi[r];
    }
    probe_rmse(cc, packed_mel, kp.secretKey, packed_expected, "packed_mel");
  }
  auto packed_cheb =
      cc->EvalChebyshevFunction(fn, packed_mel, CHEB_A, CHEB_B, CHEB_DEGREE);
  if (ZERO_CENTER) {
    // The ciphertext was shifted by POWER_OFFSET before evaluation.  Remove
    // the exact shift image, not the polynomial's extrapolated value at zero;
    // for x^(1/4), the latter has a large endpoint approximation error.
    const double offset_image =
        std::pow(std::max(POWER_OFFSET, 0.0), POWER_ALPHA);
    std::vector<C> zero_v(S, C(offset_image, 0.0));
    auto zero_pt = cc->MakeCKKSPackedPlaintext(zero_v);
    packed_cheb = cc->EvalSub(packed_cheb, zero_pt);
  }

  // Keep the real output in the original first-80 slots.  Extract the
  // imaginary output from the shifted region and rotate it back.  With the
  // OpenFHE rotation convention, +80 maps destination slot r to source slot
  // r+80.  The next DCT linear transform performs the required rescale.
  auto cr = packed_cheb;
  auto ci = cc->EvalRotate(packed_cheb, 80 * static_cast<int>(ROT_STRIDE));
  if (!bench) {
    std::vector<C> e_cr_valid(e_cr.begin(), e_cr.begin() + 80);
    std::vector<C> e_ci_valid(e_ci.begin(), e_ci.begin() + 80);
    probe_rmse(cc, cr, kp.secretKey, e_cr_valid, "cheb_re");
    probe_rmse(cc, ci, kp.secretKey, e_ci_valid, "cheb_im");
  }
  CT cii = std::make_shared<CiphertextImpl<DCRTPoly>>(*ci);
  mul_i_inplace(cii);
  auto c = cc->EvalAdd(cr, cii);
  double t_cheb = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count() - t_proj - t_mel;
  if (!bench) {
    std::vector<C> e_c_valid(e_c.begin(), e_c.begin() + 80);
    probe_rmse(cc, c, kp.secretKey, e_c_valid, "recombine_c");
  }
  auto mr_val = decrypt_values(cc, mr, kp.secretKey);
  auto mi_val = decrypt_values(cc, mi, kp.secretKey);
  auto c_val_pre = decrypt_values(cc, c, kp.secretKey);
  print_stats("mr", mr_val);
  print_stats("mi", mi_val);
  print_stats("c_re", c_val_pre);
  std::vector<C> c_im(N);
  for (size_t i = 0; i < N; ++i) c_im[i] = C(c_val_pre[i].imag(), 0.0);
  print_stats("c_im", c_im);
  for (size_t k = 0; k < 13; ++k) {
    C sum = 0; double abs_sum = 0;
    for (size_t n = 0; n < 80; ++n) {
      C term = D[k][n] * c_val_pre[n];
      sum += term; abs_sum += std::abs(term);
    }
    printf("dct_row k=%zu value_abs=%.6e abs_terms=%.6e ratio=%.6e\n",
           k, std::abs(sum), abs_sum, abs_sum / std::max(std::abs(sum), 1e-30));
  }
  auto out = mv(cc, c, dct_pts);
  double t_dct = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count() - t_proj - t_mel - t_cheb;
  if (!bench) probe_rmse(cc, out, kp.secretKey, e_out, "dct_out");
  Plaintext c_dec;
  cc->Decrypt(c, kp.secretKey, &c_dec);
  c_dec->SetLength(S);
  auto c_val = c_dec->GetCKKSPackedValue();
  std::vector<C> dct_from_decoded_c(N);
  for (size_t r = 0; r < N; ++r)
    for (size_t j = 0; j < N; ++j)
      dct_from_decoded_c[r] += D[r][j] * c_val[j * LANES];
  if (!bench) probe_rmse(cc, out, kp.secretKey, dct_from_decoded_c, "dct_vs_decoded");
  auto sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
          .count();
  Plaintext dec;
  cc->Decrypt(out, kp.secretKey, &dec);
  dec->SetLength(S);
  auto v = dec->GetCKKSPackedValue();
  if (const char* output_bin = std::getenv("E52_OUTPUT_BIN")) {
    // Layout: lane-major [LANES][13][real, imag], float32.  Keeping both
    // components preserves the complex-packing result for downstream code;
    // real-input callers will observe an approximately zero imaginary part.
    std::ofstream out_file(output_bin, std::ios::binary);
    if (!out_file) {
      std::fprintf(stderr, "failed to open E52_OUTPUT_BIN=%s\n", output_bin);
      return 2;
    }
    for (size_t lane = 0; lane < LANES; ++lane)
      for (size_t k = 0; k < 13; ++k) {
        const C z = v[k * LANES + lane];
        const float pair[2] = {static_cast<float>(z.real()),
                               static_cast<float>(z.imag())};
        out_file.write(reinterpret_cast<const char*>(pair), sizeof(pair));
      }
    if (!out_file) {
      std::fprintf(stderr, "failed to write E52_OUTPUT_BIN=%s\n", output_bin);
      return 2;
    }
  }
  std::vector<double> qr(N), qi(N);
  for (size_t r = 0; r < N; r++) {
    qr[r] = ey[r].real() * ey[r].real();
    qi[r] = ey[r].imag() * ey[r].imag();
  }
  std::vector<double> er(N), ei(N);
  for (size_t r = 0; r < N; r++)
    for (size_t j = 0; j < N; j++) {
      er[r] += 8.0 * (B[r][j] * qr[j]).real();
      ei[r] += 8.0 * (B[r][j] * qi[j]).real();
    }
  std::vector<double> zr(N), zi(N);
  for (size_t r = 0; r < N; ++r) {
    zr[r] = exact_root(er[r]);
    zi[r] = exact_root(ei[r]);
  }
  std::vector<C> ec(N), eo(N);
  for (size_t r = 0; r < N; r++)
    ec[r] = C(zr[r], zi[r]);
  for (size_t r = 0; r < N; r++) {
    for (size_t j = 0; j < N; j++)
      eo[r] += D[r][j] * ec[j];
  }
  double se = 0, me = 0;
  for (size_t i = 0; i < 13; i++) {
    double e = std::abs(v[i * LANES] - eo[i]);
    se += e * e;
    me = std::max(me, e);
  }
  printf("runtime_s=%.3f level=%zu rmse=%.6e maxerr=%.6e\n", sec,
         out->GetLevel(), std::sqrt(se / 13), me);
  const double real_frames = 2.0 * static_cast<double>(LANES);
  printf("throughput_fps=%.6f\n", real_frames / sec);
  printf("bench_stage_s projection=%.6f mel=%.6f cheb=%.6f dct=%.6f\n",
         t_proj, t_mel, t_cheb, t_dct);
  printf("bench_stage_fps projection=%.6f mel=%.6f cheb=%.6f dct=%.6f\n",
         real_frames / t_proj, real_frames / t_mel, real_frames / t_cheb,
         real_frames / t_dct);
  for (int i = 0; i < 4; i++)
    printf("slot%d=(%.6g,%.6g) exp=(%.6g,%.6g)\n", i, v[i * LANES].real(),
           v[i * LANES].imag(),
           eo[i].real(), eo[i].imag());
  if (std::getenv("E52_DUMP_OUTPUT")) {
    for (size_t i = 0; i < 13; ++i) {
      const C got = v[i * LANES];
      const double err = std::abs(got - eo[i]);
      printf("decrypted[%zu]=(%.17g,%.17g) expected=(%.17g,%.17g) abs_err=%.9e\n",
             i, got.real(), got.imag(), eo[i].real(), eo[i].imag(), err);
    }
  }
  return 0;
}
