#include "openfhe.h"
#include "math/chebyshev.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace lbcrypto;
using C = std::complex<double>;
using CT = Ciphertext<DCRTPoly>;
static constexpr size_t N = 512, WIN = 400, BLOCK = 1024, FRAMES = 64;
static constexpr size_t BLOCKS = FRAMES / 2;
static constexpr size_t MV_N1 = 16;
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
  for (size_t g = 0; g < 3; ++g) {
    std::vector<std::vector<C>> M(N, std::vector<C>(N));
    for (size_t i = 0; i < N; ++i) M[i][i] = 1.0;
    for (size_t s = 0; s < 3; ++s) M = mm(stages[g * 3 + s], M);
    if (g == 0)
      for (size_t j = 0; j < WIN; ++j) {
        double w = .5 - .5 * cos(2 * M_PI * j / WIN);
        for (size_t i = 0; i < N; ++i) M[i][j] *= w;
      }
    out.push_back(std::move(M));
  }
  std::vector<std::vector<C>> br(N, std::vector<C>(N));
  for (size_t k = 0; k < N; ++k) {
    size_t q = 0;
    for (int b = 0; b < 9; ++b) q = (q << 1) | ((k >> b) & 1);
    br[k][q] = 1.0;
  }
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
  for (int r = 0; r < 80; r++) {
    int lo = 1 + r * 3, mid = lo + 3, hi = mid + 3;
    for (int k = lo; k < hi && k < 257; k++)
      M[r][k] = (k <= mid) ? double(k - lo) / (mid - lo)
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
static void probe(const CryptoContext<DCRTPoly> &cc, const CT &ct,
                  const PrivateKey<DCRTPoly> &sk, const char *n) {
  Plaintext p;
  cc->Decrypt(ct, sk, &p);
  p->SetLength(4);
  auto v = p->GetCKKSPackedValue();
  printf("%s level=%zu s0=(%.5g,%.5g) s1=(%.5g,%.5g)\n", n, ct->GetLevel(),
         v[0].real(), v[0].imag(), v[1].real(), v[1].imag());
}
static void probe_rmse(const CryptoContext<DCRTPoly> &cc, const CT &ct,
                       const PrivateKey<DCRTPoly> &sk,
                       const std::vector<C> &expected, const char *name) {
  Plaintext p;
  cc->Decrypt(ct, sk, &p);
  p->SetLength(expected.size());
  auto v = p->GetCKKSPackedValue();
  double se = 0.0, me = 0.0;
  for (size_t i = 0; i < expected.size(); ++i) {
    double e = std::abs(v[i] - expected[i]);
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
  p->SetLength(N);
  return p->GetCKKSPackedValue();
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
    if (used_i[i]) baby[i] = cc->EvalFastRotation(x, i, 2 * cc->GetRingDimension(), dig);
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
      inner = cc->EvalRotate(inner, (int)(k * n1));
    if (first) {
      acc = inner;
      first = false;
    } else
      cc->EvalAddInPlace(acc, inner);
  }
  cc->RescaleInPlace(acc);
  return acc;
}
static CT mv_signed(const CryptoContext<DCRTPoly> &cc, CT x,
                    const std::vector<Plaintext> &pts,
                    const std::vector<int> &ds, int offset) {
  CT acc;
  bool first = true;
  for (int d : ds) {
    CT rot = d ? cc->EvalRotate(x, d) : x;
    CT term = cc->EvalMult(rot, pts[d + offset]);
    if (first) {
      acc = term;
      first = false;
    } else {
      cc->EvalAddInPlace(acc, term);
    }
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
  p.SetRingDim(1 << 16);
  p.SetSecurityLevel(HEStd_128_classic);
  p.SetScalingTechnique(FIXEDMANUAL);
  p.SetCKKSDataType(COMPLEX);
  auto cc = GenCryptoContext(p);
  for (auto f : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE})
    cc->Enable(f);
  auto kp = cc->KeyGen();
  cc->EvalMultKeyGen(kp.secretKey);
  std::vector<int32_t> rk;
  for (int i = 1; i < (int)N; i++)
    rk.push_back(i);
  for (int i = 1; i <= 12; ++i)
    rk.push_back(-i);
  cc->EvalRotateKeyGen(kp.secretKey, rk);
  uint32_t cj = 2 * cc->GetRingDimension() - 1;
  auto ck = cc->EvalAutomorphismKeyGen(kp.secretKey, {cj});
  size_t S = cc->GetRingDimension() / 2;
  const bool bench = std::getenv("E52_BENCH") != nullptr;
  std::vector<C> x(S);
  for (size_t b = 0; b < BLOCKS; ++b)
    for (size_t h = 0; h < 2; ++h)
      for (size_t i = 0; i < WIN; i++)
        x[b * BLOCK + h * N + i] =
            C(.001 * sin(.01 * i + .07 * b),
              .001 * cos(.013 * i + .11 * b));
  auto ct = cc->Encrypt(kp.publicKey, cc->MakeCKKSPackedPlaintext(x));
  auto F = fft3();
  auto P = proj(), B = mel(), D = dct();
  P = F[0];
  for (size_t s = 1; s < F.size(); ++s) P = mm(F[s], P);
  for (auto &row : B)
    for (auto &z : row)
      z *= 0.125;
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
            q[b * BLOCK + h * N + r] = M[rr][(rr + d) % N];
          }
      v[d] = cc->MakeCKKSPackedPlaintext(q);
    }
    return v;
  };
  std::vector<std::vector<Plaintext>> fps;
  for (auto &M : F) fps.push_back(enc_bsgs(M));
  auto mp = enc_bsgs(B);
  std::vector<int> dct_ds;
  std::vector<Plaintext> dct_direct(N);
  for (int d = -12; d <= 79; ++d) {
    int dm = (d + static_cast<int>(N)) % static_cast<int>(N);
    if (std::find(dct_ds.begin(), dct_ds.end(), dm) == dct_ds.end())
      dct_ds.push_back(dm);
    std::vector<C> q(S);
    for (size_t b = 0; b < BLOCKS; ++b)
      for (size_t h = 0; h < 2; ++h)
        for (size_t r = 0; r < N; ++r) {
          int n = static_cast<int>(r) + d;
          if (n >= 0 && n < static_cast<int>(N))
            q[b * BLOCK + h * N + r] = D[r][n];
        }
    dct_direct[dm] = cc->MakeCKKSPackedPlaintext(q);
  }
  auto fn = [](double z) { return std::pow(std::max(z, 0.0), .25); };
  std::vector<C> ey(N);
  for (size_t r = 0; r < N; r++)
    for (size_t j = 0; j < N; j++)
      ey[r] += P[r][j] * x[j];
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
  auto ptxt_fn = [&](const std::vector<double>& in) {
    return EvalChebyshevFunctionPtxt(fn, in, 0, 1, 63);
  };
  auto e_cr_real = ptxt_fn(e_mr_real);
  auto e_ci_real = ptxt_fn(e_mi_real);
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
  for (auto &fp : fps) y = mv(cc, y, fp);
  double t_proj = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
  if (!bench) probe_rmse(cc, y, kp.secretKey, ey, "projection");
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
  auto mc = conj(m);
  auto mr = cc->EvalAdd(m, mc);
  auto mi = cc->EvalSub(m, mc);
  mul_i_inplace(mi);
  for (auto &c : mi->GetElements())
    c = c.Times(-1);
  if (!bench) probe_rmse(cc, mr, kp.secretKey, e_mr, "split2_re");
  if (!bench) probe_rmse(cc, mi, kp.secretKey, e_mi, "split2_im");
  auto cr = cc->EvalChebyshevFunction(fn, mr, 0, 1, 63);
  auto ci = cc->EvalChebyshevFunction(fn, mi, 0, 1, 63);
  if (!bench) probe_rmse(cc, cr, kp.secretKey, e_cr, "cheb_re");
  if (!bench) probe_rmse(cc, ci, kp.secretKey, e_ci, "cheb_im");
  CT cii = std::make_shared<CiphertextImpl<DCRTPoly>>(*ci);
  mul_i_inplace(cii);
  auto c = cc->EvalAdd(cr, cii);
  double t_cheb = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count() - t_proj - t_mel;
  if (!bench) probe_rmse(cc, c, kp.secretKey, e_c, "recombine_c");
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
  auto out = mv_signed(cc, c, dct_direct, dct_ds, 0);
  double t_dct = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count() - t_proj - t_mel - t_cheb;
  if (!bench) probe_rmse(cc, out, kp.secretKey, e_out, "dct_out");
  Plaintext c_dec;
  cc->Decrypt(c, kp.secretKey, &c_dec);
  c_dec->SetLength(N);
  auto c_val = c_dec->GetCKKSPackedValue();
  std::vector<C> dct_from_decoded_c(N);
  for (size_t r = 0; r < N; ++r)
    for (size_t j = 0; j < N; ++j) dct_from_decoded_c[r] += D[r][j] * c_val[j];
  if (!bench) probe_rmse(cc, out, kp.secretKey, dct_from_decoded_c, "dct_vs_decoded");
  auto sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
          .count();
  Plaintext dec;
  cc->Decrypt(out, kp.secretKey, &dec);
  dec->SetLength(13);
  auto v = dec->GetCKKSPackedValue();
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
  auto zr = ptxt_fn(er);
  auto zi = ptxt_fn(ei);
  std::vector<C> ec(N), eo(N);
  for (size_t r = 0; r < N; r++)
    ec[r] = C(zr[r], zi[r]);
  for (size_t r = 0; r < N; r++) {
    for (size_t j = 0; j < N; j++)
      eo[r] += D[r][j] * ec[j];
  }
  double se = 0, me = 0;
  for (size_t i = 0; i < 13; i++) {
    double e = std::abs(v[i] - eo[i]);
    se += e * e;
    me = std::max(me, e);
  }
  printf("runtime_s=%.3f level=%zu rmse=%.6e maxerr=%.6e\n", sec,
         out->GetLevel(), std::sqrt(se / 13), me);
  printf("bench_stage_s projection=%.6f mel=%.6f cheb=%.6f dct=%.6f\n",
         t_proj, t_mel, t_cheb, t_dct);
  for (int i = 0; i < 4; i++)
    printf("slot%d=(%.6g,%.6g) exp=(%.6g,%.6g)\n", i, v[i].real(), v[i].imag(),
           eo[i].real(), eo[i].imag());
  return 0;
}
