# ppaudio-he-frontend

Reference C++ implementation of a CKKS encrypted log-mel speech front end for the
ICASSP 2027 submission, *Multiplicative Depth Budgets for Encrypted Log-Mel Front-Ends
Are Set by Dynamic Range, Not Approximability*.

The main executable, `cpp_complex_e2e`, evaluates the following chain on encrypted
input:

```text
windowed STFT projection -> modulus square -> mel filterbank
-> power compressor u^alpha (Chebyshev approximation) -> DCT
```

The prototype uses complex CKKS packing, a three-stage radix-8 DFT, conjugation-based
real/imaginary extraction, sparse plaintext linear transforms, and a 13-coefficient
DCT output. It is intended for reproducible experiments and implementation diagnosis,
not as production software.

## Scope

This repository contains:

- the OpenFHE C++ encrypted front-end prototypes;
- unit tests for complex CKKS operations, inverse-two placement, linear transforms, and DCT;
- plaintext and decrypted-HE benchmark scripts used for the downstream screen.

It does not contain framing, model training, bootstrapping, client/server serialization,
or production error handling. The default C++ benchmark uses deterministic synthetic
input. Real LibriSpeech input can be supplied through the binary-input interface below.

The current C++ timing path is a degree-63 prototype with alpha = 3/4. The paper's
accuracy-selected degree-2047 configuration has not been run end to end in this
repository; its cost is estimated from separately measured operator prices.

## Requirements

- OpenFHE with complex CKKS support (`CKKSDataType` and `SetCKKSDataType(COMPLEX)`);
- CMake >= 3.16;
- a C++17 compiler;
- OpenMP-enabled OpenFHE;
- `NATIVE_SIZE=64` and `MATHBACKEND=4` in the OpenFHE build;
- at least 8 GB of available memory for the full `N = 2^16` prototype.

### OpenFHE version

The paper records OpenFHE 1.4.0. The local build used during development was based on
OpenFHE 1.2.4, but the source tree used by this repository requires the newer complex
CKKS APIs. Confirm the version before compiling. A clean installation can be made as
follows:

```bash
git clone --branch v1.4.0 https://github.com/openfheorg/openfhe-development.git
cd openfhe-development
cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/.local \
      -DWITH_OPENMP=ON \
      -DBUILD_SHARED=ON \
      -DBUILD_UNITTESTS=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_BENCHMARKS=OFF \
      -DNATIVE_SIZE=64
cmake --build build -j"$(nproc)"
cmake --install build
```

Then configure this project with either:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$HOME/.local
```

or, if necessary:

```bash
cmake -S . -B build -DOpenFHE_DIR=$HOME/.local/lib/OpenFHE
```

The OpenFHE include directory reported by CMake must exist. A build tree that has not
been installed may export paths pointing at a nonexistent `/usr/local` prefix.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

To build only the main program:

```bash
cmake --build build --target cpp_complex_e2e -j"$(nproc)"
```

The build is out of source and places executables in `build/`. OpenFHE compiler flags
include `-Wall -Werror`; a compiler upgrade may therefore expose new warnings.

## Executables

| Target | Purpose |
| --- | --- |
| `cpp_complex_e2e` | Full complex-packed chain: radix-8 projection, conjugate split, square, mel, power compressor, and DCT. Reports per-stage RMSE, max error, levels, and timing. |
| `cpp_complex_ops` | Smoke test for multiplication by `i`, inverse-two handling, conjugation, real/imaginary split, square, and recombination. |
| `cpp_complex_smoke` | Full-slot complex ciphertext/plaintext multiplication test. It returns failure when the maximum error exceeds `1e-5`. |
| `cpp_dct_small` | Signed-diagonal DCT test on fresh and deliberately deep ciphertexts. |
| `cpp_inv2_placement` | Compares inverse-two placement strategies and checks DCRT and ciphertext round trips. |
| `cpp_linear_smoke` | Comparison against OpenFHE's built-in linear-transform implementation. |
| `cpp_mv_unit` | Unit test for the hand-written BSGS matrix-vector routine and DCT diagonal layout. |
| `cpp_real_dct` | Historical REAL-data-type comparison path. Its recorded high-order DCT error is large and it is not the current method. |

All programs write diagnostics to stdout. Key generation and plaintext diagonal
encoding can dominate process time; the reported encrypted-evaluation time does not
necessarily include those costs.

## Main-program parameters

The main source is `src/cpp_complex_e2e.cpp`. Important defaults are:

| Parameter | Value |
| --- | --- |
| CKKS ring dimension | `2^16` (32768 complex slots) |
| Multiplicative-depth setting | `20` |
| Scaling modulus | `59` bits |
| Security level | `HEStd_128_classic` |
| CKKS data type | `COMPLEX` |
| FFT length | `512` |
| Window and hop | `400` samples and `160` samples |
| Default frame layout | `64` frames, 32 blocks, 1024 slots per block |
| Mel channels | `80` |
| DCT coefficients retained | `13` |
| Compressor | `u^(3/4)` |
| Chebyshev degree | `63` |
| BSGS baby-step count | `16` |

Useful environment variables:

- `E52_BENCH=1`: skip the detailed stage probes;
- `E52_INPUT_BIN=...`: read little-endian float32 input;
- `E52_MEL_NORM_FILE=...`: load 80 public per-channel calibration scales;
- `E52_POWER_ALPHA=...`: override the compressor exponent;
- `E52_POWER_OFFSET=...`: override the compressor offset;
- `E52_CHEB_DEGREE=...`: override the Chebyshev degree;
- `E52_ZERO_CENTER=1`: enable the real-speech zero-centering path;
- `E52_WRITE_MEL_NORMS=...`: write calibration scales for a supplied input.

The `level=` field printed by OpenFHE is an internal level index. The paper reports
remaining levels below the configured ceiling; do not interpret the two quantities as
the same number.

## Python benchmark scripts

The scripts in `scripts/` support the plaintext/decrypted-HE downstream screen:

- `generate_uniform_dct_dataset.py` generates a plaintext feature dataset matching the
  uniform filterbank used by the C++ prototype;
- `benchmark_plain_models.py` trains and evaluates temporal CNN models and their
  hard/soft ensembles;
- `evaluate_he_test400.py` compares decrypted HE features with the plaintext reference
  and evaluates the saved downstream models.

These scripts require Python 3, NumPy, PyTorch, and (for dataset generation) SoundFile.
They are evaluation utilities, not part of the C++ build.

## Reproduction

Build and run the main timing prototype:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cpp_complex_e2e -j"$(nproc)"
E52_BENCH=1 ./build/cpp_complex_e2e
```

The historical reference run on an Intel Xeon used the same `N = 2^16` configuration:

| Quantity | Reference |
| --- | ---: |
| Encrypted-evaluation runtime | `43.800 s` |
| RMSE on the first 13 DCT coefficients | `8.436337e-03` |
| Maximum error | `8.973544e-03` |
| Projection stage | `22.710301 s` |
| Mel stage | `9.747905 s` |
| Compressor stage | `5.058857 s` |
| DCT stage | `6.220999 s` |

Runtime is machine- and thread-dependent. The reference runtime excludes context
construction, key generation, encryption, decryption, and some CPU-side preprocessing.

### Real-speech smoke benchmark

For a supplied 128-frame LibriSpeech binary input, first generate public mel scales:

```bash
E52_INPUT_BIN=/path/libri_128frames_complex_f32.bin \
E52_WRITE_MEL_NORMS=/path/mel_norms_cpp.txt \
./build/cpp_complex_e2e
```

Then run the calibrated chain:

```bash
E52_BENCH=1 E52_ZERO_CENTER=1 E52_CHEB_DEGREE=511 \
E52_POWER_ALPHA=0.75 E52_POWER_OFFSET=1e-6 \
E52_MEL_NORM_FILE=/path/mel_norms_cpp.txt \
E52_INPUT_BIN=/path/libri_128frames_complex_f32.bin \
./build/cpp_complex_e2e
```

The recorded Xeon run was `42.993 s` for 128 frames, with endpoint RMSE
`7.66e-03` and maximum error `1.81e-02`. This is a front-end consistency check,
not a downstream accuracy result.

## Known limitations

1. The C++ prototype uses 80 uniformly spaced three-bin triangular filters, not the
   HTK mel filterbank used by the paper's speech experiments. Its timing therefore
   represents sparse linear-transform cost, not final feature quality.
2. Reported RMSE is computed on the first 13 DCT coefficients, not all ciphertext slots.
3. The default input is synthetic. Real speech requires calibration scales from the same
   public data distribution.
4. The radix-8 factorization is faster but less accurate than the dense projection path;
   the two historical measurements also use different frame layouts.
5. The degree-2047 accuracy target has not been run end to end here. Do not combine its
   estimated cost with the degree-63 timing as if they were one measurement.
6. `cpp_real_dct` is a failed historical comparison with high-order coefficient error;
   it is retained for diagnosis only.
7. The current repository has not been fully revalidated against every OpenFHE release.
   Confirm the installed headers, libraries, and complex CKKS APIs before reporting new
   measurements.

## License

The licensing status is currently TBD.
