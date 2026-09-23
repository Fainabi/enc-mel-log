#!/usr/bin/env python3
"""Plaintext reference for the current complex-packed C++ frontend."""
import argparse
import json
from pathlib import Path

import numpy as np
import soundfile as sf


N, WIN, HOP, LANES = 512, 400, 160, 64


def frame_feature(a, b, norms, alpha, offset):
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(WIN) / WIN)
    x = np.zeros(N, dtype=np.complex128)
    x[:WIN] = w * (a + 1j * b)
    z = np.fft.fft(x)
    p = np.zeros(N, dtype=np.complex128)
    p[0] = z[0]
    for k in range(1, N // 2):
        p[2 * k - 1] = 0.5 * (z[k] + z[N - k])
        p[2 * k] = -0.5j * z[k] + 0.5j * z[N - k]
    p[N - 1] = z[N // 2]
    q = 4 * p.real * p.real + 1j * 4 * p.imag * p.imag
    mel = np.zeros(80, dtype=np.complex128)
    for r in range(80):
        lo, mid, hi = 1 + 3 * r, 4 + 3 * r, 7 + 3 * r
        for k in range(lo, hi):
            weight = (k - lo) / 3 if k <= mid else (hi - k) / 3
            mel[r] += weight * q[k] * (0.125 / norms[r])
    shift = offset ** alpha
    c = np.power(np.maximum(2 * mel.real + offset, 0), alpha) - shift
    ci = np.power(np.maximum(2 * mel.imag + offset, 0), alpha) - shift
    dct = np.zeros(13, dtype=np.complex128)
    for k in range(13):
        for r in range(80):
            dct[k] += (np.cos(np.pi * k * (r + 0.5) / 80) *
                       (norms[r] ** alpha)) * (c[r] + 1j * ci[r])
    return dct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segments", required=True)
    ap.add_argument("--norms", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--alpha", type=float, default=0.25)
    ap.add_argument("--offset", type=float, default=1e-4)
    args = ap.parse_args()
    norms = np.loadtxt(args.norms)[:80]
    segs = json.loads(Path(args.segments).read_text())
    out = np.empty((len(segs), 100, 13, 2), dtype=np.float32)
    for si, s in enumerate(segs):
        x, sr = sf.read(s["path"], dtype="float32")
        assert sr == 16000
        base = s["segment"] * 100 * HOP
        frames = np.stack([x[base + j * HOP:base + j * HOP + WIN]
                           for j in range(100)])
        for j in range(0, 100, 2):
            y = frame_feature(frames[j], frames[j + 1], norms,
                              args.alpha, args.offset)
            out[si, j, :, 0] = y.real
            out[si, j + 1, :, 0] = y.imag
    np.save(args.output, out)
    print(f"output={args.output} shape={out.shape}")


if __name__ == "__main__":
    main()
