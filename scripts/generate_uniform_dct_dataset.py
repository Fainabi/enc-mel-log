#!/usr/bin/env python3
"""Generate plaintext features matching cpp_complex_e2e's uniform filterbank."""
import argparse
from pathlib import Path

import numpy as np
import soundfile as sf


SR, WIN, HOP, FFT, FRAMES = 16000, 400, 160, 512, 100


def make_features(x, norms, alpha, offset):
    nseg = ((len(x) - WIN) // HOP + 1) // FRAMES
    if nseg <= 0:
        return np.empty((0, FRAMES, 13), dtype=np.float32)
    starts = np.arange(nseg * FRAMES)[:, None] * HOP
    frames = x[starts + np.arange(WIN)[None, :]].reshape(nseg * FRAMES, WIN)
    frames = frames * (0.5 - 0.5 * np.cos(2 * np.pi * np.arange(WIN) / WIN))
    z = np.fft.rfft(frames, n=FFT, axis=1)
    p = np.empty((len(frames), FFT), dtype=np.float64)
    p[:, 0] = z[:, 0].real
    p[:, -1] = z[:, -1].real
    for k in range(1, FFT // 2):
        p[:, 2 * k - 1] = z[:, k].real
        p[:, 2 * k] = z[:, k].imag
    weights = np.zeros((80, FFT), dtype=np.float64)
    for r in range(80):
        lo, mid, hi = 1 + 3 * r, 4 + 3 * r, 7 + 3 * r
        for k in range(lo, hi):
            weights[r, k] = (k - lo) / 3 if k <= mid else (hi - k) / 3
    mel = (p * p) @ weights.T
    c = np.maximum(mel / norms[None, :] + offset, 0.0) ** alpha
    c -= offset ** alpha
    dct = c * norms[None, :] ** alpha
    basis = np.cos(np.pi * np.arange(13)[:, None] *
                   (np.arange(80)[None, :] + 0.5) / 80.0)
    out = dct @ basis.T
    return out.reshape(nseg, FRAMES, 13).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--feat", required=True, help="feat.npz containing u")
    ap.add_argument("--norms", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--alpha", type=float, default=0.25)
    ap.add_argument("--offset", type=float, default=1e-4)
    args = ap.parse_args()

    u = np.load(args.feat)["u"]
    norms = np.loadtxt(args.norms)[:80].astype(np.float64)
    files = sorted(Path(args.data).glob("*/*/*.flac"))
    out = np.empty((len(u), FRAMES, 13), dtype=np.float32)
    pos = 0
    for uid, path in enumerate(files):
        rows = np.flatnonzero(u == uid)
        if not len(rows):
            continue
        x, sr = sf.read(path, dtype="float32")
        if sr != SR:
            raise ValueError(f"unexpected sample rate for {path}: {sr}")
        f = make_features(x, norms, args.alpha, args.offset)
        if len(f) != len(rows):
            raise ValueError(f"segment mismatch uid={uid}: {len(f)} != {len(rows)}")
        out[rows] = f
        pos += len(rows)
        if uid % 250 == 0:
            print(f"utterances={uid} rows={pos}", flush=True)
    np.save(args.output, out)
    print(f"output={args.output} shape={out.shape}")


if __name__ == "__main__":
    main()
