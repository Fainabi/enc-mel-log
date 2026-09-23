#!/usr/bin/env python3
"""Run independent CKKS frontend jobs concurrently and collect features.

Each input block must contain one E52 frame batch (400*64 float32 samples, or
twice that for complex input).  The C++ worker writes 64 lanes x 13 DCT bins x
two float components in lane-major order.
"""
import argparse
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import subprocess

import numpy as np


def run_one(args):
    binary, inp, out, env = args
    e = os.environ.copy()
    e.update(env)
    e.update({"E52_BENCH": "1", "E52_INPUT_BIN": str(inp),
              "E52_OUTPUT_BIN": str(out)})
    subprocess.run([binary], env=e, stdout=subprocess.DEVNULL,
                   stderr=subprocess.PIPE, check=True, text=True)
    return str(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--inputs", required=True,
                    help="directory containing sorted *.bin frame blocks")
    ap.add_argument("--output", required=True,
                    help=".npy output with shape [blocks,64,13,2]")
    # Each worker owns an OpenFHE context and key material; using every CPU
    # thread would overcommit memory.  Increase explicitly after measuring.
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--alpha", default="0.25")
    ap.add_argument("--offset", default="1e-4")
    ap.add_argument("--degree", default="511")
    ap.add_argument("--mel-norm-file", required=True)
    args = ap.parse_args()

    inputs = sorted(Path(args.inputs).glob("*.bin"))
    if not inputs:
        raise SystemExit(f"no input blocks found in {args.inputs}")
    out_dir = Path(args.output).with_suffix("")
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    env = {
        "E52_POWER_ALPHA": args.alpha,
        "E52_POWER_OFFSET": args.offset,
        "E52_ZERO_CENTER": "1",
        "E52_CHEB_DEGREE": args.degree,
        "E52_MEL_NORM_FILE": args.mel_norm_file,
    }
    for i, inp in enumerate(inputs):
        jobs.append((args.binary, inp, out_dir / f"{i:06d}.out", env))

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        outputs = list(pool.map(run_one, jobs))
    feat = np.stack([np.fromfile(p, dtype=np.float32).reshape(64, 13, 2)
                     for p in outputs])
    np.save(args.output, feat)
    print(f"blocks={len(inputs)} workers={args.workers} output={args.output} "
          f"shape={feat.shape}")


if __name__ == "__main__":
    main()
