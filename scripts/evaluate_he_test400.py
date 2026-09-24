#!/usr/bin/env python3
"""Compare packed HE features with plaintext reference and run downstream ID."""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


def dct13(x):
    n = np.arange(80, dtype=np.float64)
    k = np.arange(13, dtype=np.float64)[:, None]
    return np.einsum("...n,kn->...k", x,
                     np.cos(np.pi * k * (n + 0.5) / 80.0))


def power_dct(mel, norms, alpha, offset):
    shift = offset ** alpha
    c = np.maximum(mel / norms + offset, 0.0) ** alpha - shift
    return dct13(c * norms ** alpha).astype(np.float32)


class Net(nn.Module):
    def __init__(self, n_spk):
        super().__init__()
        def block(i, o):
            return nn.Sequential(nn.Conv2d(i, o, 3, padding=1),
                                 nn.BatchNorm2d(o), nn.ReLU(),
                                 nn.MaxPool2d(2))
        self.body = nn.Sequential(block(1, 32), block(32, 64),
                                   block(64, 128), nn.AdaptiveAvgPool2d(1),
                                   nn.Flatten())
        self.head = nn.Linear(128, n_spk)

    def forward(self, x):
        return self.head(self.body(x))


def train_eval(train_x, train_y, test_x, test_y, n_spk, seed):
    mu, sd = train_x.mean(), train_x.std()
    tr = torch.from_numpy(((train_x - mu) / sd).astype(np.float32))[:, None]
    te = torch.from_numpy(((test_x - mu) / sd).astype(np.float32))[:, None]
    yt = torch.from_numpy(train_y.astype(np.int64))
    torch.manual_seed(seed)
    net = Net(n_spk)
    opt = torch.optim.Adam(net.parameters(), 1e-3)
    gen = torch.Generator().manual_seed(seed)
    for _ in range(8):
        order = torch.randperm(len(tr), generator=gen)
        net.train()
        for i in range(0, len(order), 256):
            b = order[i:i + 256]
            opt.zero_grad()
            nn.functional.cross_entropy(net(tr[b]), yt[b]).backward()
            opt.step()
    net.eval()
    with torch.no_grad():
        pred = net(te).argmax(1).numpy()
    return float(np.mean(pred == test_y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feat", required=True)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--he", required=True)
    ap.add_argument("--plain", required=True)
    ap.add_argument("--train-dct", required=True,
                    help="uniform-filterbank plaintext DCT tensor for all rows")
    ap.add_argument("--norms", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--alpha", type=float, default=0.25)
    ap.add_argument("--offset", type=float, default=1e-4)
    args = ap.parse_args()

    z = np.load(args.feat)
    X, y, u = z["X"], z["y"], z["u"]
    n_spk = int(z["n_spk"])
    segs = json.loads(Path(args.segments).read_text())
    test_u = np.array([s["utterance"] for s in segs], dtype=np.int64)
    test_u = np.unique(test_u)
    te = np.isin(u, test_u)
    tr = ~te
    norms = np.loadtxt(args.norms)[:80].astype(np.float64)

    train_dct_all = np.load(args.train_dct, mmap_mode="r")
    if train_dct_all.shape[0] != len(X):
        raise ValueError("train DCT row count does not match feat.npz")
    train_dct = train_dct_all[tr]
    plain = np.load(args.plain)[:, :, :, 0].astype(np.float32)
    he = np.load(args.he)
    he = he.astype(np.float64)
    frames = np.empty((he.shape[0] * he.shape[1] * 2, 13), dtype=np.float64)
    for b in range(he.shape[0]):
        for lane in range(he.shape[1]):
            frames[b * 128 + 2 * lane] = he[b, lane, :, 0]
            frames[b * 128 + 2 * lane + 1] = he[b, lane, :, 1]
    he_dct = frames[:len(segs) * 100].reshape(len(segs), 100, 13)
    plain_dct = plain
    test_y = np.array([s["speaker"] for s in segs], dtype=np.int64)
    # speaker ids in the cache are already dense, but preserve this invariant.
    if test_y.max() >= n_spk:
        raise ValueError("test speaker id outside cached label range")

    delta = he_dct.astype(np.float64) - plain_dct.astype(np.float64)
    result = {
        "shape_he": list(he.shape),
        "n_test_segments": len(segs),
        "n_train_segments": int(tr.sum()),
        "rmse": float(np.sqrt(np.mean(delta * delta))),
        "max_abs": float(np.max(np.abs(delta))),
        "rmse_even": float(np.sqrt(np.mean(delta[:, 0::2] * delta[:, 0::2]))),
        "rmse_odd": float(np.sqrt(np.mean(delta[:, 1::2] * delta[:, 1::2]))),
        "per_coeff_rmse": [float(np.sqrt(np.mean(delta[:, :, k] ** 2)))
                           for k in range(13)],
    }
    # The 400 test segments are selected explicitly and are not rows of X in
    # general, so use their speaker labels and the separately generated plain
    # reference as the held-out evaluation set.
    accs = []
    for seed in range(5):
        accs.append(train_eval(train_dct, y[tr], plain_dct,
                               test_y, n_spk, seed))
    he_accs = []
    for seed in range(5):
        he_accs.append(train_eval(train_dct, y[tr], he_dct,
                                  test_y, n_spk, seed))
    result.update({
        "plaintext_test_accs": accs,
        "plaintext_test_mean": float(np.mean(accs)),
        "plaintext_test_std": float(np.std(accs)),
        "he_test_accs": he_accs,
        "he_test_mean": float(np.mean(he_accs)),
        "he_test_std": float(np.std(he_accs)),
    })
    Path(args.output).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
