#!/usr/bin/env python3
"""Five-seed downstream speaker-ID check for one exact power compressor."""
import json
import os
import time

import numpy as np
import torch
import torch.nn as nn


DATA = os.environ.get("E27_DATA", "feat.npz")
OUT = os.environ.get("E27_TASK_OUT", "alpha_task.json")
ALPHA = float(os.environ.get("E27_TASK_ALPHA", "0.75"))
THREADS = int(os.environ.get("E27_THREADS", "48"))
SEEDS = [0, 1, 2, 3, 4]


def build(n):
    def block(i, o):
        return nn.Sequential(nn.Conv2d(i, o, 3, padding=1), nn.BatchNorm2d(o),
                             nn.ReLU(), nn.MaxPool2d(2))
    return nn.Sequential(block(1, 32), block(32, 64), block(64, 128),
                         block(128, 128), nn.AdaptiveAvgPool2d(1), nn.Flatten(),
                         nn.Linear(128, n))


def run(A, y, tr, te, n_spk, seed):
    mu, sd = A[tr].mean(), A[tr].std()
    atr = torch.from_numpy(((A[tr] - mu) / sd).astype(np.float32)).unsqueeze(1)
    ate = torch.from_numpy(((A[te] - mu) / sd).astype(np.float32)).unsqueeze(1)
    yt = torch.from_numpy(y[tr])
    torch.manual_seed(seed)
    net = build(n_spk)
    opt = torch.optim.Adam(net.parameters(), 1e-3)
    loss = nn.CrossEntropyLoss()
    gen = torch.Generator().manual_seed(seed)
    for _ in range(8):
        net.train()
        order = torch.randperm(len(atr), generator=gen)
        for i in range(0, len(order), 256):
            b = order[i:i + 256]
            opt.zero_grad()
            loss(net(atr[b]), yt[b]).backward()
            opt.step()
    net.eval()
    good = 0
    with torch.no_grad():
        for i in range(0, len(ate), 512):
            good += (net(ate[i:i + 512]).argmax(1).numpy() == y[te][i:i + 512]).sum()
    return float(good) / len(y[te])


def main():
    torch.set_num_threads(THREADS)
    z = np.load(DATA)
    X, y, u, n_spk = z["X"], z["y"], z["u"], int(z["n_spk"])
    rng = np.random.default_rng(0)
    ids = np.unique(u)
    rng.shuffle(ids)
    te = np.isin(u, ids[:len(ids) // 5])
    tr = ~te
    median = float(np.median(X[tr][X[tr] > 0]))
    hi = float((X[tr] / median).max())
    V = (X / median / hi).astype(np.float64)
    A = np.power(V, ALPHA).astype(np.float32)
    t0 = time.time()
    accs = [run(A, y, tr, te, n_spk, seed) for seed in SEEDS]
    result = {"alpha": ALPHA, "seeds": SEEDS, "threads": THREADS,
              "n_train": int(tr.sum()), "n_test": int(te.sum()),
              "median": median, "hi": hi, "accs": accs,
              "mean": float(np.mean(accs)), "std": float(np.std(accs)),
              "elapsed_s": time.time() - t0}
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
