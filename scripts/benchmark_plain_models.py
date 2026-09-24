#!/usr/bin/env python3
"""Small architecture/normalization screen for the uniform DCT task."""
import argparse
import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


class Net(nn.Module):
    def __init__(self, n_spk, kind, channels=13):
        super().__init__()
        if kind == "temporal1d_delta":
            kind = "temporal1d"
        def block(i, o, pool):
            return nn.Sequential(nn.Conv2d(i, o, 3, padding=1),
                                 nn.BatchNorm2d(o), nn.ReLU(),
                                 nn.MaxPool2d(pool))
        if kind == "temporal1d":
            self.temporal = True
            self.temporal_dual = False
            self.attn = False
            self.body = nn.Sequential(
                nn.Conv1d(channels, 64, 5, padding=2), nn.BatchNorm1d(64), nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(64, 128, 5, padding=2), nn.BatchNorm1d(128), nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(128, 256, 3, padding=1), nn.BatchNorm1d(256), nn.ReLU(),
                nn.AdaptiveAvgPool1d(1), nn.Flatten())
            self.head = nn.Linear(256, n_spk)
            return
        if kind == "temporal1d_dual":
            self.temporal = True
            self.attn = False
            self.body = nn.Sequential(
                nn.Conv1d(channels, 64, 5, padding=2), nn.BatchNorm1d(64), nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(64, 128, 5, padding=2), nn.BatchNorm1d(128), nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(128, 256, 3, padding=1), nn.BatchNorm1d(256), nn.ReLU())
            self.avg = nn.AdaptiveAvgPool1d(1)
            self.max = nn.AdaptiveMaxPool1d(1)
            self.head = nn.Linear(512, n_spk)
            self.temporal_dual = True
            return
        if kind == "temporal1d_attn":
            self.temporal = True
            self.temporal_dual = False
            self.attn = True
            self.body = nn.Sequential(
                nn.Conv1d(channels, 64, 5, padding=2), nn.BatchNorm1d(64), nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(64, 128, 5, padding=2), nn.BatchNorm1d(128), nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(128, 256, 3, padding=1), nn.BatchNorm1d(256), nn.ReLU())
            self.score = nn.Conv1d(256, 1, 1)
            self.head = nn.Linear(256, n_spk)
            return
        self.attn = False
        self.temporal = False
        self.temporal_dual = False
        if kind == "current":
            body = [block(1, 32, 2), block(32, 64, 2), block(64, 128, 2)]
        elif kind == "freq_preserve":
            body = [block(1, 32, (2, 1)), block(32, 64, (2, 1)),
                    block(64, 128, (2, 2)), block(128, 128, (2, 2))]
        elif kind == "temporal":
            body = [block(1, 32, (2, 1)), block(32, 64, (2, 1)),
                    block(64, 128, (2, 1))]
        else:
            raise ValueError(kind)
        self.body = nn.Sequential(*body, nn.AdaptiveAvgPool2d(1),
                                   nn.Flatten())
        self.head = nn.Linear(128, n_spk)

    def forward(self, x):
        if self.temporal:
            z = self.body(x[:, 0].transpose(1, 2))
            if self.attn:
                w = torch.softmax(self.score(z), dim=-1)
                return self.head((z * w).sum(dim=-1))
            if self.temporal_dual:
                return self.head(torch.cat([self.avg(z).flatten(1),
                                            self.max(z).flatten(1)], dim=1))
            return self.head(z)
        return self.head(self.body(x))


def one(train_x, train_y, test_x, test_y, n_spk, kind, norm, seed, epochs,
        return_pred=False, test_extra=None):
    if norm == "global":
        mu, sd = train_x.mean(), train_x.std()
    else:
        mu = train_x.mean(axis=(0, 1), keepdims=True)
        sd = train_x.std(axis=(0, 1), keepdims=True)
    tr = torch.from_numpy(((train_x - mu) / (sd + 1e-8)).astype(np.float32))[:, None]
    te = torch.from_numpy(((test_x - mu) / (sd + 1e-8)).astype(np.float32))[:, None]
    te_extra = None
    if test_extra is not None:
        te_extra = torch.from_numpy(((test_extra - mu) / (sd + 1e-8)).astype(np.float32))[:, None]
    yt = torch.from_numpy(train_y.astype(np.int64))
    torch.manual_seed(seed)
    net = Net(n_spk, kind, train_x.shape[2])
    opt = torch.optim.Adam(net.parameters(), 1e-3, weight_decay=1e-5)
    gen = torch.Generator().manual_seed(seed)
    for _ in range(epochs):
        order = torch.randperm(len(tr), generator=gen)
        net.train()
        for i in range(0, len(order), 256):
            b = order[i:i + 256]
            opt.zero_grad()
            nn.functional.cross_entropy(net(tr[b]), yt[b]).backward()
            opt.step()
    net.eval()
    with torch.no_grad():
        logits = net(te)
        pred = logits.argmax(1).numpy()
        extra_logits = net(te_extra) if te_extra is not None else None
    if return_pred:
        return float(np.mean(pred == test_y)), pred, logits.numpy(), (
            extra_logits.numpy() if extra_logits is not None else None)
    return float(np.mean(pred == test_y))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feat", required=True)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--dct", required=True)
    ap.add_argument("--he")
    ap.add_argument("--output", required=True)
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--kind", default="all")
    ap.add_argument("--norm", default="all")
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--ensemble", action="store_true")
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    z = np.load(args.feat)
    y, u, n_spk = z["y"], z["u"], int(z["n_spk"])
    segs = json.loads(Path(args.segments).read_text())
    test_u = np.unique([s["utterance"] for s in segs])
    tr = ~np.isin(u, test_u)
    X = np.load(args.dct, mmap_mode="r")
    rows = []
    for s in segs:
        rows.append(np.flatnonzero(u == s["utterance"])[s["segment"]])
    test = X[np.asarray(rows)]
    he_test = None
    if args.he:
        raw = np.load(args.he).astype(np.float32)
        frames = np.empty((raw.shape[0] * raw.shape[1] * 2, 13), dtype=np.float32)
        for b in range(raw.shape[0]):
            for lane in range(raw.shape[1]):
                frames[b * 128 + 2 * lane] = raw[b, lane, :, 0]
                frames[b * 128 + 2 * lane + 1] = raw[b, lane, :, 1]
        he_test = frames[:len(segs) * 100].reshape(len(segs), 100, 13)
    if args.kind == "temporal1d_delta":
        def augment(a):
            d1 = np.diff(a, axis=1, prepend=a[:, :1])
            d2 = np.diff(d1, axis=1, prepend=d1[:, :1])
            return np.concatenate([a, d1, d2], axis=2)
        X = augment(np.asarray(X))
        test = augment(test)
        if he_test is not None:
            he_test = augment(he_test)
    results = {}
    norms = ("global", "per_coeff") if args.norm == "all" else (args.norm,)
    kinds = ("current", "freq_preserve", "temporal", "temporal1d", "temporal1d_dual", "temporal1d_attn") if args.kind == "all" else (args.kind,)
    seeds = [int(x) for x in args.seeds.split(",") if x]
    for norm in norms:
        for kind in kinds:
            test_y = np.asarray([s["speaker"] for s in segs])
            runs = [one(X[tr], y[tr], test, test_y, n_spk, kind, norm,
                        seed, args.epochs, args.ensemble, he_test)
                        for seed in seeds]
            vals = [r[0] if args.ensemble else r for r in runs]
            results[f"{kind}+{norm}"] = {
                "accs": vals, "mean": float(np.mean(vals)),
                "std": float(np.std(vals))}
            if args.ensemble:
                votes = np.stack([r[1] for r in runs])
                ens = np.apply_along_axis(lambda x: np.bincount(x, minlength=n_spk).argmax(),
                                          0, votes)
                results[f"{kind}+{norm}"]["ensemble_acc"] = float(np.mean(ens == test_y))
                probs = np.stack([np.exp(r[2] - r[2].max(1, keepdims=True)) /
                                  np.exp(r[2] - r[2].max(1, keepdims=True)).sum(1, keepdims=True)
                                  for r in runs])
                soft = probs.mean(0).argmax(1)
                results[f"{kind}+{norm}"]["soft_ensemble_acc"] = float(np.mean(soft == test_y))
                if args.he:
                    he_logits = np.stack([r[3] for r in runs])
                    he_probs = np.exp(he_logits - he_logits.max(2, keepdims=True))
                    he_probs /= he_probs.sum(2, keepdims=True)
                    he_soft = he_probs.mean(0).argmax(1)
                    results[f"{kind}+{norm}"]["he_soft_ensemble_acc"] = float(np.mean(he_soft == test_y))
            print(kind, norm, vals, flush=True)
    Path(args.output).write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
