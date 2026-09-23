#!/usr/bin/env python3
"""Prepare an utterance-disjoint real-speech block set for the CKKS frontend."""
import argparse
import json
from pathlib import Path

import numpy as np
import soundfile as sf


SR, WIN, HOP, FRAMES, LANES = 16000, 400, 160, 100, 64


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="LibriSpeech/dev-clean")
    ap.add_argument("--feat", required=True, help="matching feat.npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-speaker", type=int, default=10)
    args = ap.parse_args()

    out = Path(args.out)
    blocks = out / "blocks"
    blocks.mkdir(parents=True, exist_ok=True)
    z = np.load(args.feat)
    y, u = z["y"], z["u"]
    utterance_ids = np.unique(u)
    rng = np.random.default_rng(0)
    rng.shuffle(utterance_ids)
    test_u = set(utterance_ids[:len(utterance_ids) // 5].tolist())

    files = sorted(Path(args.data).glob("*/*/*.flac"))
    if len(files) != len(utterance_ids):
        raise RuntimeError(f"file count {len(files)} != feat utterances {len(utterance_ids)}")
    selected = []
    per_spk = {}
    for uid, path in enumerate(files):
        if uid not in test_u:
            continue
        spk = int(y[np.flatnonzero(u == uid)[0]])
        if per_spk.get(spk, 0) >= args.per_speaker:
            continue
        x, sr = sf.read(path, dtype="float32")
        if sr != SR or len(x) < WIN + HOP * (FRAMES - 1):
            continue
        nseg = (len(x) - WIN) // HOP + 1
        nseg //= FRAMES
        for seg in range(min(nseg, args.per_speaker - per_spk.get(spk, 0))):
            frames = np.stack([x[seg * FRAMES * HOP + j * HOP:
                               seg * FRAMES * HOP + j * HOP + WIN]
                               for j in range(FRAMES)])
            selected.append({"speaker": spk, "utterance": uid,
                             "path": str(path), "segment": seg,
                             "frames": frames})
            per_spk[spk] = per_spk.get(spk, 0) + 1
        if len(per_spk) == len(set(y.tolist())) and all(
                v >= args.per_speaker for v in per_spk.values()):
            break
    if len(selected) != len(set(y.tolist())) * args.per_speaker:
        raise RuntimeError(f"selected {len(selected)} segments, per-speaker={per_spk}")

    # Complex packing: two real frames share one complex lane.
    all_frames = np.concatenate([s["frames"] for s in selected], axis=0)
    total = len(all_frames)
    block_count = (total + 2 * LANES - 1) // (2 * LANES)
    manifest = []
    for b in range(block_count):
        block = np.zeros((WIN, LANES, 2), dtype=np.float32)
        for lane in range(LANES):
            for comp in range(2):
                fi = b * 2 * LANES + 2 * lane + comp
                if fi >= total:
                    continue
                block[:, lane, comp] = all_frames[fi]
                sidx, fidx = divmod(fi, FRAMES)
                manifest.append({"block": b, "lane": lane, "component": comp,
                                 "segment_index": sidx, "frame": fidx,
                                 "speaker": selected[sidx]["speaker"],
                                 "utterance": selected[sidx]["utterance"],
                                 "segment": selected[sidx]["segment"]})
        block.tofile(blocks / f"{b:06d}.bin")
    with (out / "manifest.jsonl").open("w") as f:
        for row in manifest:
            f.write(json.dumps(row) + "\n")
    with (out / "segments.json").open("w") as f:
        json.dump([{k: v for k, v in s.items() if k != "frames"} for s in selected], f)
    print(f"segments={len(selected)} frames={total} blocks={block_count} output={out}")


if __name__ == "__main__":
    main()
