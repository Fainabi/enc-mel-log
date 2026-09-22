# Downstream task check for the power compressor

This benchmark measures the downstream effect of the exact floating-point
power transform, before CKKS approximation and encryption overhead. It uses
the cached LibriSpeech `dev-clean` features, the fixed utterance-level split,
the same 40-way speaker-ID CNN, eight training epochs, and five seeds used by
the exact-log reference.

| input transform | mean accuracy | std. dev. | per-seed accuracy |
|---|---:|---:|---|
| exact log reference | 0.98399 | 0.00675 | 0.97369, 0.98992, 0.99132, 0.97873, 0.98629 |
| exact `x^0.75` | 0.69247 | 0.12808 | 0.63476, 0.75679, 0.47439, 0.75007, 0.84635 |

The `x^0.75` row is a task-level measurement of the compressor semantics. It
is not yet an end-to-end encrypted-task result: the CKKS degree-511 frontend
has been validated against the real-valued frontend on a 128-frame sample,
while running encrypted inference over the complete downstream dataset would
require a separate batching and decryption pipeline.

