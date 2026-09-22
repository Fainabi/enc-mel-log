# CKKS validation for the quarter-power compressor

The encrypted frontend was run on the same 128-frame real-speech input used
by the existing endpoint benchmark. The CKKS parameters, degree, packing, and
linear stages were unchanged; only `E52_POWER_ALPHA=0.25` and the public
offset were varied.

| offset | degree | runtime (s) | throughput (frame/s) | endpoint RMSE | endpoint max error |
|---:|---:|---:|---:|---:|---:|
| `1e-6` | 511 | 41.817 | 3.061 | `4.9998e-02` | `1.0711e-01` |
| `1e-4` | 511 | 41.085 | 3.115 | `4.5145e-03` | `8.1277e-03` |

Both runs consumed the same multiplicative depth as the existing degree-511
path. The larger error at `1e-6` was initially caused by
an implementation bug: zero-centering subtracted the polynomial's value at
zero instead of the exact public value `offset^(1/4)`. After replacing that
subtraction with the exact offset image, `1e-4` gives a usable real-speech
frontend result. The logs are retained next to this record.
