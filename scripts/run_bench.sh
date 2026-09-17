#!/usr/bin/env bash
# 跑一次 cpp_complex_e2e，把 stdout/stderr 存进带时间戳的日志。
#
# 用法:
#   scripts/run_bench.sh [build_dir] [log_dir]
#
# 默认 build_dir = <repo>/build, log_dir = <repo>/logs。
# 环境变量 E52_BENCH 置为任意非空值时, cpp_complex_e2e 会跳过每一级的
# RMSE probe, 只输出最终 runtime/RMSE 与分阶段计时 (源码 cpp_complex_e2e.cpp
# 里的 `bench` 开关)。本脚本默认不设它, 保留完整的逐级诊断输出;
# 想要纯计时就 `E52_BENCH=1 scripts/run_bench.sh`。

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_DIR="${1:-${REPO_DIR}/build}"
LOG_DIR="${2:-${REPO_DIR}/logs}"
BIN="${BUILD_DIR}/cpp_complex_e2e"

if [ ! -d "${BUILD_DIR}" ]; then
  echo "error: build 目录不存在: ${BUILD_DIR}" >&2
  echo "先执行: cmake -S \"${REPO_DIR}\" -B \"${BUILD_DIR}\" -DCMAKE_BUILD_TYPE=Release && cmake --build \"${BUILD_DIR}\" -j" >&2
  exit 1
fi

if [ ! -x "${BIN}" ]; then
  echo "error: 找不到可执行文件或它不可执行: ${BIN}" >&2
  echo "先执行: cmake --build \"${BUILD_DIR}\" --target cpp_complex_e2e -j" >&2
  exit 1
fi

mkdir -p "${LOG_DIR}"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/cpp_complex_e2e_${STAMP}.log"

{
  echo "# host      : $(hostname)"
  echo "# date      : $(date -Iseconds)"
  echo "# binary    : ${BIN}"
  echo "# E52_BENCH : ${E52_BENCH:-<unset>}"
  echo "# OMP_NUM_THREADS : ${OMP_NUM_THREADS:-<unset, OpenFHE 默认>}"
  echo "#"
} > "${LOG_FILE}"

echo "running ${BIN} -> ${LOG_FILE}"
START="$(date +%s)"
set +e
"${BIN}" 2>&1 | tee -a "${LOG_FILE}"
RC="${PIPESTATUS[0]}"
set -e
END="$(date +%s)"

{
  echo "#"
  echo "# exit_code        : ${RC}"
  echo "# wrapper_elapsed_s: $((END - START))"
} >> "${LOG_FILE}"

echo "log: ${LOG_FILE}"
echo "参照值 (results/cpp_complex_fft3_64.log): runtime_s=43.800 rmse=8.436337e-03 maxerr=8.973544e-03"
exit "${RC}"
