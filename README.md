# ppaudio-he-frontend

CKKS 密文域 log-mel 语音前端的 C++ 参考实现（OpenFHE）。

---

## a) 这是什么

这个目录是论文 *Multiplicative Depth Budgets for Encrypted Log-Mel Front-Ends
Are Set by Dynamic Range, Not Approximability*（ICASSP 2027 投稿）里 C++ 部分的
全部源码：一条在 CKKS 密文上跑完的 log-mel 前端链，以及围绕它的若干单元测试和
诊断程序。

主程序 `cpp_complex_e2e` 在一条密文上依次执行：

```
加窗 STFT 投影  →  取模平方  →  mel 滤波器组  →  u^α 压缩器（Chebyshev 近似）  →  DCT
```

其中 STFT 投影用三级 radix-8 分解（512 = 8³）加复数打包实现，模平方靠共轭
automorphism 把实部虚部拆开后各自 `EvalSquare` 再重组。当前 degree-63 原型先将
mel energy 公开归一化到 [0,1]，再用 α=3/4 的幂函数做 Chebyshev 近似；mel 和 DCT 都是密文-明文线性映射，用
baby-step/giant-step 的对角线乘法实现。整条链跑完输出前 13 个 DCT 系数。

**它做什么**：给出这条链在 CKKS 下的一个能跑通的实现，测它的耗时与它对同参数明文链
的偏差，并把耗时拆到 projection / mel / cheb / DCT 四段上。

**它不做什么**：

- 不做分帧。分帧在客户端明文侧完成，这里的输入已经是排好的帧。
- 不做自举（bootstrapping）。乘法深度配的是 20，整条链在这个预算内跑完。
- 不做训练、不做下游任务评测。论文里 LibriSpeech 上的动态范围实验和 40 人说话人
  识别实验是 Python 侧做的，**不在本仓库内**。
- 不接受真实语音输入。计时和精度验证用的是程序内生成的合成正弦信号（见 f 节）。
- 不是产品代码。没有序列化、没有客户端/服务端拆分、没有错误处理。

---

## b) 依赖与部署

### 依赖清单

| 项 | 要求 | 说明 |
| --- | --- | --- |
| OpenFHE | 论文写的是 **1.4.0** | 必须是带 `CKKSDataType` / `SetCKKSDataType(COMPLEX)` 的版本，见下方「版本核实」 |
| CMake | ≥ 3.16 | `CMakeLists.txt` 里写死的下限 |
| C++ 编译器 | 支持 C++17 | 本机 `/usr/bin/c++`；之前构建 OpenFHE 用的是 GCC 12.3.1 / 14.2.1 |
| OpenMP | 必需 | 不是本项目直接调用，而是 OpenFHE 自己编译时开的 `-fopenmp`，会通过 `OpenFHE_CXX_FLAGS` 和 `OpenFHE_SHARED_LIBRARIES` 传进来 |
| MATHBACKEND | 4 | 由 OpenFHE 的构建决定，通过 `-DMATHBACKEND=4` 传进来 |
| NATIVE_SIZE | 64 | 同上 |

内存：论文记录整条链的峰值 RSS 在 `Dtot=19` 时是 3.67 GB，本仓库的
`cpp_complex_e2e`（深度 20、N=2^16、511 个旋转密钥）请按 **至少 8 GB 可用内存** 准备。

### 版本核实（2026-09-17，读文件得出，未编译验证）

论文 §5 写的是 OpenFHE 1.4.0。**本机装的不是这个版本。** 逐条证据：

| 来源 | 内容 |
| --- | --- |
| `build_complex_upstream/CMakeCache.txt` | `OpenFHE_DIR=/home/faina/projects/openfhe-development/build`、`OpenFHE_INCLUDE=/home/faina/projects/openfhe-development/src`、`OpenFHE_LIBDIR=/home/faina/projects/openfhe-development/build/lib` |
| `build_diag/CMakeCache.txt` | `OpenFHE_DIR` 同上 |
| `/home/faina/projects/openfhe-development/build/OpenFHEConfig.cmake` | `BASE_OPENFHE_VERSION 1.2.4`；`OpenFHE_INCLUDE=/usr/local/include/openfhe`；`OpenFHE_LIBDIR=/usr/local/lib`；`OpenFHE_CXX_COMPILER_VERSION 12.3.1`；`-DMATHBACKEND=4 -fopenmp`；`NATIVE_SIZE 64` |
| `openfhe-development` git | HEAD 停在 tag `v1.2.4`（commit `6bcca75`，2025-03-21） |
| 另一份 checkout `/home/faina/projects/openfhe-ckks-btp` | `CMakeLists.txt` 里 `OPENFHE_VERSION = 1.2.4` |
| `.so` 文件 | 只有 `libOPENFHE{core,pke,binfhe}.so.1.2.4`，在 `openfhe-development/build/lib`（2025-05）和 `/tmp/ofhe-native/lib`（2026-09-12） |
| 系统安装 | **没有**。`/usr/local/include` 这个目录不存在，`/usr/local/lib` 下也没有任何 OpenFHE 库 |

由此得到三条要记住的事：

1. **本机能找到的 OpenFHE 全部是 1.2.4，与论文声明的 1.4.0 不一致。**
2. **OpenFHE 没有装到系统里**，只有构建树。而构建树里那份
   `OpenFHEConfig.cmake` 导出的 `OpenFHE_INCLUDE` / `OpenFHE_LIBDIR` 指向的是
   *安装后* 的路径 `/usr/local/...`，这些路径当前不存在。所以直接
   `find_package(OpenFHE)` 拿到的头文件路径是空的——`results/` 下旧构建目录的
   `flags.make` 里 `CXX_INCLUDES = -I/usr/local/include/openfhe ...` 就是这么来的。
   **正确做法是把 OpenFHE 真正 `make install` 一遍**，见下方部署步骤。
3. **本仓库源码用到的两个 API 在 1.2.4 里不存在**（对两份 1.2.4 checkout 全文 grep 均无命中）：
   - `CCParams::SetCKKSDataType(...)` —— `src/cpp_complex_e2e.cpp:221`，8 个源文件全部用到
   - `EvalChebyshevFunctionPtxt(...)` —— `src/cpp_complex_e2e.cpp:322`（`math/chebyshev.h`）

   所以 **拿 1.2.4 编不过本仓库**，必须用带复数 CKKS 打包的版本（论文写的 1.4.0）。
   *待确认*：`results/` 里的历史产物（一个 2026-09-17 构建的
   `/tmp/cpp_complex_e2e`）却是链接到 1.2.4 的 `.so` 上的，且符号全部解析成功。
   它当时用的头文件目录（`/usr/local/include/openfhe`）现在已经不在盘上，无法复原，
   因此「实际编出那份二进制的 OpenFHE 到底是哪个版本」这一点本 README 无法定论。
   以论文的 1.4.0 为准，首次部署请自行确认。

### 从零部署

```bash
# 1. 取 OpenFHE（版本以论文为准：1.4.0）
git clone --branch v1.4.0 https://github.com/openfheorg/openfhe-development.git
cd openfhe-development

# 2. 配置。WITH_OPENMP 必须开，NATIVE_SIZE=64，MATHBACKEND 用默认的 4。
#    CMAKE_INSTALL_PREFIX 可以换成任意你有写权限的目录，比如 $HOME/.local。
cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/usr/local \
      -DWITH_OPENMP=ON \
      -DBUILD_SHARED=ON \
      -DBUILD_UNITTESTS=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_BENCHMARKS=OFF \
      -DNATIVE_SIZE=64

# 3. 编译（很慢，几十分钟量级）
cmake --build build -j"$(nproc)"

# 4. 安装。装到 /usr/local 需要 root；装到 $HOME/.local 不需要。
#    这一步不能跳过——构建树里的 OpenFHEConfig.cmake 导出的头文件路径
#    指向的是安装后的位置，不装的话 find_package 会给出一个不存在的目录。
sudo cmake --install build          # 或 cmake --install build（前缀在家目录时）

# 5. 回到本项目
cd /path/to/he-frontend
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

### 让 `find_package(OpenFHE)` 找得到

`find_package(OpenFHE REQUIRED)` 找的是 `OpenFHEConfig.cmake`。三种指法，任选其一：

```bash
# (1) 装到了标准前缀（/usr/local）：通常不用做任何事

# (2) 装到了非标准前缀
cmake -S . -B build -DCMAKE_PREFIX_PATH=$HOME/.local

# (3) 直接指到含 OpenFHEConfig.cmake 的目录
#     （安装后是 <prefix>/lib/OpenFHE，未安装的构建树是 <openfhe>/build）
cmake -S . -B build -DOpenFHE_DIR=$HOME/.local/lib/OpenFHE
```

configure 成功时会打印：

```
-- OpenFHE include : .../include/openfhe
-- OpenFHE libdir  : .../lib
-- OpenFHE version : 1.4.0
-- targets         : cpp_complex_e2e;cpp_complex_ops;...
```

如果 `OpenFHE include` 指向的目录不存在，说明 OpenFHE 只构建未安装，回到第 4 步。

**关于头文件路径的一处处理**：`src/cpp_complex_e2e.cpp` 写的是
`#include "math/chebyshev.h"`，而 `src/cpp_real_dct.cpp` 写的是
`#include "openfhe/core/math/chebyshev.h"`——同一个头文件的两种写法，需要两个不同的
include 根。本仓库的 `CMakeLists.txt` 因此同时加入了 `${OpenFHE_INCLUDE}/core` 和
`${OpenFHE_INCLUDE}` 的父目录。源码本身没有改（见下）。

### 本 README 未验证的部分

> **本目录的 CMake 配置未在本机实际验证过编译。** 没有跑过 `cmake configure`，
> 也没有跑过 `make`。首次使用请自行 `cmake -S . -B build` 一次，确认
> OpenFHE 被正确找到、8 个 target 都能配置出来，再开始编译。

---

## c) 编译本项目

out-of-source build：

```bash
cd he-frontend
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

只编主程序：

```bash
cmake --build build --target cpp_complex_e2e -j"$(nproc)"
```

产物全部落在 `build/`，名字就是 target 名。清理直接 `rm -rf build`。

注意编译选项来自 OpenFHE 的 `OpenFHE_CXX_FLAGS`，其中含 `-Wall -Werror`。换一个
编译器版本时可能因为新的告警而在 `-Werror` 上失败，这时把该告警加进
`OPENFHE_FLAGS` 的抑制列表，或者临时 `-DCMAKE_CXX_FLAGS=-Wno-error`。

---

## d) 8 个可执行文件各是干什么的

以下全部按源码逐个核对。除主程序可选的 binary waveform 入口外，输入在程序内部生成；
所有输出都走 stdout。「耗时」一栏里，密文求值时间有日志佐证的地方给了实测值，
其余是按参数（密钥生成规模、明文矩阵预计算规模）给的量级判断，标为「估」。

### 1. `cpp_complex_e2e` — 整条链的主程序

| | |
| --- | --- |
| 源码 | `src/cpp_complex_e2e.cpp` |
| 跑什么 | 完整链：三级 radix-8 DFT 投影 → 共轭拆分实虚部 → 各自平方 → 重组 → mel → 归一化到 [0,1] → 再拆分 → 63 次 Chebyshev 的 u^(3/4) → 重组 → DCT，最后解密前 13 个系数并与浮点 reference 比较 |
| 输入 | 默认生成 64 个合成 complex frame。设置 `E52_INPUT_BIN` 后读取 little-endian float32：`WIN*LANES` 个值表示 real frame，`2*WIN*LANES` 个值表示每个 complex slot 的 real/imag 两个 frame |
| 环境变量 | `E52_BENCH` 跳过逐级 probe；`E52_INPUT_BIN` 使用外部 binary 输入；`E52_MEL_NORM_FILE` 读取 80 个公开 per-channel calibration scale；`E52_POWER_ALPHA`、`E52_POWER_OFFSET`、`E52_CHEB_DEGREE` 覆盖 nonlinear prototype 参数 |
| 输出 | 逐级 `stage <name> level=.. rmse=.. maxerr=..`；`input <name> min/max/mean/std/neg/over1` 统计；13 行 `dct_row k=..`；`runtime_s=.. level=.. rmse=.. maxerr=..`；`bench_stage_s projection=.. mel=.. cheb=.. dct=..`；4 行 `slot<i>=(..) exp=(..)` |
| 耗时 | **密文求值 43.800 s**（有日志：`../results/cpp_complex_fft3_64.log`）。这个数字只覆盖求值，**不含** 建上下文、KeyGen（511 个正向 + 12 个反向旋转密钥 + 共轭 automorphism 密钥）、以及 CPU 侧的 512×512 稠密复矩阵连乘（`mm()`，`:29-36`，共调用 13 次）与全部明文对角线编码。整个进程的实际运行时间明显长于 43.8 s，本仓库没有测过这个总时间 |
| 备注 | 输出的 `level=` 是 OpenFHE 的内部层索引，不是论文里说的「剩余层数」 |

### 2. `cpp_complex_ops` — 复数原语的 smoke test

| | |
| --- | --- |
| 源码 | `src/cpp_complex_ops.cpp` |
| 跑什么 | 逐个验证主程序依赖的四个手写原语：乘 i（乘单项式 X^(M/4)，`:30-42`）、除以 2（在每个 RNS 模数上乘 2 的模逆，`:18-29`）、共轭 automorphism、由此得到的实虚部拆分 → 各自 `EvalSquare` → 乘 i 重组 |
| 输入 | 两个槽，`x[0]=(0.1,0.2)`、`x[1]=(0.3,0.4)`（`:60-61`） |
| 输出 | 每一步一行 `<name>=(re,im),(re,im) level=..`：`in` / `mul_i` / `re` / `im` / `r2` / `i2` / `i2i` / `q` |
| 耗时 | 估：秒到数十秒。参数是 depth 20 / N=2^16（`:45-47`），成本主要在建上下文和 KeyGen；**没有生成旋转密钥**，所以比主程序快得多 |

### 3. `cpp_complex_smoke` — 复数 CKKS 是否真的可用

| | |
| --- | --- |
| 源码 | `src/cpp_complex_smoke.cpp` |
| 跑什么 | 全槽（16384 个）复数密文 × 复数明文乘法 + rescale，比较与明文逐元素乘的最大偏差 |
| 输入 | 程序内生成 `x[i]=(0.001 sin(0.01i), 0.001 cos(0.013i))`、`m[i]=(0.7+0.1 sin(0.007i), 0.2 cos(0.011i))`（`:26-30`） |
| 输出 | `plaintext m[0]=..`、`input got/expected`、`sample got/expected`、`slots=16384 max_complex_mult_error=..` |
| 退出码 | **8 个 target 里唯一有 pass/fail 退出码的**：`maxerr < 1e-5` 返回 0，否则返回 1（`:47`） |
| 耗时 | 估：数秒。参数最轻（depth 5 / N=2^15 / scale 50，`:13-16`），无旋转密钥 |
| 用途 | 换 OpenFHE 版本后第一个该跑的，用来确认 `SetCKKSDataType(COMPLEX)` 这条路在当前版本上通 |

### 4. `cpp_dct_small` — DCT 有符号对角线实现的单元测试

| | |
| --- | --- |
| 源码 | `src/cpp_dct_small.cpp` |
| 跑什么 | 把 13×80 的 DCT 矩阵编成 159 条有符号对角线（d = -79..79，`:71-79`），用 `mv_direct`（`:34-46`）求值，对两种输入各跑一次：`fresh`（新鲜密文）和 `deep`（先做 10 轮乘 1 + rescale 把密文压到第 10 层，`:86-93`）。目的是看这条 DCT 路在链的深处还准不准 |
| 输入 | 80 个槽的指数衰减复数值（`:61-63`） |
| 输出 | 每个 case 一行 `small_dct case=fresh|deep input_level=.. output_level=.. rmse=.. max=..`，外加 13 行逐系数 `k<i> got=(..) exp=(..)` |
| 耗时 | 估：分钟量级。depth 20 / N=2^16，KeyGen 要生成 ±1..±79 共 158 个旋转密钥（`:56-58`） |
| 备注 | 槽数 `S` 是写死的 32768（`:9`），不是从上下文查的；改 `SetRingDim` 必须同步改它 |

### 5. `cpp_inv2_placement` — 决定「÷2 放在哪一步」的诊断程序

| | |
| --- | --- |
| 源码 | `src/cpp_inv2_placement.cpp`（单行密排风格，整个 `main` 是一行） |
| 跑什么 | 共轭拆分 re = (z + conj z)/2 里那个 1/2 的四种放法对比，外加两个前置的正确性检查：<br>① `dcrt_roundtrip`：多项式乘 2⁻¹ 再乘 2 是否逐系数复原<br>② `ciphertext_roundtrip`：整条密文上做同样的往返是否复原<br>③ 四种放法：`scalar_after_add`（`EvalMult(.5)`）、`inv2_before_conj`、`inv2_after_add`、`inv2_times2` |
| 输入 | 两个槽，`x[0]=(0.1,0.2)`、`x[1]=(0.3,0.4)` |
| 输出 | `dcrt_roundtrip=PASS|FAIL`、`ciphertext_roundtrip=PASS|FAIL`，以及每种放法一行 `<name> level=.. scale=.. noiseDeg=.. [两个槽的值]` |
| 耗时 | 估：秒到数十秒。depth 8 / N=2^16，无旋转密钥，无 `EvalMultKeyGen` |
| 备注 | 源码里的 `printf` 用的是 `"\\n"`（转义了反斜杠），所以输出**不换行**，字面打出 `\n`。这是原始行为，本仓库没有改 |

### 6. `cpp_linear_smoke` — OpenFHE 内置线性变换的对照

| | |
| --- | --- |
| 源码 | `src/cpp_linear_smoke.cpp` |
| 跑什么 | 用 OpenFHE 自带的 `EvalLinearTransformPrecompute` / `EvalLinearTransform` 算一个 8×8 复矩阵乘向量，跟明文结果比。这是手写 BSGS `mv()` 的对照基线 |
| 输入 | 8×8 矩阵 `A[i][j] = (i==j ? 1.0 : 0.01*(i+j))`，向量 `x[i] = (i+1, 10+i)`（`:10-11`） |
| 输出 | 8 行 `<i> <got> exp <expected>`（`std::cout` 打 `std::complex`，形如 `(re,im)`） |
| 耗时 | 估：数秒到数十秒。depth 5 / N=2^15，旋转密钥只有 1..7 |

### 7. `cpp_mv_unit` — 手写 BSGS 矩阵-向量乘的单元测试

| | |
| --- | --- |
| 源码 | `src/cpp_mv_unit.cpp` |
| 跑什么 | 测 `mv()`（n1=16 的 baby-step/giant-step，`:9-41`，与主程序里的同名函数同构）两件事：<br>① 一个只有两个非零位置的挑选矩阵，输出跟直接 `EvalRotate(ct,1)` 对照，用来验 giant-step 的对齐；<br>② 13×80 的 DCT 编成**未预移位**的普通对角线，看结果 |
| 输入 | 三个槽，`x[0]=0.2`、`x[1]=0.3`、`x[17]=0.7`（`:60-62`） |
| 输出 | `out=(..),(..)  level=..`、`rotate=(..),(..)`、`dct unit=(..),(..)` |
| 耗时 | 估：分钟量级。depth 只有 4，但 KeyGen 要生成 **511 个旋转密钥**（`:55-56`），这是主要开销 |
| 参考值 | `../results/cpp_mv_unit_final.log`：`out=(1.52e-14,1.17e-14),(1,4.21e-13) level=1`，`dct unit=(0.5,..),(0.29947968,..)` |

### 8. `cpp_real_dct` — REAL 数据类型的对照路径

| | |
| --- | --- |
| 源码 | `src/cpp_real_dct.cpp` |
| 跑什么 | 用 `SetCKKSDataType(REAL)`（`:20`，**8 个里唯一一个不用 COMPLEX 的**）跑「Chebyshev 压缩器 + DCT」这一段：对 80 个实数先做 63 次 Chebyshev 的 u^(1/4)，再用显式的「旋转 + 掩码对角线累加」做 DCT，解密前 13 个系数跟明文比。此项保留为历史对照，当前可用的 complex full-chain 原型见下一节 |
| 输入 | 80 个实数，`x[i] = 0.08 + 0.0002*i`（`:31-32`） |
| 输出 | 一行 `real_dct level=.. rmse=.. maxerr=.. s0=.. exp=.. s1=.. exp=..` |
| 耗时 | 估：分钟量级。depth 20 / N=2^16 + 511 个旋转密钥 |
| 参考值 | `../results/cpp_real_dct.log`：`real_dct level=9 rmse=1.146732e+00 maxerr=1.859029e+00 s0=43.550815 exp=43.548637 s1=-0.83751568 exp=-0.40207664`。**注意这条路的 RMSE 是 1.15，很差**——k=0 那一项对得上，高阶系数对不上。这是一条被记录下来的失败对照，不是可用结果 |

---

## e) 关键参数（`src/cpp_complex_e2e.cpp`）

| 参数 | 值 | 源码行 | 含义 |
| --- | --- | --- | --- |
| ring dimension N | `1 << 16` = 65536 | `:218` | CKKS 环维数。论文 §5 同值。可用槽数 = N/2 = 32768 |
| 槽数 S | `cc->GetRingDimension()/2` = 32768 | `:235` | 复数 CKKS 下每个槽装一个复数 |
| `SetMultiplicativeDepth` | **20** | `:216` | 乘法深度预算。论文测得 N=2^16 / 128-bit 安全 / 默认 HYBRID digit 数下 L_max = 20 |
| scaling modulus 位宽 | **59** | `:217` | `SetScalingModSize(59)`。论文 §5 同值 |
| 安全等级 | `HEStd_128_classic` | `:219` | 128-bit 经典安全 |
| scaling technique | `FIXEDMANUAL` | `:220` | 手动 rescale：每个矩阵-向量乘只 rescale 一次（`mv()` 末尾，`:180`），而不是每次乘法都 rescale。论文记这比默认快 2.43× |
| CKKS 数据类型 | `COMPLEX` | `:221` | 复数打包。这是需要新版 OpenFHE 的地方 |
| DFT 长度 | `N = 512` | `:13` | 变量名也叫 `N`，与环维数同名但含义不同——这里是 DFT/矩阵的边长 |
| 窗长 WIN | `400` | `:13` | Hann 窗长度，`0.5 - 0.5 cos(2πj/400)`（`:60`） |
| BLOCK | `1024` | `:13` | 每个 block 占 1024 个槽 = 两条 512 槽的帧 lane |
| 帧数 FRAMES | `64` | `:13` | 64 帧 × 10 ms hop = 0.64 s 音频 |
| BLOCKS | `FRAMES/2 = 32` | `:14` | **槽位布局：32 blocks × 2 lanes × 512 slots = 32768 槽，正好填满** |
| BSGS baby-step 数 `MV_N1` | `16` | `:15` | `mv()` 的 n1。giant-step 数 = 512/16 = 32 |
| radix-8 分组 | 3 组，每组 3 级蝶形 | `:54-57` | 9 级 radix-2 蝶形按 3 级一组合并成 1 个矩阵，等价于三级 radix-8。加上最后的「位反转 + 复数转实数交织」矩阵（`:65-78`），共 4 个明文矩阵 |
| mel 滤波器数 | **80** | `:83` | `for (int r = 0; r < 80; r++)` |
| mel 滤波器形状 | `lo = 1+3r, mid = lo+3, hi = mid+3` | `:84` | **均匀 3-bin 间隔的三角形，不是 HTK mel**。见 f 节 |
| mel 全局缩放 | `× (0.125 / MEL_NORM)`，`MEL_NORM=1e-2` | `src/cpp_complex_e2e.cpp:280` | 将 mel 输出归一化到公开的 Chebyshev 定义域 |
| 压缩器指数 α | **0.75** | `src/cpp_complex_e2e.cpp:28` | `std::pow(std::max(z, 0.0), POWER_ALPHA)`，即 u^(3/4) |
| Chebyshev 次数 d | **63** | `:322`（明文参考）、`:370`、`:371`（密文） | 论文选定的工作点是 d=2047，**本源码是 63**。见 f 节 |
| Chebyshev 定义域 | `[0, 1]` | `src/cpp_complex_e2e.cpp:26-27` | |
| DCT 保留系数 | **13** | `:93`（矩阵）、`:388`、`:414`、`:436`（比较与打印） | 13 × 80 的 DCT-II 核，`cos(πk(n+0.5)/80)`（`:95`） |
| 旋转密钥 | 正向 1..511，反向 -1..-12 | `:228-232` | 反向那 12 个给 DCT 的有符号对角线用（d 从 -12 到 79，`:280`） |
| 共轭 automorphism 索引 | `2*ringDim - 1` | `:233` | 实虚部拆分靠它 |
| 启用的功能 | `PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE` | `:223-224` | |

---

## f) 已知局限

**照实列，不要在引用这些数字时省掉这一节。**

1. **mel 不是标准 HTK mel。** `src/cpp_complex_e2e.cpp:81-90` 里的滤波器组是 80 个
   **在频率 bin 上均匀间隔 3 个 bin** 的三角形（`lo = 1+3r, mid = lo+3, hi = mid+3`，
   r=0..79 覆盖 bin 1..243），而论文 §5 描述的前端用的是 80 个 HTK mel 通道。两者不是同一
   个滤波器组。这里的均匀版本只是保留了「稀疏但对角线很多」这个成本结构（论文 §6：
   mel 矩阵 2.45% 稠密却有 257 条对角线里的 178 条非零），足以支撑耗时结论，
   但不能用来评特征质量。

2. **RMSE 只算前 13 个 DCT 系数。** 最终那行 `rmse=` 的求和范围是
   `for (size_t i = 0; i < 13; i++)`（`:436`），不是全部 512 个输出槽，也不是
   32768 个槽。它衡量的是「实际关心的那 13 个 MFCC 系数」的偏差，不是整个密文的偏差。

3. **默认 benchmark 是合成输入；真实语音需要显式 calibration。** `E52_INPUT_BIN` 可以
   注入之前工作的 LibriSpeech `dev-clean` frame，但真实 mel energy 的跨 channel 范围远大于
   合成信号，必须使用同一公开数据上的 `E52_MEL_NORM_FILE`。当前保持 `alpha=0.75`，
   使用 per-channel scale、zero-centering 和 degree-511 时测得 endpoint RMSE
   `7.66e-03`、max error `1.81e-02`；这只是前端一致性检查，不是下游任务精度。

4. **radix-8 分解版的精度比非分解版差约四个半数量级。** 非分解（稠密投影）路径
   在 32 帧布局下测得 `rmse=2.377179e-07`（`../results/cpp_complex_32frame_bench.log`），
   三级 radix-8 分解版在 64 帧布局下是 `rmse=8.436337e-03`
   （`../results/cpp_complex_fft3_64.log`）。比值约 3.5×10⁴，即约 4.5 个数量级。
   分解换来的是明文乘法 512→38、旋转 44→16 的算子数下降，代价就是这个精度损失。
   两次测量的帧数布局不同（32 vs 64），不是严格的同条件对照。

5. **论文选定的工作点 d=2047 从未整链跑通过。** 本仓库源码写死 d=63
   （`:322`、`:370`、`:371`）。论文 §5.3 选的工作点是 d=2047 / 深度 16，那个配置的
   成本（表 2 的 50.2 s/s）是**分阶段算子单价推出来的**，不是端到端测出来的。
   整条链只在 d=63 上实际跑过。这两个数不能混用，也不能相加。

6. **`cpp_real_dct` 这条对照路径是失败的。** 记录的 RMSE 是 1.15
   （`../results/cpp_real_dct.log`），高阶系数完全对不上。保留它是为了记录
   REAL 数据类型路线试过且不行，不是可用结果。

7. **本仓库未在本机验证编译。** 见 b 节末尾。本机可见的 OpenFHE 全是 1.2.4，
   缺源码用到的 API。

---

## g) 复现

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cpp_complex_e2e -j"$(nproc)"
scripts/run_bench.sh
```

`scripts/run_bench.sh` 会把 `cpp_complex_e2e` 的完整输出存到
`logs/cpp_complex_e2e_<YYYYmmdd_HHMMSS>.log`，并在头部记下主机名、时间、二进制路径、
`E52_BENCH` 与 `OMP_NUM_THREADS` 的取值。只要计时不要逐级诊断就
`E52_BENCH=1 scripts/run_bench.sh`。

**跑出来的数字该和什么对得上**，参照日志 `../results/cpp_complex_fft3_64.log`：

| 输出项 | 参照值 |
| --- | --- |
| `runtime_s=` | `43.800`（= 论文的 68.44 s per second of audio，43.800 / 0.64） |
| `rmse=` | `8.436337e-03` |
| `maxerr=` | `8.973544e-03` |
| `level=` | `15` |
| `bench_stage_s projection=` | `22.710301`（35.48 s/s） |
| `bench_stage_s mel=` | `9.747905`（15.23 s/s） |
| `bench_stage_s cheb=` | `5.058857`（7.90 s/s） |
| `bench_stage_s dct=` | `6.220999`（9.72 s/s） |
| `slot0=` | `(7.71655, 7.72113)`，期望 `(7.71952, 7.71494)` |

几点要先说清楚，否则对不上会以为是 bug：

- **RMSE 和各级 level 是确定性的**，同一份 OpenFHE、同一套参数下应当完全一致
  （输入是确定性生成的，CKKS 的噪声来自随机加密，所以 RMSE 会有末位量级的抖动，
  但不该动到有效数字的前两位）。
- **runtime_s 跟机器强相关。** 参照值来自 Intel Xeon Gold 6230R、OpenFHE 默认
  OpenMP 配置。换机器、换线程数，绝对值必然不同；论文里的所有结论建立在
  **同一配置下的相对差异**上，不是绝对值上。
- **runtime_s 只覆盖密文求值**，不含建上下文、KeyGen、编码、加密、解密。
  `run_bench.sh` 额外记的 `wrapper_elapsed_s` 是整个进程的时间，两者差很多是正常的。
- 更多背景（分段归一化、为什么 68.44 而不是别的数、哪些日志是弃用的探索记录）
  见本目录的 `BENCHMARK.md`。原始主结果日志仍在 `../results/cpp_complex_fft3_64.log`。

`results/` 里其它 `.log` 大多是探索期的失败布局和不同帧数的变体，除了上面点名的
几个，不要拿来当参照。

### 真实语音 smoke benchmark

之前工作的 LibriSpeech `dev-clean` 样本可按 128 个连续 frame 组织成一个
`400 x 64` complex binary 输入。先用同一个输入生成 C++ projection 顺序下的公开
calibration scale：

```bash
E52_INPUT_BIN=/path/libri_128frames_complex_f32.bin \
E52_WRITE_MEL_NORMS=/path/mel_norms_cpp.txt \
./build/cpp_complex_e2e
```

再运行真实输入 benchmark：

```bash
E52_BENCH=1 E52_ZERO_CENTER=1 E52_CHEB_DEGREE=511 \
E52_POWER_ALPHA=0.75 E52_POWER_OFFSET=1e-6 \
E52_MEL_NORM_FILE=/path/mel_norms_cpp.txt \
E52_INPUT_BIN=/path/libri_128frames_complex_f32.bin \
./build/cpp_complex_e2e
```

在 Xeon zpf 上该样本的记录为 `42.993 s`、`2.977 frame/s`，endpoint RMSE
`7.66e-03`、max error `1.81e-02`。这个配置保持论文的 `alpha=0.75`；相比默认
degree-63 synthetic prototype，它额外使用 degree-511 和 zero-centering 来处理真实
mel dynamic range。

---

## h) 授权

授权方式待定（TBD）。
