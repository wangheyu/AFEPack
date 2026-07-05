# AFEPack 中文安装说明

本文说明当前 AFEPack 源码树的安装、依赖检查、示例编译和测试方式。

## 1. 安装方式概览

当前构建系统基于 Autotools：

- 默认使用串行 `gcc/g++` 编译，不依赖 MPI。
- `make` 会编译库、工具、模板动态库以及 `examples` 下的全部示例程序。
- `make test` 会先执行完整编译，再运行测试。
- MPI 支持为可选功能，仅在配置时显式传入 `--enable-mpi` 后启用。
- 3D 示例测试依赖 `gmsh`。默认 `--enable-gmsh-tests=auto`，缺少 `gmsh` 时跳过 3D 测试；如果传入 `--enable-gmsh-tests=yes`，缺少 `gmsh` 会导致 `configure` 失败并给出安装建议。

如果源码包中已经包含 `configure`，普通用户不需要运行 `autoreconf`。只有在修改过 `configure.ac`、`Makefile.am` 或相关 Autotools 输入文件后，源码维护者才需要重新生成配置脚本：

```sh
autoreconf -fi
```

## 2. Debian/Ubuntu 系统依赖

最小串行构建依赖：

```sh
sudo apt update
sudo apt install build-essential libopenblas-dev liblapacke-dev libboost-serialization-dev
```

如果需要从维护者源码重新生成 `configure` 和 `Makefile.in`，还需要：

```sh
sudo apt install autoconf automake libtool
```

如果需要运行 3D/gmsh 示例测试：

```sh
sudo apt install gmsh
```

如果需要启用 MPI 构建：

```sh
sudo apt install libopenmpi-dev openmpi-bin
```

## 3. 推荐的串行 g++ 构建流程

建议使用源码树外构建目录，避免生成文件污染源码目录：

```sh
mkdir -p build
cd build
CC=gcc CXX=g++ ../configure --prefix="$HOME/.local/afepack"
make -j"$(nproc)"
make test
make install
```

说明：

- `CC=gcc CXX=g++` 用于明确选择非 MPI 编译器。
- 默认不进入 `src/mpi`，也不安装 MPI 头文件。
- `make test` 会先编译所有目标，再运行 `make check`。
- 如果没有安装 `gmsh`，3D 示例测试会显示为 `SKIP`，这不是失败。
- 如果希望缺少 `gmsh` 时直接报错，请在配置时加入 `--enable-gmsh-tests=yes`。

## 4. 自定义依赖安装路径

如果 CBLAS、LAPACKE、OpenBLAS 或 Boost 安装在非系统目录，需要在 `configure` 时指定头文件和库目录。

OpenBLAS 作为组合提供者时：

```sh
CC=gcc CXX=g++ ../configure \
  --prefix="$HOME/.local/afepack" \
  --with-cblas-includedir=/path/to/deps/include \
  --with-lapacke-includedir=/path/to/deps/include \
  --with-openblas-libdir=/path/to/deps/lib \
  --with-boost=/path/to/deps \
  --with-boost-libdir=/path/to/deps/lib
```

如果 CBLAS 和 LAPACKE 来自不同库目录，或者 OpenBLAS 不提供 LAPACKE 符号，请不要使用 `--with-openblas-libdir`，改用：

```sh
CC=gcc CXX=g++ ../configure \
  --prefix="$HOME/.local/afepack" \
  --with-cblas-includedir=/path/to/cblas/include \
  --with-cblas-libdir=/path/to/cblas/lib \
  --with-lapacke-includedir=/path/to/lapacke/include \
  --with-lapacke-libdir=/path/to/lapacke/lib \
  --with-boost=/path/to/boost \
  --with-boost-libdir=/path/to/boost/lib
```

## 5. 配置选项

常用选项：

```text
--prefix=DIR                    安装前缀，默认为 /usr/local
--enable-debug                  使用 -O0 -g 等调试编译选项
--enable-mpi                    启用 MPI 库和 MPI 头文件安装，默认关闭
--enable-gmsh-tests=auto|yes|no 控制 gmsh 相关示例测试，默认 auto
--with-cblas-includedir=DIR     CBLAS 头文件目录
--with-cblas-libdir=DIR         CBLAS 库目录
--with-lapacke-includedir=DIR   LAPACKE 头文件目录
--with-lapacke-libdir=DIR       LAPACKE 库目录
--with-openblas-libdir=DIR      OpenBLAS 库目录
--with-boost=DIR                Boost 安装前缀
--with-boost-libdir=DIR         Boost 库目录
```

查看完整选项：

```sh
../configure --help
```

## 6. 启用 MPI

默认构建是非 MPI 串行构建。如果确实需要 MPI 版本：

```sh
mkdir -p build-mpi
cd build-mpi
CC=gcc CXX=g++ ../configure \
  --prefix="$HOME/.local/afepack-mpi" \
  --enable-mpi
make -j"$(nproc)"
make test
make install
```

启用 MPI 后，`configure` 会查找 `mpicxx`、`mpic++` 或 `mpiCC`，并验证 MPI C++ 编译器可以编译 C++20 程序。如果检查失败，脚本会提示安装：

```sh
sudo apt install libopenmpi-dev openmpi-bin
```

如果系统上有多个 MPI 实现，建议显式指定：

```sh
MPICXX=/path/to/mpicxx CC=gcc CXX=g++ ../configure --enable-mpi
```

## 7. 测试说明

运行：

```sh
make test
```

当前测试覆盖：

- `examples/poisson_equation`
- `examples/elliptic_equation`
- `examples/TaiChiYinYang`
- `examples/poisson_equation_3D`
- `examples/DysonSphere`

2D 测试会使用源码树内的 `external/easymesh/easymesh` 生成网格。

3D 测试会使用：

- `gmsh`
- `tools/gmsh2mesh`
- `tools/mesh_refine`

如果 `gmsh` 不存在且使用默认 `--enable-gmsh-tests=auto`，3D 测试会跳过并返回 Automake 的 `SKIP` 状态。如果需要 CI 环境强制验证 3D 测试，请配置：

```sh
../configure --enable-gmsh-tests=yes
```

## 8. 安装后的使用

安装后通常会得到：

- `libAFEPack` 库文件
- `AFEPack` 头文件
- 有限元模板文件和模板动态库
- `tools` 下的辅助工具
- `external/easymesh/easymesh`
- `examples` 示例和示例 Makefile

运行安装后的示例时，通常需要设置：

```sh
export LD_LIBRARY_PATH="$HOME/.local/afepack/lib:${LD_LIBRARY_PATH}"
export AFEPACK_PATH="$HOME/.local/afepack/include/AFEPack"
export AFEPACK_TEMPLATE_PATH="$AFEPACK_PATH/template/triangle"
```

3D/tetrahedron 示例可使用：

```sh
export AFEPACK_TEMPLATE_PATH="$AFEPACK_PATH/template/tetrahedron"
```

示例目录中生成的 `Makefile.sample` 和 `run.sh.sample` 可以作为安装后编译、运行用户程序的参考。

## 9. 常见诊断和处理

### 缺少 C++20 编译器

现象：

```text
Missing required dependency: a C++ compiler with C++20 support
C++20 support is required
```

处理：

```sh
sudo apt install build-essential
CC=gcc CXX=g++ ../configure
```

如果系统默认 `g++` 版本过旧，请安装更新版本的 GCC，并通过 `CXX=/path/to/g++` 指定。

### 找不到 cblas.h 或 lapacke.h

现象：

```text
Missing required dependency: cblas.h
Missing required dependency: lapacke.h
```

处理：

```sh
sudo apt install libopenblas-dev liblapacke-dev
```

如果依赖在自定义目录，传入：

```sh
--with-cblas-includedir=DIR
--with-lapacke-includedir=DIR
```

### CBLAS/LAPACKE 链接失败

现象：

```text
CBLAS link check failed
LAPACKE/CBLAS link check failed
```

处理：

- 确认头文件和库来自同一套 BLAS/LAPACK 安装。
- 如果使用系统包，重新安装 `libopenblas-dev liblapacke-dev`。
- 如果使用自定义目录，同时指定 include 和 lib 路径：

```sh
--with-cblas-includedir=DIR
--with-cblas-libdir=DIR
--with-lapacke-includedir=DIR
--with-lapacke-libdir=DIR
```

### OpenBLAS 冲突

现象：

```text
Dependency conflict: --with-openblas-libdir was used, but libopenblas ...
does not link both cblas_dgemm and LAPACKE_dgetrf
```

原因：

当前配置方式把 `--with-openblas-libdir` 视为“OpenBLAS 同时提供 CBLAS 和 LAPACKE”的组合模式。如果该 OpenBLAS 构建没有 LAPACKE 符号，链接检查会失败。

处理：

- 使用带 LAPACKE 支持的 OpenBLAS 构建；或
- 不传 `--with-openblas-libdir`，改用匹配的 `--with-cblas-libdir` 和 `--with-lapacke-libdir`；或
- 使用发行版提供的匹配包：

```sh
sudo apt install libopenblas-dev liblapacke-dev
```

### Boost 或 Boost.Serialization 冲突

现象：

```text
Boost headers not found or too old
Boost.Serialization link check failed
```

处理：

```sh
sudo apt install libboost-serialization-dev
```

如果 Boost 在自定义目录，确保头文件和库来自同一前缀：

```sh
--with-boost=/path/to/boost
--with-boost-libdir=/path/to/boost/lib
```

### gmsh 缺失

现象：

```text
gmsh not found; 3D example run tests will be skipped
```

默认 `auto` 模式下这是警告，`make test` 仍可通过。若需要运行 3D 示例测试：

```sh
sudo apt install gmsh
```

如果配置了 `--enable-gmsh-tests=yes`，缺少 `gmsh` 会使 `configure` 失败，这是预期行为。

### MPI 编译器缺失或不兼容

现象：

```text
--enable-mpi requires an MPI C++ compiler
MPI C++20 compile check failed
```

处理：

```sh
sudo apt install libopenmpi-dev openmpi-bin
```

如果系统存在多个 MPI，请指定：

```sh
MPICXX=/path/to/mpicxx ../configure --enable-mpi
```

## 10. 清理

只清理编译产物：

```sh
make clean
```

清理 `configure` 生成的文件，使构建目录回到未配置状态：

```sh
make distclean
```

如果使用源码树外构建目录，也可以直接删除整个构建目录。
