# AFEPack
（Adaptive Finite Element Package）

自适应有限元计算软件包，是北京大学李若教授开发和长期维护，面向数值计算科研人员的计算软件包。它以 C++ 编写，要求使用者具有 Linux 平台下 C++ 编程能力。目前我们在李老师指导下，将主要维护工作转移到 GitHub 这里。欢迎有能力的朋友参与维护和修改（想成为合作者直接联系我）。

## 编译安装

### 依赖

AFEPack 的编译依赖如下：

- C/C++ 编译器（gcc/g++）
- MPI（OpenMPI 或 MPICH，含 mpicc/mpicxx）
- [Boost](https://www.boost.org/)（需包含 serialization 头文件）
- [OpenBLAS](https://www.openblas.net/)（提供 BLAS 和 LAPACK 接口）
- GNU Autotools（autoconf、automake、libtool）

在 Ubuntu / Debian 上安装依赖：

```bash
sudo apt install build-essential automake autoconf libtool
sudo apt install libopenblas-dev libboost-all-dev mpich
```

### 编译

```bash
cd AFEPack
autoreconf -fi
./configure
make -j$(nproc)
```

`./configure` 默认使用 `mpicxx` / `mpicc` 作为编译器。如果需要指定 MPI 路径或其他编译选项：

```bash
./configure --prefix=/path/to/install \
            CC=/usr/bin/mpicc \
            CXX=/usr/bin/mpicxx
```

### 安装

```bash
make install
```

`make install` 会将头文件（`include/AFEPack/`）、库文件（`libAFEPack.so`、`libAFEPack_mpi.so`）、模板数据以及 easymesh 工具安装到目标目录（默认 `/usr/local`）。

安装后，确保 `lib` 目录在 `LD_LIBRARY_PATH` 中：

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

### 编译运行例子

例子位于 `examples/` 目录下，每个子目录包含源码和 `Makefile.sample`。安装 AFEPack 后：

```bash
cd examples/poisson_equation
cp Makefile.sample Makefile
make
./main D.d
```

其余例子（`examples/step-7/`、`examples/moving_mesh/` 等）用法类似。

> 注：`examples/` 下部分额外例子（`debug/`、`FDM/`、`Poisson_withTrilinos/`、`template_element/` 等）仍使用旧版 Makefile，需单独编译。

### 关于 deal.II 和 Trilinos

当前版本的 AFEPack 已将 deal.II 和 Trilinos 依赖移去。如果你需要 deal.II 接口或 Trilinos 线性求解器支持，请使用旧版 AFEPack 或参考 Git 历史中 `configure.ac` 里被注释的部分自行添加。
