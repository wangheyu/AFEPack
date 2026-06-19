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

安装后，需设置以下环境变量：

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/bin:$PATH
export AFEPACK_TEMPLATE_PATH=/usr/local/include/AFEPack/template/triangle
```

`AFEPACK_TEMPLATE_PATH` 需指向所用几何类型对应的模板子目录（如 `triangle`、`tetrahedron` 等），多个路径可用 `:` 分隔。

### 编译运行例子

安装后的例子位于 `share/doc/AFEPack/example/`。以 poisson_equation 为例：

```bash
cd share/doc/AFEPack/example/poisson_equation
cp Makefile.sample Makefile
make
easymesh D          # 生成网格 D.n, D.s, D.e
./main D            # 传入网格文件前缀（不含扩展名）
```

> **注意**：部分例子（如 `poisson_equation_3D`、`laplacian_evp3D`）无需生成网格，直接 `./main <输入文件>` 即可。参见各子目录下的 `README` 和 `run.sh.sample`。

> 注：`examples/` 下部分额外例子（`debug/`、`FDM/`、`Poisson_withTrilinos/`、`template_element/` 等）仍使用旧版 Makefile，需单独编译。

### 关于 deal.II 和 Trilinos

当前版本的 AFEPack 已将 deal.II 和 Trilinos 依赖移去。如果你需要 deal.II 接口或 Trilinos 线性求解器支持，请使用旧版 AFEPack 或参考 Git 历史中 `configure.ac` 里被注释的部分自行添加。
