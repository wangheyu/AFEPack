# AFEPack
（Adaptive Finite Element Package）

自适应有限元计算软件包，是北京大学李若教授开发和长期维护，面向数值计算科研人员的计算软件包。它以C++编写，要求使用者具有Linux平台下C++编程能力。目前我们在李老师指导下，将主要维护工作转移到Github这里。目前首先维持它在Ubuntu 18.04+deal.II-8.1.0平台上稳定运行。然后逐步将边界条件处理用iterator操作代替指针操作，使得8.1.0以后的deal.II版本可以继续支持。欢迎有能力的朋友参与维护和修改（想成为合作者直接联系我）。

## 缺省安装
李老师修改了缺省安装策略以方便初学者。如果你用的是Ubuntu 16.04，可以直接采用源中的dealii-8.1.0，跳过下面编译dealii的部分。如果是Ubuntu 18.04及以上，那么你需要手工编译dealii-8.1.0，缺省编译也很简单：

安装如下依赖包：

sudo apt-get install cmake

sudo apt-get install g++

sudo apt-get install libboost-all-dev

sudo apt-get install libtbb-dev

sudo apt-get install automake

对deal.II-8.1.0的代码做如下修改：

deal.II/include/deal.II/lac/sparsity_pattern.h, 在最开始增加一行：

`#include <algorithm>`

deal.II/source/base/parameter_handler.cc, line 1278, 如下修改：

return (p.get_optional<std::string>("value")); -> return bool(p.get_optional<std::string>("value"));

然后在deal.II的目录：

mkdir build

cd build

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 ..

可以自己选择安装目的地。

make

sudo make install

然后可以编译AFEPack了， 将AFEPack在自己目录解开， sudo mv到安装目录（以下为/usr/local/AFEPack）：

cd /usr/local/AFEPack

aclocal

autoconf

automake

./configure --with-dealii="/usr/local/dealii-8.1.0"

make

结束。

## 带Trilinos接口的AFEPack

Trilinos是一个并行计算包，它里面包含了大量的线性求解器和预处理器。如果在dealii中链接了trilinos，那么我们就可以在AFEPack中使用这些求解器，从而增加AFEPack的计算能力。这个操作稍微有点复杂，不建议初学者尝试。

安装如下包：

sudo apt-get install cmake

sudo apt-get install g++

sudo apt-get install liblapack-dev

sudo apt-get install libboost-all-dev

sudo apt-get install libmumps-dev

sudo apt-get install libtbb-dev

sudo apt-get install trilinos-all-dev

sudo apt-get install libsuitesparse-dev

sudo apt-get install libarpack2-dev

sudo apt-get install automake

做如下链接修改， 这里不确定这样做是否最好， 但至少可以继续下去。
（或者定义到你自己指定的一个目录下）
sudo ln -s /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so /usr/local/lib/libscalapack.so

sudo ln -s /usr/lib/x86_64-linux-gnu/libptscotch-6.so /usr/lib/x86_64-linux-gnu/libptscotch.so

sudo ln -s /usr/lib/x86_64-linux-gnu/libptscotcherr-6.so /usr/lib/x86_64-linux-gnu/libptscotcherr.so

sudo ln -s /usr/lib/x86_64-linux-gnu/libscotch-6.so /usr/lib/x86_64-linux-gnu/libscotch.so

sudo ln -s /usr/lib/x86_64-linux-gnu/libscotcherr-6.so /usr/lib/x86_64-linux-gnu/libscotcherr.so

对deal.II-8.1.0的代码做如下修改：

deal.II/include/deal.II/lac/sparsity_pattern.h, 在最开始增加一行：

`#include <algorithm>`

deal.II/source/base/parameter_handler.cc, line 1278, 如下修改：

return (p.get_optional<std::string>("value")); -> return bool(p.get_optional<std::string>("value"));

然后在deal.II的目录：

mkdir build

cd build

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DDEAL_II_WITH_PETSC=OFF -DHDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi ..

这里屏蔽了PETSC，因为最新的PETSC和之前变化太大，除非你自己安装一个很低版本的。 反正源里的这个不行。 同理， 如果你有P4EST， 也屏蔽掉。 继续：

make

sudo make install

设置一下deal.ii的库链接：

sudo ln -s /usr/local/dealii-8.1.0/lib/* /usr/local/lib

然后可以编译AFEPack了， 将AFEPack在自己目录解开， sudo mv到安装目录（以下为/usr/local/AFEPack）：

cd /usr/local/AFEPack

aclocal

autoconf

automake

CC="/usr/bin/mpicc" CXX="/usr/bin/mpic++" CFLAGS="-I/usr/include/trilinos" CPPFLAGS="-I/usr/include/trilinos" CXXFLAGS="-I/usr/include/trilinos" ./configure --with-mpi=yes --with-dealii="/usr/local/dealii-8.1.0"

make

结束。

在examples中增加了一个withoutBilinearOperator的例子， 只要手工设置边界条件，暂时都不会有问题。 之前包的设置和安装的一主要目的是能使用Trilinos做并行，同时可以用它的AMGPreconditioner做二次元的AMG预处理（李老师的AMGSolver暂时只能处理一次元， 但效率很高）。见examples下增加了一个Poisson_withTrilinos的例子。

此方案未经严格测试， 大家小心使用， 有问题及时反馈给我，谢谢！

此外， 在

git@github.com:vickor/Navier_Stokes下有一个用Taylor-Hood元计算NS方程的例子。

## Debian / Ubuntu
此方案已经在 Debian 9.5 和 Ubuntu 18.04 测试过。

首先安装如下依赖
```
sudo apt install liblapack-dev cmake gcc gcc-c++ automake autoconf mpich
```

AFEPack 依赖特定版本的 boost 和 deal.II，以下说明如下安装。

### boost-1.50.0
```
wget https://phoenixnap.dl.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
tar -xjf boost_1_50_0.tar.bz2
cd boost_1_50_0 && ./bootstrap.sh --prefix=/path/to/install
# 根据设置的安装目录决定是否需要 sudo
./b2 && sudo ./b2 install
```

### deal.II-8.1.0
配置 deal.II 为启用 MPI，关闭 threads（目前用不到）。
```
tar -xzf dealii-8.1.0.tar.gz
cd deal.II
mkdir build && cd build

# 根据 boost 安装路径修改
export BOOST_DIR=/path/to/boost/install
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install \
      -DDEAL_II_WITH_MPI=on \
      -DDEAL_II_WITH_THREADS=off \
      -DCMAKE_C_COMPILER=/usr/bin/mpicc \
      -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx \
      -DCMAKE_Fortran_COMPILER=/usr/bin/mpifort \
      ..

make -j 4
# 根据设置的安装目录决定是否需要 sudo
sudo make install
```

**若 deal.II 报 boost 的编译错误，请在出错的语句 return 后加上**

```
static_cast<bool>
```

### AFEPack
AFEPack 依赖为 deal.II 和 boost，编译 MPI 版本的 AFEPack 时需要任意一种 MPI
实现。
```
CC=mpicc CXX=mpicxx ./configure --prefix=/path/to/install \
    --with-boost=/path/to/boost/install \
    --with-dealii=/path/to/dealii/install
make
# 根据情况决定是否需要 sudo
sudo make install
```
`make install` 会将头文件，库文件，所有 template 安装到目标文件夹。随 AFEPack
还附赠一个 easymesh。

在 `examples` 文件夹下可以继续编译测试例子，直接使用 `make` 即可 。
```
cd examples && make
```
