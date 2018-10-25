# AFEPack
（Adaptive Finite Element Package）

按照作者李若教授的指示上传一个AFEPack的版本在Github。目前首先维持它在Ubuntu 18.04+deal.II-8.1.0平台上稳定运行。然后逐步将边界条件处理用iterator操作代替指针操作，使得8.1.0以后的deal.II版本可以继续支持。欢迎有能力的同学参与维护和修改。

总结一下目前在Ubuntu 18.04上的安装方案， dealii采用8.1.0， 并行库依赖openmpi和trilinos. 首先建议将源改成aliyun. 请注意这里的安装方案基于这个源提供的AFEPack包，而不是李若老师的CVS-snapshot。二者的主要区别在安装配置文件，也就是configure.in中。接下去如果有我处理不了的bug，我也会忽悠李老师亲自来处理。

安装如下包：

sudo apt-get install cmake

sudo apt-get install g++

sudo apt-get install liblapack-dev

sudo apt-get install libboost-all-dev

sudo apt-get install libmumps-dev

sudo apt-get install libtbb-dev

sudo apt-get install trilinos-all-dev

sudo apt-get install libsuitesparse-dev

sudo apt-get install automake

做如下链接修改， 这里不确定这样做是否最好， 但至少可以继续下去。

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

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DDEAL_II_WITH_PETSC=OFF ..

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

./configure

如果configure出问题， 自己再修改configure.in文件， 主要要检查一下deal.II库和头文件的位置， 以及c++的std标准。

make

中间会出错， 增加链接：

sudo ln -s /usr/local/AFEPack/library/include/ /usr/local/include/AFEPack

sudo ln -s /usr/local/AFEPack/library/lib/*.so /usr/local/lib/

继续：

make

结束。

注意： 

1. 建议在/etc/ld.so.conf.d下增加一个.conf文件， 里面内容为：

/usr/local/lib

/usr/local/AFEPack/library

/usr/local/dealii-8.1.0/lib

或相应的库的位置。 然后

sudo ldconfig

2. 目前AFEPack中的Boundary部分是不能使用的，所以example下的Poisson例子已经不能运行。 我增加了一个withoutBilinearOperator的例子， 只要手工设置边界条件，暂时都不会有问题。 之前包的设置和安装的一主要目的是能使用Trilinos做并行，同时可以用它的AMGPreconditioner做二次元的AMG预处理（李老师的AMGSolver暂时只能处理一次元， 但效率很高）。见examples下增加了一个Poisson_withTrilinos的例子。

3. 此方案未经严格测试， 大家小心使用， 有问题及时反馈给我，谢谢！

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
