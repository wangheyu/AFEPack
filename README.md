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
