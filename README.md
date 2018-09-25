# AFEPack
（Adaptive Finite Element Package）

按照作者李若教授的指示上传一个AFEPack的版本在Github。目前首先维持它在Ubuntu 18.04+deal.II-8.1.0平台上稳定运行。然后逐步将边界条件处理用iterator操作代替指针操作，使得8.1.0以后的deal.II版本可以继续支持。欢迎有能力的同学参与维护和修改。
    
suggestion: change the repositories to Aliyun.

sudo apt-get install cmake
sudo apt-get install g++
sudo apt-get install liblapack-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libmumps-dev
sudo apt-get install libtbb-dev
sudo apt-get install trilinos-all-dev
sudo apt-get install libsuitesparse-dev
sudo apt-get install automake


sudo ln -s /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so /usr/local/lib/libscalapack.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libptscotch-6.so /usr/lib/x86_64-linux-gnu/libptscotch.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libptscotcherr-6.so /usr/lib/x86_64-linux-gnu/libptscotcherr.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libscotch-6.so /usr/lib/x86_64-linux-gnu/libscotch.so
sudo ln -s /usr/lib/x86_64-linux-gnu/libscotcherr-6.so /usr/lib/x86_64-linux-gnu/libscotcherr.so

deal.II/include/deal.II/lac/sparsity_pattern.h, add 
#include <algorithm>
in the front.

deal.II/source/base/parameter_handler.cc, line 1278, modify
return (p.get_optional<std::string>("value")); -> return bool(p.get_optional<std::string>("value"));

mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/dealii-8.1.0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DDEAL_II_WITH_PETSC=OFF ..

make
make install

cd /usr/local/AFEPack
aclocal
autoconf
automake

sudo ln -s /usr/local/AFEPack/library/include/ /usr/local/include/AFEPack
sudo ln -s /usr/local/AFEPack/library/lib/*.so /usr/local/lib/
