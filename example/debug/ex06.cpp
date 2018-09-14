/**
 * @file   ex06.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Sat Oct  3 20:39:17 2009
 * 
 * @brief  测试 HGeometryForest 的 readMesh 函数。
 * 
 * 
 */

#include <mpi.h>
#include <AFEPack/mpi/MPI_HGeometry.h>

#define DIM 2
#define DOW 2

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  AFEPack::MPI::HGeometryForest<DIM,DOW> forest(MPI_COMM_WORLD);
  forest.readMesh(argv[1]);

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */

