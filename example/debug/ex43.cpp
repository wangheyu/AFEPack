/**
 * @file   ex43.cpp
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Sun Dec  5 21:37:37 2010
 * 
 * @brief  测试修正根单元动态负载平衡：载入。
 * 
 * 一个附带的作用是将一个二进制输出的数据目录转化成为 dx 文件用于可视
 * 化。
 * 
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_ULoadBalance.h>

#define DIM 2
#define DOW DIM

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.load_data(argv[1], bvs);

  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  char filename[256];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  sprintf(filename, "mesh%d.dx", rank);
  ir_mesh.regularMesh().writeOpenDXData(filename);
  MPI_Finalize();

  return 0;
}

/**
 * end of file
 * 
 */

