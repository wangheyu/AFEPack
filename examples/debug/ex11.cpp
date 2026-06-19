/**
 * @file   ex11.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Wed Oct 14 10:01:34 2009
 * 
 * @brief  测试动态负载平衡：载入后数据的正确性。
 * 
 * 
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#define DIM 2
#define DOW 2

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.load_data("tmp", bvs);

  u_int round = (argc>=2)?atoi(argv[1]):1;
  for (u_int i = 0;i < round;++ i) {
    ir_mesh.randomRefine(20.0); /// 加密 20% 的单元
  }
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

