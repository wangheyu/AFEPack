/**
 * @file   ex14.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Wed Oct 21 20:19:48 2009
 * 
 * @brief  测试动态负载平衡：载入并输出。
 * 
 * 
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#define DIM 2
#define DOW DIM

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  Migration::initialize(forest.communicator());
  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.load_data("tmp", bvs);

  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  hlb.config(ir_mesh);

  int n_new_rank = atoi(argv[1]);
  hlb.partition(n_new_rank);
  hlb.save_data(argv[2], bvs);

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

