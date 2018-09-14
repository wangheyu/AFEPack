/**
 * @file   ex44.cpp
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Mon Dec  6 15:24:21 2010
 * 
 * @brief  测试修正根单元的动态负载平衡：载入并输出。
 * 
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
  Migration::initialize(forest.communicator());
  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.load_data(argv[1], bvs); /// 第一个参数是旧数据目录

  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  hlb.config(ir_mesh);

  int n_new_rank = atoi(argv[2]); /// 第二个参数是新分区个数
  hlb.partition(n_new_rank);
  hlb.save_data(argv[3], bvs); /// 第三个参数是新数据目录

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

