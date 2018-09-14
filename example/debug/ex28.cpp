/**
 * @file   ex28.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Dec  3 14:54:08 2009
 * 
 * @brief  测试去除一层背景单元。
 * 
 * 
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#define DIM 2
#define DOW 2

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  forest_t forest(MPI_COMM_WORLD);

  u_int global_round = 1;
  if (argc >= 3) {
    global_round = atoi(argv[2]); /// 第二个参数是全局加密的轮数
  }
  u_int random_round = 1;
  if (argc >= 4) {
    random_round = atoi(argv[3]); /// 第二个参数是随机加密的轮数
  }
  int n_new_rank = forest.n_rank();
  if (argc >= 5) {
    n_new_rank = atoi(argv[4]); /// 第三个参数是新分区的分区个数
  }

  forest.readMesh(argv[1]); /// 第一个参数是网格文件名
  forest.renumerateRootElement();

  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  ir_mesh.globalRefine(global_round);
  for (u_int i = 0;i < random_round;++ i) {
    ir_mesh.randomRefine(5.0); /// 加密 5% 的单元
  }
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  /// 移除一层背景单元
  ir_mesh.eraseRootElement();
  forest.eraseRootElement();

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.config(ir_mesh);
  hlb.partition(n_new_rank);
  Migration::initialize(forest.communicator());

  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);
  hlb.save_data("tmp", bvs);

  char filename[256];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  sprintf(filename, "old_mesh%d.dx", rank);
  ir_mesh.regularMesh().writeOpenDXData(filename);

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */

