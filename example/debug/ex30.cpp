/**
 * @file   ex30.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Dec 28 13:52:10 2009
 * 
 * @brief  和ex29配合，利用去除背景单元进行分区。
 * 
 *         用法： mpirun -np # ./ex30 src dst M N
 *
 *               #: 旧分区个数
 *             src: 源数据目录
 *             dst: 产生的数据存储的目录
 *               M: 全局加密的轮数(缺省为 1)
 *               N: 新分区的个数(缺省为旧分区个数)
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

  if (argc < 4) {
    system("head -n 16 ex30.cpp");
    MPI_Finalize();
    return 1;
  }

  u_int global_round = 1;
  if (argc >= 4) {
    global_round = atoi(argv[3]); /// 第三个参数是全局加密的轮数
  }
  int n_new_rank = forest.n_rank();
  if (argc >= 5) {
    n_new_rank = atoi(argv[4]); /// 第四个参数是新分区的分区个数
  }

  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  AFEPack::MPI::load_mesh(argv[1], forest, ir_mesh);

  ir_mesh.globalRefine(global_round);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  /// 移除背景单元
  ir_mesh.eraseRootElement(global_round);
  forest.eraseRootElement(global_round);
  forest.renumerateRootElement();

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.config(ir_mesh);
  hlb.partition(n_new_rank);
  Migration::initialize(forest.communicator());

  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);
  hlb.save_data(argv[2], bvs);

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */

