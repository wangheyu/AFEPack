/**
 * @file   ex29.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Dec 28 13:42:45 2009
 * 
 * @brief  和ex30配合，利用去除背景单元进行分区。
 * 
 *         用法： mpirun -np # ./ex29 src M dst N [L]
 *
 *               #: 旧分区个数
 *             src: 源数据文件名(AFEPack内部格式网格数据文件)
 *               M: 全局加密的轮数
 *             dst: 数据存储目录
 *               N: 新分区的个数
 *               L: 移除背景单元轮数(缺省为M)
 * 
 *  一个典型的用法是：对一个 Easymesh 产生的网格，用 #=1 情形运行，获
 *  得一个二进制形式的分区了的数据，这样可以对数据进行初始分区。
 *
 *  我们在背景单元移除以后进行基于希尔伯特空间曲线填充的背景单元排序，
 *  使得分区相对处于比较合理的情形。
 *
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_ULoadBalance.h>

#ifdef DIM3
#define DIM 3
#endif

#ifndef DIM
#define DIM 2
#endif

#define DOW DIM

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  forest_t forest(MPI_COMM_WORLD);

  if (argc < 5) {
    system("head -n 20 ex29.cpp");
    MPI_Finalize();
    return 1;
  }

  u_int global_round = atoi(argv[2]); /// 第二个参数是全局加密的轮数
  int n_new_rank = atoi(argv[4]); /// 第四个参数是新分区的分区个数
  u_int n_erase_round = global_round;
  if (argc >= 6) {
    n_erase_round = atoi(argv[5]); /// 第五个参数是移除背景单元的轮数
  }

  forest.readMesh(argv[1]); /// 第一个参数是网格文件名

  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  ir_mesh.globalRefine(global_round);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  /// 移除背景单元
  ir_mesh.eraseRootElement(n_erase_round);
  forest.eraseRootElement(n_erase_round);
  forest.renumerateRootElement();

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.config(ir_mesh);
  hlb.partition(n_new_rank);
  Migration::initialize(forest.communicator());

  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);
  hlb.save_data(argv[3], bvs); /// 第三个参数是数据存储目录

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */

