/**
 * @file   ex42.cpp
 * @author rli <rli@pelz>
 * @date   Sat Dec  4 15:37:56 2010
 * 
 * @brief  测试动态负载平衡中的修正背景单元功能。
 * 
 * 此处载入一个非常稀疏的网格，然后对其中的第一个单元进行反复加密，使得
 * 这个单元内子单元个数非常多，然后进行存储。假设编译得到的可执行文件
 * 名为 ex42，此例使用方法为：
 *
 *   mpirun -np 1 ./ex42 {mesh_file round n_partition data_dir}
 *
 * 必须使用一个进程运行，读入EasyMesh类型的数据文件 mesh_file(不要加后
 * 缀)，然后对第一个单元加密 round (正整数) 轮，最后将数据拆分为
 * n_partition 个分区，存储在 data_dir 目录中。
 *
 * 为了将过去的例子修改成为现在的负载平衡状态，只需要将头文件
 *
 *   <AFEPack/mpi/MPI_LoadBalance.h>
 *
 * 修改为
 *
 *   <AFEPack/mpi/MPI_ULoadBalance.h>
 *
 * 即可。
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_ULoadBalance.h>

#define DIM 2
#define DOW 2

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  forest_t forest(MPI_COMM_WORLD);

  forest.readMesh(argv[1]); /// 第一个参数是网格文件名
  u_int round = atoi(argv[2]); /// 第二个参数是局部加密的轮数
  int n_new_rank = atoi(argv[3]); /// 第三个参数是新分区的分区个数

  forest.renumerateRootElement();

  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  for (int i = 0;i < round;++ i) {
    ir_mesh.semiregularize();
    ir_mesh.regularize(false);

    u_int n_ele = ir_mesh.regularMesh().n_geometry(DIM);
    std::vector<bool> flag(n_ele, false);
    flag[0] = true;
    ir_mesh.regularMesh().refineLabelled(flag);
  }
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.config(ir_mesh);
  hlb.partition(n_new_rank);
  Migration::initialize(forest.communicator());

  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);
  hlb.save_data(argv[4], bvs); /// 第四个参数是数据存储目录

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

