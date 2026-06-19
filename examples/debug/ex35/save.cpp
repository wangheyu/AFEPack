/**
 * @file   save.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Mar 22 10:33:34 2010
 * 
 * @brief  测试周期区域情形下的自适应和序列化：存储
 * 
 * 程序需要使用一个进程运行，参数表为：
 *
 * 第一个参数：网格文件名；
 * 第二个参数：周期描述文件名；
 * 第三个参数：随机加密轮数；
 * 第四个参数：存储的分区数；
 *
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>
#include <AFEPack/mpi/MPI_PeriodHandler.h>

#define DIM 2
#define DOW 2

typedef AFEPack::MPI::Periodic::PointDistance<DOW> pd_t;
typedef AFEPack::MPI::HGeometryForest<DIM,DOW,pd_t> forest_t;
typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
typedef AFEPack::MPI::HLoadBalance<forest_t> load_balancer_t;
typedef AFEPack::MPI::BirdViewSet<forest_t> birdview_set_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  forest_t forest(MPI_COMM_WORLD);

  u_int round = 1;
  if (argc >= 4) {
    round = atoi(argv[3]); /// 第三个参数是随机加密的轮数
  }
  int n_new_rank = forest.n_rank();
  if (n_new_rank != 1) {
    std::cout << "This program can only run with ONE process!" 
              << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (argc >= 5) {
    n_new_rank = atoi(argv[4]);
  }

  forest.readMesh(argv[1]); /// 第一个参数是网格文件名
  forest.renumerateRootElement();

  AFEPack::MPI::Periodic::periodizeForest(forest, argv[2]); /// 第二个参数是周
                                                   /// 期描述文件

  ir_mesh_t ir_mesh(forest);
  for (u_int i = 0;i < round;++ i) {
    ir_mesh.randomRefine(20.0); /// 加密 5% 的单元
  }
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  load_balancer_t hlb(forest);
  hlb.config(ir_mesh);
  hlb.partition(n_new_rank);
  Migration::initialize(forest.communicator());

  birdview_set_t bvs;
  bvs.add(ir_mesh);
  hlb.save_data("tmp", bvs);

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */

