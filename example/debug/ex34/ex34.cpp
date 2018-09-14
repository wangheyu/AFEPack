/**
 * @file   ex34.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Mar 22 10:33:34 2010
 * 
 * @brief  测试周期区域情形下的网格生成和自适应。
 * 
 * 程序需要使用一个进程运行，参数表为：
 *
 * 第一个参数：网格文件名；
 * 第二个参数：周期描述文件名；
 * 第三个参数：随机加密轮数；
 *
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_PeriodHandler.h>

#define DIM 2
#define DOW 2

typedef AFEPack::MPI::Periodic::PointDistance<DOW> pd_t;
typedef AFEPack::MPI::HGeometryForest<DIM,DOW,pd_t> forest_t;
typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  forest_t forest(MPI_COMM_WORLD);

  u_int round = 1;
  if (argc >= 4) {
    round = atoi(argv[3]); /// 第三个参数是随机加密的轮数
  }

  forest.readMesh(argv[1]); /// 第一个参数是网格文件名
  forest.renumerateRootElement();

  /// 第二个参数是周期描述文件
  AFEPack::MPI::Periodic::periodizeForest(forest, argv[2]); 

  ir_mesh_t ir_mesh(forest);
  for (u_int i = 0;i < round;++ i) {
    ir_mesh.randomRefine(20.0); /// 加密 5% 的单元
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

