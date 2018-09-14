/**
 * @file   ex07.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Sun Oct  4 14:19:45 2009
 * 
 * @brief 测试并行版本下的 match_geometry, sync_is_used, regularize 函
 * 数。
 * 
 * 
 */

#include <AFEPack/mpi/MPI_HGeometry.h>

#define DIM 2
#define DOW 2

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  forest.readMesh(argv[1]); /// 第一个参数是网格文件名

  IrregularMesh<DIM,DOW> ir_mesh(forest);
  u_int round = atoi(argv[2]); /// 第二个参数是随机加密的轮数
  for (u_int i = 0;i < round;++ i) {
    ir_mesh.randomRefine(5.0); /// 加密 5% 的单元
  }

  while (1) {
    ir_mesh.semiregularize();

    bool is_operated = false;
    AFEPack::MPI::HGeometryMatcher<forest_t> matcher(forest);
    is_operated |= matcher.match_geometry();
    is_operated |= matcher.sync_is_used();

    if (! is_operated) break;
  }
  ir_mesh.regularize();

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

