/**
 * @file   ex41.cpp
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Mon May 31 12:06:29 2010
 * 
 * @brief  此文件和 ex29 完全一致，只是用来检测重排背景单元时候使用了
 *         一个坐标变换。
 * 
 *         用法： mpirun -np # ./ex41 src M N
 *
 *               #: 旧分区个数
 *             src: 源数据文件名(AFEPack内部格式网格数据文件)
 *               M: 全局加密的轮数(缺省为 1)
 *               N: 新分区的个数(缺省为旧分区个数)
 * 
 * 此处使用的变换是二维的极坐标变换，考察在圆环区域下的效果。
 *
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#define DIM 2
#define DOW 2

#define EPS (1.0e-12)

void map(const double * xyz, double * rtheta) {
  const double& x = xyz[0];
  const double& y = xyz[1];
  const double& z = xyz[2];

  double& r = rtheta[0];
  double& theta = rtheta[1];
  double& null = rtheta[2];

  r = sqrt(x*x + y*y);
  if (x > EPS) {
    theta = atan(y/x);
  } else if (x < -EPS) {
    theta = atan(y/x) + M_PIl;
  } else {
    theta = ((y > 0)?(0.5):(-0.5))*M_PIl;
  }
  null = 0.0;
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  forest_t forest(MPI_COMM_WORLD);

  if (argc < 3) {
    system("head -n 20 ex29.cpp");
    MPI_Finalize();
    return 1;
  }

  u_int global_round = 1;
  if (argc >= 3) {
    global_round = atoi(argv[2]); /// 第二个参数是全局加密的轮数
  }
  int n_new_rank = forest.n_rank();
  if (argc >= 4) {
    n_new_rank = atoi(argv[3]); /// 第三个参数是新分区的分区个数
  }

  forest.readMesh(argv[1]); /// 第一个参数是网格文件名

  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  ir_mesh.globalRefine(global_round);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  /// 移除背景单元
  ir_mesh.eraseRootElement(global_round);
  forest.eraseRootElement(global_round);
  forest.renumerateRootElement(&map);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.config(ir_mesh);
  hlb.partition(n_new_rank);
  Migration::initialize(forest.communicator());

  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);
  hlb.save_data("tmp", bvs);

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */

