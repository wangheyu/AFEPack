/**
 * @file   ex12.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct 20 10:48:20 2009
 * 
 * @brief  测试并行的有限元函数的序列化：存储部分
 * 
 * 
 */

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <AFEPack/Serialization.h>
#include <AFEPack/Operator.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#define DOW 2
#define DIM DOW

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

double _u_(const double * p)
{
  return sin(12*p[0])*sin(12*p[1]);
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  forest.readMesh(argv[1]); /// 第一个参数是网格文件名
  u_int round = 1;
  if (argc >= 3) {
    round = atoi(argv[2]); /// 第二个参数是随机加密的轮数
  }
  int n_new_rank = forest.n_rank();
  if (argc >= 4) {
    n_new_rank = atoi(argv[3]); /// 第三个参数是新分区的分区个数
  }

  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  for (u_int i = 0;i < round;++ i) {
    ir_mesh.randomRefine(20.0); /// 加密 5% 的单元
  }
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  TemplateGeometry<DIM> tri;
  tri.readData("triangle.tmp_geo");
  CoordTransform<DIM,DIM> tri_ct;
  tri_ct.readData("triangle.crd_trs");
  TemplateDOF<DIM> tri_td(tri);
  tri_td.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> tri_bf(tri_td);
  tri_bf.readData("triangle.1.bas_fun");

  TemplateGeometry<DIM> twin_tri;
  twin_tri.readData("twin_triangle.tmp_geo");
  CoordTransform<DIM,DIM> twin_tri_ct;
  twin_tri_ct.readData("twin_triangle.crd_trs");
  TemplateDOF<DIM> twin_tri_td(twin_tri);
  twin_tri_td.readData("twin_triangle.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> twin_tri_bf(twin_tri_td);
  twin_tri_bf.readData("twin_triangle.1.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > template_element(2);
  template_element[0].reinit(tri, tri_td, tri_ct, tri_bf);
  template_element[1].reinit(twin_tri, twin_tri_td, twin_tri_ct, twin_tri_bf);

  RegularMesh<DIM,DOW>& mesh = ir_mesh.regularMesh();
  FEMSpace<double,DIM,DOW> fem_space(mesh, template_element);
	
  int n_element = mesh.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++) {
    if (mesh.geometry(DIM, i).n_vertex() == 3) {
      fem_space.element(i).reinit(fem_space, i, 0);
    } else {
      fem_space.element(i).reinit(fem_space, i, 1);
    }
  }

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  hlb.config(ir_mesh);

  hlb.partition(n_new_rank);

  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);

  FEMFunction<double,DIM,DOW> u_h(fem_space);
  Operator::L2Interpolate(_u_, u_h);

  Migration::initialize(forest.communicator());
  Migration::clear_data_buffer(forest);
  Migration::data_id_t id = Migration::register_data_name("u_h");
  Migration::export_fe_func(u_h, id);

  hlb.save_data("tmp", bvs);

  MPI_Finalize();
  return 0;
}

/**
 * end of file
 * 
 */
