/**
 * @file   ex13.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct 20 11:33:11 2009
 * 
 * @brief  测试并行的有限元函数的序列化：载入部分
 * 
 * 
 */

#include <AFEPack/mpi/MPI_HGeometry.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>

#define DIM 2
#define DOW 2

typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;

int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  AFEPack::MPI::BirdView<forest_t> ir_mesh(forest);
  AFEPack::MPI::BirdViewSet<forest_t> bvs;
  bvs.add(ir_mesh);

  AFEPack::MPI::HLoadBalance<forest_t> hlb(forest);
  Migration::initialize(forest.communicator());
  hlb.load_data("tmp", bvs);

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

  FEMFunction<double,DIM,DOW> u_h(fem_space);
  Migration::data_id_t id = Migration::name_to_id("u_h");
  Migration::import_fe_func(u_h, id);
  Migration::clear_data_buffer(forest);

  char filename[256];
  sprintf(filename, "u_h%d.dx", forest.rank());
  u_h.writeOpenDXData(filename);
  MPI_Finalize();

  return 0;
}

/**
 * end of file
 * 
 */

