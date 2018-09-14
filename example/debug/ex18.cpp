/**
 * @file   ex16.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Oct 23 22:30:39 2009
 * 
 * @brief  测试 MPI::DOF 名空间的自由度匹配代码
 * 
 * 
 */

#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/mpi/MPI_DOF.h>

#define DIM 2
#define DOW DIM

int main(int argc, char * argv[])
{
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef FEMSpace<double,DIM,DOW> fe_space_t;  
  typedef AFEPack::MPI::DOF::GlobalIndex<forest_t, fe_space_t> global_index_t;

  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  forest.readMesh(argv[1]);
  ir_mesh_t ir_mesh(forest);

  int round = 0;
  if (argc >= 3) round = atoi(argv[2]);

  ir_mesh.globalRefine(round);
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

  std::vector<TemplateElement<double,DIM,DIM> > tmp_ele(1);
  tmp_ele[0].reinit(tri, tri_td, tri_ct, tri_bf);

  RegularMesh<DIM,DOW>& mesh = ir_mesh.regularMesh();
  fe_space_t fem_space(mesh, tmp_ele);

  u_int n_ele = mesh.n_geometry(DIM);
  fem_space.element().resize(n_ele);
  for (int i = 0;i < n_ele;i ++) {
    fem_space.element(i).reinit(fem_space, i, 0);
  }

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  global_index_t global_index(forest, fem_space);
  global_index.build();

  MPI_Finalize();

  return 0;
}

/**
 * end of file
 * 
 */

