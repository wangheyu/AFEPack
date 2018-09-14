/**
 * @file   ex16.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Oct 23 22:30:39 2009
 * 
 * @brief  测试 MPI::FaceData 名空间的功能
 * 
 * 
 */

#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/mpi/MPI_FaceData.h>

#define DIM 2
#define DOW DIM

int main(int argc, char * argv[])
{
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef DGFEMSpace<double,DIM,DOW> fe_space_t;  
  typedef AFEPack::MPI::FaceData::FESpace<double,fe_space_t> data_packer_t;
  typedef data_packer_t::packer_t packer_t;
  typedef AFEPack::MPI::FaceData::Syncer<forest_t,packer_t> syncer_t;

  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  forest.readMesh(argv[1]);
  ir_mesh_t ir_mesh(forest);

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

  TemplateGeometry<DIM-1> interval;
  interval.readData("interval.tmp_geo");
  CoordTransform<DIM-1,DIM> int22d;
  int22d.readData("interval.to2d.crd_trs");
  
  std::vector<TemplateDGElement<DIM-1,DIM> > dg_tmp_ele(1);
  dg_tmp_ele[0].reinit(interval, int22d);

  RegularMesh<DIM,DOW>& mesh = ir_mesh.regularMesh();
  fe_space_t fem_space(mesh, tmp_ele, dg_tmp_ele);

  u_int n_ele = mesh.n_geometry(DIM);
  fem_space.element().resize(n_ele);
  for (int i = 0;i < n_ele;i ++) {
    fem_space.element(i).reinit(fem_space, i, 0);
  }

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  u_int n_dg_ele = mesh.n_geometry(DIM-1);
  fem_space.dgElement().resize(n_dg_ele);
  for (int i = 0;i < n_dg_ele;i ++) {
    fem_space.dgElement(i).reinit(fem_space, i, 0);
  }
  fem_space.buildDGElement();

  syncer_t syncer(forest, data_packer_t::get_packer(fem_space));

  fe_space_t::DGElementIterator
    the_dg_ele = fem_space.beginDGElement(),
    end_dg_ele = fem_space.endDGElement();
  for (;the_dg_ele != end_dg_ele;++ the_dg_ele) {
    double * p_data = syncer.get_send_buffer(*the_dg_ele);
    (*p_data) = forest.rank();
  }
  syncer.sync();

  MPI_Finalize();

  return 0;
}

/**
 * end of file
 * 
 */

