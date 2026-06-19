/**
 * @file   ex15.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Oct 22 21:03:19 2009
 * 
 * @brief  测试 MPI::FaceData 名空间的功能。
 * 
 * 
 */

#include <AFEPack/mpi/MPI_FaceData.h>

#define DIM 2
#define DOW DIM

int main(int argc, char * argv[])
{
  typedef AFEPack::MPI::HGeometryForest<DIM,DOW> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef RegularMesh<DIM,DOW> mesh_t;
  typedef AFEPack::MPI::FaceData::Mesh<double,mesh_t> data_packer_t;
  typedef data_packer_t::packer_t packer_t;
  typedef AFEPack::MPI::FaceData::Syncer<forest_t,packer_t> syncer_t;

  MPI_Init(&argc, &argv);

  forest_t forest(MPI_COMM_WORLD);
  forest.readMesh(argv[1]);
  ir_mesh_t ir_mesh(forest);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);

  syncer_t syncer(forest, data_packer_t::get_packer(ir_mesh));

  RegularMesh<DIM,DOW>& mesh = ir_mesh.regularMesh();
  u_int n_ele = mesh.n_geometry(DIM);
  for (u_int i = 0;i < n_ele;++ i) {
    GeometryBM& ele = mesh.geometry(DIM, i);
    u_int n_bnd = ele.n_boundary();
    for (u_int j = 0;j < n_bnd;++ j) {
      u_int bnd_idx = ele.boundary(j);
      GeometryBM& bnd = mesh.geometry(DIM-1, bnd_idx);
      double * p_data = syncer.get_send_buffer(bnd);
      (*p_data) = forest.rank();
    }
  }

  syncer.sync();

  MPI_Finalize();

  return 0;
}

/**
 * end of file
 * 
 */

