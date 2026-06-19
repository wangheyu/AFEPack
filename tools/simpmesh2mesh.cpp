/**
 * @file   simpmesh2mesh.cpp
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Wed Aug 31 10:42:33 2011
 * 
 * @brief  将三维四面体网格从SimplestMesh格式转换为Mesh格式。
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>

#define DIM 3

int main(int argc, char * argv[]) {
  std::vector<TemplateGeometry<DIM> > tg(1);
  tg[0].readData("tetrahedron.tmp_geo");

  SimplestMesh<DIM> simp_mesh;
  simp_mesh.setTemplateGeometry(tg);

  std::ifstream nis(argv[1]);
  u_int n_vtx;
  nis >> n_vtx;
  simp_mesh.point().resize(n_vtx);
  for (u_int i = 0;i < n_vtx;++ i) {
    int vtx_idx, vtx_bm;
    nis >> vtx_idx;
    nis >> simp_mesh.point(vtx_idx - 1);
    nis >> vtx_bm;
  }
  nis.close();

  std::ifstream eis(argv[2]);
  u_int n_ele;
  eis >> n_ele;
  simp_mesh.element().resize(n_ele);
  for (u_int i = 0;i < n_ele;++ i) {
    int ele_idx;
    eis >> ele_idx;
    simp_mesh.element(i).template_element = 0;
    std::vector<int>& ele_vtx = simp_mesh.element(i).vertex;
    ele_vtx.resize(4);
    for (u_int j = 0;j < 4;++ j) {
      eis >> ele_vtx[j];
      ele_vtx[j] -= 1;
    }
  }
  eis.close();

  Mesh<DIM> mesh;
  simp_mesh.generateMesh(mesh);
  mesh.writeData(argv[3]);

  return 0;
}

/**
 * end of file
 * 
 */
