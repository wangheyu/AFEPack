/**
 * @file   mesh2opendx.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Jan  1 09:42:35 2009
 * 
 * @brief  将内部网格数据文件格式转化为 Open DX 的数据格式用于显示
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

int main(int argc, char * argv[])
{
  if (argc != 4) {
    std::cout << "Usage: " << argv[0]
              << " dimension input_mesh output_mesh" << std::endl;
    return 1;
  }
  int dimension = atoi(argv[1]);
  if (dimension == 2) {
#define DIM 2
    HGeometryTree<DIM> h_geometry_tree;
    h_geometry_tree.readMesh(argv[2]);
    IrregularMesh<DIM> irregular_mesh(h_geometry_tree);
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
    regular_mesh.writeOpenDXData(argv[3]);
#undef DIM
  } else if (dimension == 3) {
#define DIM 3
    HGeometryTree<DIM> h_geometry_tree;
    h_geometry_tree.readMesh(argv[2]);
    IrregularMesh<DIM> irregular_mesh(h_geometry_tree);
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
    regular_mesh.writeOpenDXData(argv[3]);
#undef DIM
  }
  return 0;
}

/**
 * end of file
 * 
 */
