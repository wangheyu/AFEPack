/**
 * @file   mesh2easymesh.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Jan  1 09:51:40 2009
 * 
 * @brief  将内部网格数据文件格式转化为 EasyMesh 网格数据文件，输入文
 *         件中只能有三角形。
 * 
 * 
 */


#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/EasyMesh.h>

int main(int argc, char * argv[])
{
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " input_mesh output_mesh"
              << std::endl;
    return 1;
  }
#if 0
  HGeometryTree<2> h_tree;
  h_tree.readMesh(argv[1]);
  IrregularMesh<2> ir_mesh(h_tree);
  ir_mesh.regularize();
  RegularMesh<2>& mesh = ir_mesh.regularMesh();
  mesh.writeEasyMesh(argv[2]);
#else
  Mesh<2> mesh;
  mesh.readData(argv[1]);
  RegularMesh<2> * p_mesh = (RegularMesh<2> *)(&mesh);
  p_mesh->writeEasyMesh(argv[2]);
#endif
  return 0;
}

/**
 * end of file
 * 
 */
