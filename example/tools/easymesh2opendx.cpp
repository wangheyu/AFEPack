/**
 * @file   easymesh2opendx.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Jan  1 09:38:31 2009
 * 
 * @brief  将一个 EasyMesh 网格数据文件转换为可以用 Open DX 显示的格式
 * 
 * 
 */

#include <iostream>
#include <AFEPack/Geometry.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>

int main(int argc, char * argv[])
{
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " easymesh_filename opendx_filename"
              << std::endl;
    return 1;
  }
  EasyMesh mesh;
  mesh.readData(argv[1]);
  mesh.writeOpenDXData(argv[2]);
  return 0;
}

/**
 * end of file
 * 
 */


