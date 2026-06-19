/**
 * @file   dbmesh2mesh.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Jan  1 10:17:37 2009
 * 
 * @brief  将 DBMesh 数据文件转换为内部数据文件格式
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/DBMesh.h>

int main(int argc, char * argv[])
{
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " dbmesh-file mesh-file"
              << std::endl;
    exit(1);
  }
  DBMesh dbmesh;
  dbmesh.readData(argv[1]);
  Mesh<2,2> mesh;
  dbmesh.generateMesh(mesh);
  mesh.writeData(argv[2]);
  return 0;
}

/**
 * end of file
 * 
 */


