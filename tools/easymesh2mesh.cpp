/**
 * @file   easymesh2mesh.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Jan  1 10:28:19 2009
 * 
 * @brief  
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/EasyMesh.h>

#ifndef D
#define D 2
#endif

int main(int argc, char * argv[])
{
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " input_mesh output_mesh"
              << std::endl;
    return 1;
  }
  TriangleMesh<D> mesh;
  mesh.readData(argv[1]);
  dynamic_cast<Mesh<2,D> *>(&mesh)->writeData(argv[2]);
  return 0;
}

/**
 * end of file
 * 
 */
