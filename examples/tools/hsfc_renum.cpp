/**
 * @file   hsfc_renum.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Mon Oct 19 12:55:40 2009
 * 
 * @brief  
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

int main(int argc, char * argv[])
{
  if (argc != 4) {
    std::cout << "Usage: " << argv[0]
              << " dimension input_mesh output_mesh"
              << std::endl;
  return 1;
  }
  int dimension = atoi(argv[1]);
  if (dimension == 2) {
#define DIM 2
    Mesh<DIM,DIM> mesh;
    mesh.readData(argv[2]);
    mesh.renumerateElementHSFC();
    mesh.writeData(argv[3]);
#undef DIM
  } else if (dimension == 3) {
#define DIM 3
    Mesh<DIM,DIM> mesh;
    mesh.readData(argv[2]);
    mesh.renumerateElementHSFC();
    mesh.writeData(argv[3]);
#undef DIM
  }
  return 0;
}

/*
 * end of file
 * ***************************************************************************/
