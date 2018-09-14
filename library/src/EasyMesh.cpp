/**
 * @file   EasyMesh.cpp
 * @author Robert Lie
 * @date   Fri Jun 23 15:34:07 2003
 * 
 * @brief  
 * 
 * 
 */

#include "EasyMesh.templates.h"

AFEPACK_OPEN_NAMESPACE

#define DOW 2
template class TriangleMesh<DOW>;
#undef DOW

#define DOW 3
template class TriangleMesh<DOW>;
#undef DOW

void EasyMesh::writeOpenDXData(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(12);
  os.setf(std::ios::scientific, std::ios::floatfield);
  u_int n_node = n_point();
	
  os << "object 1 class array type float rank 1 shape 2 item " 
     << n_node << " data follows\n";
  for (u_int i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }
	
  u_int n_element = n_geometry(2);
  os << "\nobject 2 class array type int rank 1 shape 3 item "
     << n_element << " data follows\n";
  for (u_int i = 0;i < n_element;i ++) {
    os << geometry(2,i).vertex(0) << "\t"
       << geometry(2,i).vertex(1) << "\t"
       << geometry(2,i).vertex(2) << "\t\n";
  }
  os << "attribute \"element type\" string \"triangles\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction-2d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "end\n";
  os.close();
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
