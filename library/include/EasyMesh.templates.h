/**
 * @file   EasyMesh.templates.h
 * @author Robert Lie
 * @date   Fri Jun 23 15:58:51 2003
 * 
 * @brief  
 * 
 * 
 */

#include "EasyMesh.h"

AFEPACK_OPEN_NAMESPACE

template <int DOW>
void TriangleMesh<DOW>::readData(const std::string& filename)
{
  int i, j;
  int n_node, n_side, n_element;
  char text[64];

  std::cout << "Reading easymesh data file ..." << std::endl;	
  std::cout << "\treading node data ..." << std::flush;
  std::ifstream is((filename + ".n").c_str());
  is >> n_node >> n_element >> n_side;
  is.getline(text, 64);
  this->Mesh<2,DOW>::point().resize(n_node);
  this->Mesh<2,DOW>::geometry(0).resize(n_node);
  for (i = 0;i < n_node;i ++) {
    is >> j
       >> this->Mesh<2,DOW>::point(i)
       >> this->Mesh<2,DOW>::boundaryMark(0,i);
    this->Mesh<2,DOW>::geometry(0,i).index() = j;
    this->Mesh<2,DOW>::geometry(0,i).vertex().resize(1,j);
    this->Mesh<2,DOW>::geometry(0,i).boundary().resize(1,j);
  }
  is.close();
  std::cout << " OK!" << std::endl;
	
  std::cout << "\treading side data ..." << std::flush;
  is.open((filename + ".s").c_str());
  is >> i;
  Assert(i == n_side, ExcMeshData("in side file: side number error"));
  this->Mesh<2,DOW>::geometry(1).resize(n_side);
  for (i = 0;i < n_side;i ++) {
    Geometry& g = this->Mesh<2,DOW>::geometry(1,i);
    g.vertex().resize(2);
    is >> g.index()
       >> g.vertex(0) >> g.vertex(1)
       >> j >> j
       >> this->Mesh<2,DOW>::boundaryMark(1,i);
    this->Mesh<2,DOW>::geometry(1,i).boundary() = this->Mesh<2,DOW>::geometry(1,i).vertex();
  }
  is.close();
  std::cout << " OK!" << std::endl;
	
  std::cout << "\treading element data ..." << std::flush;
  is.open((filename + ".e").c_str());
  is >> i;
  Assert(i == n_element, ExcMeshData("in element file: element number error"));
  is >> i;
  Assert(i == n_node, ExcMeshData("in element file: node number error"));
  is >> i;
  Assert(i == n_side, ExcMeshData("in element file: side number error"));
  is.getline(text, 64);
  this->Mesh<2,DOW>::geometry(2).resize(n_element);
  for (i = 0;i < n_element;i ++) {
    Geometry& g = this->Mesh<2,DOW>::geometry(2,i);
    g.vertex().resize(3);
    g.boundary().resize(3);
    is >> g.index()
       >> g.vertex(0) >> g.vertex(1) >> g.vertex(2)
       >> j >> j >> j
       >> g.boundary(0) >> g.boundary(1) >> g.boundary(2);
    this->Mesh<2,DOW>::boundaryMark(2,i) = 0;
#if (DOW==2)
    const Point<DOW>& p0 = this->Mesh<2,DOW>::point(g.vertex(0));
    const Point<DOW>& p1 = this->Mesh<2,DOW>::point(g.vertex(1));
    const Point<DOW>& p2 = this->Mesh<2,DOW>::point(g.vertex(2));
    if ((p1[0] - p0[0])*(p2[1] - p0[1]) - (p1[1] - p0[1])*(p2[0] - p0[0]) < 0) {
      std::cerr << "vertices of element " << i 
		<< " is not correctly ordered." << std::endl;
      j = g.vertex(2);
      g.vertex(2) = g.vertex(1);
      g.vertex(1) = j;
      j = g.boundary(2);
      g.boundary(2) = g.boundary(1);
      g.boundary(1) = j;
    }
#endif
  }
  is.close();
  std::cout << " OK!" << std::endl;
}


template <int DOW>
void TriangleMesh<DOW>::writeData(const std::string& filename) const
{
  int i, j, k;
	
  std::cout << "Writing easymesh data file ..." << std::endl;	
  std::cout << "\tpreparing data ..." << std::flush;
  int n_node = this->Mesh<2,DOW>::n_geometry(0);
  int n_side = this->Mesh<2,DOW>::n_geometry(1);
  int n_element = this->Mesh<2,DOW>::n_geometry(2);
  std::vector<std::vector<int> > side_neighbour(2,
						std::vector<int>(n_side, -1));
  std::vector<std::vector<int> > element_neighbour(3,
 						   std::vector<int>(n_element, -1));

  for (i = 0;i < n_element;i ++) {
    const GeometryBM& the_element = this->Mesh<2,DOW>::geometry(2, i);
    for (j = 0;j < 3;j ++) {
      k = the_element.boundary(j);
      const GeometryBM& the_side = this->Mesh<2,DOW>::geometry(1, k);
      if (the_side.vertex(0) == the_element.vertex((j+1)%3)) {
	Assert(the_side.vertex(1) == the_element.vertex((j+2)%3), ExcInternalError());
	side_neighbour[0][k] = i;
      }
      else if (the_side.vertex(0) == the_element.vertex((j+2)%3)) {
	Assert(the_side.vertex(1) == the_element.vertex((j+1)%3), ExcInternalError());
	side_neighbour[1][k] = i;
      }
    }
  }
  for (i = 0;i < n_element;i ++) {
    const GeometryBM& the_element = this->Mesh<2,DOW>::geometry(2, i);
    for (j = 0;j < 3;j ++) {
      k = the_element.boundary(j);
      const GeometryBM& the_side = this->Mesh<2,DOW>::geometry(1, k);
      if (the_side.vertex(0) == the_element.vertex((j+1)%3)) {
	element_neighbour[j][i] = side_neighbour[1][k];
      }
      else if (the_side.vertex(0) == the_element.vertex((j+2)%3)) {
	element_neighbour[j][i] = side_neighbour[0][k];
      }
      else {
	Assert(false, ExcInternalError()); // something must be wrong!
      }
    }
  }
  std::cout << " OK!" << std::endl;

  std::cout << "\twriting node data ..." << std::flush;
  std::ofstream os((filename+".n").c_str());
  os.precision(12);
  os.setf(std::ios::scientific, std::ios::floatfield);
  os << n_node << "\t"
     << n_element << "\t"
     << n_side << " **(Nnd, Nee, Nsd)**\n";
  for (i = 0;i < n_node;i ++) {
    j = this->Mesh<2,DOW>::geometry(0, i).vertex(0);
    os << i << "\t"
       << this->Mesh<2,DOW>::point(j) << "\t"
       << this->Mesh<2,DOW>::boundaryMark(0, i) << "\n";
  }
  os << "---------------------------------------------------------------\n"
     << "   n:  x                          y                            \n";
  os.close();
  std::cout << " OK!" << std::endl;
	
  std::cout << "\twriting side data ..." << std::flush;
  os.open((filename+".s").c_str());
  os << n_side << "\n";
  for (i = 0;i < n_side;i ++) {
    const GeometryBM& the_side = this->Mesh<2,DOW>::geometry(1, i);
    os << i << "\t"
       << the_side.vertex(0) << "\t"
       << the_side.vertex(1) << "\t"
       << side_neighbour[0][i] << "\t"
       << side_neighbour[1][i] << "\t"
       << this->Mesh<2,DOW>::boundaryMark(1, i) << "\n";
  }
  os << "---------------------------------------------------------------\n"
     << "   s:    c      d       ea       eb                            \n";
  os.close();
  std::cout << " OK!" << std::endl;

  std::cout << "\twriting element data ..." << std::flush;
  os.open((filename+".e").c_str());
  os << n_element << "\t"
     << n_node << "\t"
     << n_side << " **(Nee, Nnd, Nsd)**\n";
  for (i = 0;i < n_element;i ++) {
    const GeometryBM& the_element = this->Mesh<2,DOW>::geometry(2, i);
    os << i << "\t"
       << the_element.vertex(0) << "\t"
       << the_element.vertex(1) << "\t"
       << the_element.vertex(2) << "\t"
       << element_neighbour[0][i] << "\t"
       << element_neighbour[1][i] << "\t"
       << element_neighbour[2][i] << "\t"
       << the_element.boundary(0) << "\t"
       << the_element.boundary(1) << "\t"
       << the_element.boundary(2) << "\n";
  }
  os << "---------------------------------------------------------------\n"
     << "   e:  i,   j,   k,   ei,   ej,   ek,   si,   sj,   sk         \n";
  os.close();
  std::cout << " OK!" << std::endl;
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

