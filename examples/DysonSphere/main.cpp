////////////////////////////////////////////////////////////////////////////////////////////
// main4.cpp :
//

#include <stdio.h>
#include <dlfcn.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>

#include <AFEPack/Vector.h>
#include <AFEPack/SparsityPattern.h>
#include <AFEPack/SparseMatrix.h>

#include <AFEPack/Geometry.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>

#define DIM 3
const Point<DIM> center (0.0, 0.0, 0.0);

double tetrahedron_volume(const double * v0, 
			  const double * v1, 
			  const double * v2, 
			  const double * v3)
{
  return ((v1[0] - v0[0])*(v2[1] - v0[1])*(v3[2] - v0[2])
	  + (v1[1] - v0[1])*(v2[2] - v0[2])*(v3[0] - v0[0])
	  + (v1[2] - v0[2])*(v2[0] - v0[0])*(v3[1] - v0[1])
	  - (v1[0] - v0[0])*(v2[2] - v0[2])*(v3[1] - v0[1])
	  - (v1[1] - v0[1])*(v2[0] - v0[0])*(v3[2] - v0[2])
	  - (v1[2] - v0[2])*(v2[1] - v0[1])*(v3[0] - v0[0]))/6.;
};



int main(int argc, char * argv[])
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " mesh refine_time [--test-once]" << std::endl;
    return 1;
  }
  bool test_once = (argc > 3 && strcmp(argv[3], "--test-once") == 0);

  HGeometryTree<DIM> h_tree;
  h_tree.readMesh(argv[1]);
  IrregularMesh<DIM> irregular_mesh(h_tree);
  int refine_time = atoi(argv[2]);
  irregular_mesh.globalRefine(refine_time);

  double l = 0.1;
  do {
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
    regular_mesh.writeOpenDXData("Cube.dx");
    if (!test_once) {
      std::cout << "\nPress ENTER to continue or CTRL+C to stop ..." << std::flush;
      getchar();
    }

    Indicator<DIM> indicator(regular_mesh);
    for (int i = 0;i < regular_mesh.n_geometry(DIM);i ++) {
      int n_vertex = regular_mesh.geometry(DIM,i).n_vertex();
      int vertex[4];
      switch(n_vertex) {
      case 4:
	vertex[0] = 0;
	vertex[1] = 1;
	vertex[2] = 2;
	vertex[3] = 3;
	break;
      case 5:
	vertex[0] = 0;
	vertex[1] = 1;
	vertex[2] = 3;
	vertex[3] = 4;
	break;
      case 7:
	vertex[0] = 0;
	vertex[1] = 1;
	vertex[2] = 2;
	vertex[3] = 3;
	break;
      default:
	assert(false);
      }
      Point<DIM> p((regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[0]))[0] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[1]))[0] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[2]))[0] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[3]))[0])/4.,
		   (regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[0]))[1] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[1]))[1] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[2]))[1] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[3]))[1])/4.,
		   (regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[0]))[2] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[1]))[2] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[2]))[2] +
		    regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[3]))[2])/4.);      
      
      double area = tetrahedron_volume(
				       regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[0])),
				       regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[1])),
				       regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[2])),
				       regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[3])));

      if(distance(p, center) < 0.05
	 ||
	 (distance(p, center) > 0.95 && distance(p, center) < 1.05)){
	indicator[i] = area;
      }
    };

    MeshAdaptor<DIM> mesh_adaptor(irregular_mesh);
    mesh_adaptor.convergenceOrder() = 0.;
    mesh_adaptor.refineStep() = 1;
    mesh_adaptor.setIndicator(indicator);
    mesh_adaptor.tolerence() = 1.0e-08;
    mesh_adaptor.adapt();
    if (test_once) break;
  } while (1);	
};

#undef DIM

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
