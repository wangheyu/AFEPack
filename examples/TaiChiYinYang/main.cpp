////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <stdio.h>
#include <dlfcn.h>

#include <iostream>
#include <fstream>

// #include <base/exceptions.h>
// #include <lac/vector.h>
// #include <lac/sparsity_pattern.h>
// #include <lac/sparse_matrix.h>

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

#define DIM 2
const Point<DIM> center (0.0, 0.0);
const double r = 1.;

int main(int argc, char * argv[])
{
  HGeometryTree<DIM> h_tree;
  h_tree.readEasyMesh(argv[1]);
  IrregularMesh<DIM> irregular_mesh(h_tree);
  int refine_time = atoi(argv[2]);
  irregular_mesh.globalRefine(refine_time);

  do {
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
    regular_mesh.writeOpenDXData("D.dx");
    std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
    getchar();

    Indicator<DIM> indicator(regular_mesh);


    double r_half = 0.5 * r;
    Point<DIM> center_upper (0.0, r_half);
    Point<DIM> center_lower (0.0, -r_half);
    double r_half3 = 0.125 * r;
    
    for (int i = 0;i < regular_mesh.n_geometry(2);i ++) {
      Point<DIM>& p0 = regular_mesh.point(regular_mesh.geometry(2,i).vertex(0));
      Point<DIM>& p1 = regular_mesh.point(regular_mesh.geometry(2,i).vertex(1));
      Point<DIM>& p2 = regular_mesh.point(regular_mesh.geometry(2,i).vertex(2));
      Point<DIM> p((p0[0] + p1[0] + p2[0])/3., (p0[1] + p1[1] + p2[1])/3.);
      double area = (p1[0] - p0[0])*(p2[1] - p0[1])
	- (p2[0] - p0[0])*(p1[1] - p0[1]);

      if((distance(p, center) < r
	  &&
	  ((distance(p, center_lower) < r_half3)
	   ||
	   (p[1] > 0. && p[0] < 0. && distance(p, center_upper) > r_half3)
	   ||
	   (p[0] > 0. && distance(p, center_upper) > r_half3 && distance(p, center_upper) < r_half)
	   ||
	   (p[1] < 0. && p[0] < 0. && distance(p, center_lower) > r_half)))
	 ||
	 fabs(distance(p, center) - r) < 0.001 ){
	indicator[i] = area;	
      }
    };

    MeshAdaptor<DIM> mesh_adaptor(irregular_mesh);
    mesh_adaptor.convergenceOrder() = 0.;
    mesh_adaptor.refineStep() = 1;
    mesh_adaptor.setIndicator(indicator);
    mesh_adaptor.tolerence() = 1.0e-06;
    mesh_adaptor.adapt();
  } while (1);	
};

#undef DIM

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
