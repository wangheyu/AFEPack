////////////////////////////////////////////////////////////////////////////////////////////
// main4.cpp :
//

#include <stdio.h>
#include <dlfcn.h>

#include <iostream>
#include <fstream>

#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include <AFEPack/Geometry.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>

#define DIM 3

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
	HGeometryTree<DIM> h_tree;
	h_tree.readMesh("cube.mesh");
	IrregularMesh<DIM> irregular_mesh(h_tree);
	irregular_mesh.globalRefine(2);

	double l = 0.1;
	do {
		irregular_mesh.semiregularize();
		irregular_mesh.regularize(false);
		RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
		regular_mesh.writeOpenDXData("cube.dx");
		std::cout << "\nPress ENTER to continue or CTRL+C to stop ..." << std::flush;
		getchar();

		Indicator<DIM> indicator(regular_mesh);
		AFEPack::Point<DIM> c0(0.25, 0.25, 0.5 - 1.5*l);
		AFEPack::Point<DIM> c1(0.25, 0.25, 0.5 + 1.5*l);
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
					Assert(false, ExcInternalError());
			}
			double d0 = 2*l;
			double d1 = 2*l;
			for (int j = 0;j < 4;j ++) {
			    AFEPack::Point<DIM>& p = regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[j]));
				double l1 = (p - c0).length();
				d0 = std::min(d0, l1);
				l1 = (p - c1).length();
				d1 = std::min(d1, l1);
			}
			double area = tetrahedron_volume(
					regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[0])),
					regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[1])),
					regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[2])),
					regular_mesh.point(regular_mesh.geometry(DIM,i).vertex(vertex[3])));
			if (d0 < l)
				indicator[i] = area;// * (l - d0);
			else if (d1 < l)
				indicator[i] = area;// * (l - d1);
			else
				indicator[i] = 0.;
		};

		if (l > 0.002) l /= 1.5;
		
		MeshAdaptor<DIM> mesh_adaptor(irregular_mesh);
		mesh_adaptor.convergenceOrder() = 0.;
		mesh_adaptor.refineStep() = 1;
		mesh_adaptor.setIndicator(indicator);
		mesh_adaptor.tolerence() = 1.0e-12;
		mesh_adaptor.adapt();
	} while (1);	
};

#undef DIM

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
