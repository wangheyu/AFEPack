/******************************************************************************
 * easymesh_refine.cpp : by R. Lie, Jan 27, 2003
 * uniform refine a easymesh data file to get another mesh
 * ****************************************************************************
Usage: 
	easymesh_refine input_mesh output_mesh [refine_round = 1]
*/

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

int main(int argc, char * argv[])
{
	if (argc != 3 && argc != 4 && argc != 5) {
		std::cout << "Usage: " << argv[0]
			<< "input_mesh output_mesh [refine_round = 1 [is_renumerate = 0|1]]"
			<< std::endl;
		return 1;
	}
	HGeometryTree<2> h_geometry_tree;
	h_geometry_tree.readEasyMesh(argv[1]);
	int refine_round = 1;
	if (argc >= 4)
		refine_round = atoi(argv[3]);
	bool is_renumerate = false;
	if (argc >= 5) {
		if (atoi(argv[4]) == 1) {
			is_renumerate = true;
		}
	}
	IrregularMesh<2> irregular_mesh(h_geometry_tree);
	irregular_mesh.globalRefine(refine_round);
	irregular_mesh.semiregularize();
	irregular_mesh.regularize(is_renumerate);
	RegularMesh<2>& regular_mesh = irregular_mesh.regularMesh();
	regular_mesh.writeEasyMesh(argv[2]);
	return 0;
};

/*
 * end of file
 * ***************************************************************************/
