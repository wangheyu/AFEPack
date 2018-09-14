/**
 * @file   mesh_refine.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Thu Jan  1 10:35:09 2009
 * 
 * @brief  
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

int main(int argc, char * argv[])
{
	if (argc != 4 && argc != 5 && argc != 6) {
		std::cout << "Usage: " << argv[0]
			<< " dimension input_mesh output_mesh [refine_round=1 [is_renum]]"
			<< std::endl;
		return 1;
	}
        int dimension = atoi(argv[1]);
        if (dimension == 2) {
#define DIM 2
          HGeometryTree<DIM> h_geometry_tree;
          h_geometry_tree.readMesh(argv[2]);
          int refine_round = 1;
          if (argc >= 5)
            refine_round = atoi(argv[4]);
          bool is_renum = true;
          if (argc >= 6)
            is_renum = false;
          IrregularMesh<DIM> irregular_mesh(h_geometry_tree);
          irregular_mesh.globalRefine(refine_round);
          irregular_mesh.semiregularize();
          irregular_mesh.regularize(is_renum);
          RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
          regular_mesh.writeData(argv[3]);
#undef DIM
        } else if (dimension == 3) {
#define DIM 3
          HGeometryTree<DIM> h_geometry_tree;
          h_geometry_tree.readMesh(argv[2]);
          int refine_round = 1;
          if (argc >= 5)
            refine_round = atoi(argv[4]);
          bool is_renum = true;
          if (argc >= 6)
            is_renum = false;
          IrregularMesh<DIM> irregular_mesh(h_geometry_tree);
          irregular_mesh.globalRefine(refine_round);
          irregular_mesh.semiregularize();
          irregular_mesh.regularize(is_renum);
          RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();
          regular_mesh.writeData(argv[3]);
#undef DIM
        }
	return 0;
};

/*
 * end of file
 * ***************************************************************************/
