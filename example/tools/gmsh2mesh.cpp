/**
 * @file   gmsh2mesh.cpp
 * @author Robert Lie
 * @date   Tue Jun 22 16:21:37 2004
 * 
 * @brief  translate the gmsh mesh file format to internal mesh file format
 * 
 * 
 */

#include <AFEPack/Geometry.h>
#include <AFEPack/GmshMesh.h>

int main(int argc, char * argv[])
{
	if (argc != 3) {
		std::cout << "Usage: " << argv[0]
			<< " gmsh-file mesh-file"
			<< std::endl;
		exit(1);
	}
	GmshMesh gmsh;
	gmsh.readData(argv[1]);
	Mesh<3,3> mesh;
	gmsh.generateMesh(mesh);
	mesh.writeData(argv[2]);
	return 0;
};
