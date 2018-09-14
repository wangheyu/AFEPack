////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <stdio.h>
#include <dlfcn.h>

#include <iostream>
#include <fstream>

#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>

int main(int argc, char * argv[])
{
	Mesh<2,2> mesh;
	mesh.readData("L.mesh");
	mesh.renumerateElement();

	TemplateGeometry<2>				triangle_template_geometry;
	triangle_template_geometry.readData("triangle.tmp_geo");
	CoordTransform<2,2>				triangle_coord_transform;
	triangle_coord_transform.readData("triangle.crd_trs");
	TemplateDOF<2>					triangle_template_dof(triangle_template_geometry);
	triangle_template_dof.readData("triangle.1.tmp_dof");
	BasisFunctionAdmin<double,2,2>	triangle_basis_function(triangle_template_dof);
	triangle_basis_function.readData("triangle.1.bas_fun");

	TemplateGeometry<2>				rectangle_template_geometry;
	rectangle_template_geometry.readData("rectangle.tmp_geo");
	CoordTransform<2,2>				rectangle_coord_transform;
	rectangle_coord_transform.readData("rectangle.crd_trs");
	TemplateDOF<2>					rectangle_template_dof(rectangle_template_geometry);
	rectangle_template_dof.readData("rectangle.1.tmp_dof");
	BasisFunctionAdmin<double,2,2>	rectangle_basis_function(rectangle_template_dof);
	rectangle_basis_function.readData("rectangle.1.bas_fun");

	std::vector<TemplateElement<double,2,2> > template_element(2);
	template_element[0].reinit(triangle_template_geometry,
		triangle_template_dof,
		triangle_coord_transform,
		triangle_basis_function);
	template_element[1].reinit(rectangle_template_geometry,
		rectangle_template_dof,
		rectangle_coord_transform,
		rectangle_basis_function);

	FEMSpace<double,2> fem_space(mesh, template_element);
	
	int n_element = mesh.n_geometry(2);
	fem_space.element().resize(n_element);
	for (int i = 0;i < n_element;i ++) {
		int n_vertex = mesh.geometry(2,i).n_vertex();
		if (n_vertex == 3)
			fem_space.element(i).reinit(fem_space,i,0);
		else if (n_vertex == 4)
			fem_space.element(i).reinit(fem_space,i,1);
		else
			Assert(false, ExcNotImplemented());
	}

	fem_space.buildElement();
	fem_space.buildDof();
	fem_space.buildDofBoundaryMark();

	StiffMatrix<2,double> stiff_matrix(fem_space);
	stiff_matrix.algebricAccuracy() = 3;
	stiff_matrix.build();
	
	std::ofstream os("stiff_matrix.gnuplot");
	stiff_matrix.get_sparsity_pattern().print_gnuplot(os);
	return 0;
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
