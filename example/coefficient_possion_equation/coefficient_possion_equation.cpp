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
#include <AFEPack/EasyMesh.h>

double u(const double *);
double f(const double *);


double a11_(const double * p)
{
	if (p[0] < p[1])
		return 10.0 + p[0];
	else
		return 2.0 + p[1];
};

double a12_(const double * p)
{
	return 1.0;
};

double a21_(const double * p)
{
	return 2.0;
};

double a22_(const double * p)
{
	if (p[0] + p[1] > 1)
		return 10.0 + p[0];
	else
		return 1.0 + p[1];
};

class Matrix : public StiffMatrix<2, double>
{
	public:
		Matrix(FEMSpace<double,2>& sp) : 
			StiffMatrix<2,double>(sp) {};
		virtual ~Matrix() {};
public:
	virtual void getElementMatrix(const Element<double,2>& element0,
		const Element<double,2>& element1,
		const ActiveElementPairIterator<2>::State state);
};

void Matrix::getElementMatrix(
		const Element<double,2>& element0,
		const Element<double,2>& element1,
		const ActiveElementPairIterator<2>::State state)
{
	int n_element_dof0 = elementDof0().size();
	int n_element_dof1 = elementDof1().size();
	double volume = element0.templateElement().volume();
	const QuadratureInfo<2>& quad_info = element0.findQuadratureInfo(algebricAccuracy());
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<2> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<std::vector<double> > > basis_gradient = element0.basis_function_gradient(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
		double Jxw = quad_info.weight(l)*jacobian[l]*volume;
		AFEPack::Point<2> q_point = element0.local_to_global(quad_info.quadraturePoint(l));
		double a11 = a11_(q_point);
		double a12 = a12_(q_point);
		double a21 = a21_(q_point);
		double a22 = a22_(q_point);
		for (int j = 0;j < n_element_dof0;j ++) {
			for (int k = 0;k < n_element_dof1;k ++) {
				elementMatrix(j,k) += Jxw*(
					a11 * basis_gradient[j][l][0] * basis_gradient[k][l][0]
					+ a12 * basis_gradient[j][l][0] * basis_gradient[k][l][1]
					+ a21 * basis_gradient[j][l][1] * basis_gradient[k][l][0]
					+ a22 * basis_gradient[j][l][1] * basis_gradient[k][l][1]
				);
			}
		}
	}
};



int main(int argc, char * argv[])
{
	EasyMesh mesh;
	mesh.readData("D");

	TemplateGeometry<2>				triangle_template_geometry;
	triangle_template_geometry.readData("triangle.tmp_geo");
	CoordTransform<2,2>				triangle_coord_transform;
	triangle_coord_transform.readData("triangle.crd_trs");
	TemplateDOF<2>					triangle_template_dof(triangle_template_geometry);
	triangle_template_dof.readData("triangle.1.tmp_dof");
	BasisFunctionAdmin<double,2,2>	triangle_basis_function(triangle_template_dof);
	triangle_basis_function.readData("triangle.1.bas_fun");

	std::vector<TemplateElement<double,2,2> > template_element(1);
	template_element[0].reinit(triangle_template_geometry,
		triangle_template_dof,
		triangle_coord_transform,
		triangle_basis_function);

	FEMSpace<double,2> fem_space(mesh, template_element);
	
	int n_element = mesh.n_geometry(2);
	fem_space.element().resize(n_element);
	for (int i = 0;i < n_element;i ++)
		fem_space.element(i).reinit(fem_space,i,0);

	fem_space.buildElement();
	fem_space.buildDof();
	fem_space.buildDofBoundaryMark();

	Matrix stiff_matrix(fem_space);
	stiff_matrix.algebricAccuracy() = 3;
	stiff_matrix.build();

	FEMFunction<double,2> solution(fem_space);
	Vector<double> right_hand_side;
	Operator::L2Discretize(&f, fem_space, right_hand_side, 3);

	BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
	BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
	boundary_admin.add(boundary);
	boundary_admin.apply(stiff_matrix, solution, right_hand_side);

	AMGSolver solver(stiff_matrix);
	solver.solve(solution, right_hand_side);	

	solution.writeOpenDXData("u.dx");
	return 0;
};

double u(const double * p)
{
	return sin(p[0]) * sin(p[1]);
};

double f(const double * p)
{
	return 1.0;
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
