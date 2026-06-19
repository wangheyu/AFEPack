////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <iostream>
#include <fstream>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>

#define PI (4.0*atan(1.0))
#define DIM 3

double u(const double *);
double f(const double *);

int main(int argc, char * argv[])
{
  Mesh<DIM> mesh;
  mesh.readData(argv[1]);

  
  TemplateGeometry<DIM> tetrahedron_template_geometry;
  tetrahedron_template_geometry.readData("tetrahedron.tmp_geo");
  CoordTransform<DIM,DIM> tetrahedron_coord_transform;
  tetrahedron_coord_transform.readData("tetrahedron.crd_trs");
  TemplateDOF<DIM> tetrahedron_template_dof(tetrahedron_template_geometry);
  tetrahedron_template_dof.readData("tetrahedron.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> tetrahedron_basis_function(tetrahedron_template_dof);
  tetrahedron_basis_function.readData("tetrahedron.1.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > template_element(1);
  template_element[0].reinit(tetrahedron_template_geometry,
			     tetrahedron_template_dof,
			     tetrahedron_coord_transform,
			     tetrahedron_basis_function);

  FEMSpace<double,DIM> fem_space;
  fem_space.reinit(mesh, template_element);  
	
  int n_element = mesh.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space,i,0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  StiffMatrix<DIM,double> stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 2;
  stiff_matrix.build();

  FEMFunction<double,DIM> solution(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(&f, fem_space, right_hand_side, 4);

  BoundaryFunction<double,DIM> boundary1(BoundaryConditionInfo::DIRICHLET, 1, &u);
  BoundaryConditionAdmin<double,DIM> boundary_admin(fem_space);
  boundary_admin.add(boundary1);
  boundary_admin.apply(stiff_matrix, solution, right_hand_side);

  AMGSolver solver;
  solver.lazyReinit(stiff_matrix);
  solver.solve(solution, right_hand_side, 1.0e-08, 200);	

  solution.writeOpenDXData("u.dx");
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 2);
  std::cerr << "\nL2 error = " << error << std::endl;

  return 0;
};

double u(const double * p)
{
  return sin(PI*p[0]) * sin(2*PI*p[1]) * cos(3*PI*p[2]);
};

double f(const double * p)
{
  return 14.0 * PI * PI * u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
