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
#include <AFEPack/EasyMesh.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>


#define PI (4.0*atan(1.0))

double u(const double *);
double f(const double *);

int main(int argc, char * argv[])
{
        Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 4);

  EasyMesh mesh;
  mesh.readData(argv[1]);

  TemplateGeometry<2>	triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<2,2>	triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<2>	triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("triangle.2.tmp_dof");
  BasisFunctionAdmin<double,2,2> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData("triangle.2.bas_fun");

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

  StiffMatrix<2,double> stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 2;
  stiff_matrix.build();

  FEMFunction<double,2> solution(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(&f, fem_space, right_hand_side, 4);

    for (int i = 0; i < fem_space.n_dof(); i++)
    {
    	FEMSpace<double, 2>::dof_info_t dof = fem_space.dofInfo(i);
    	if (dof.boundary_mark == 1)
    	{
    	    SparseMatrix<double>::iterator row_iterator = stiff_matrix.begin(i);
    	    SparseMatrix<double>::iterator row_end = stiff_matrix.end(i);
    	    double diag = row_iterator->value();
    	    double bnd_value = u(dof.interp_point); 
            right_hand_side(i) = diag * bnd_value;
    	    for ( ++row_iterator; row_iterator != row_end; ++row_iterator)
            {
            	row_iterator->value() = 0.0;
    		int k = row_iterator->column();
                SparseMatrix<double>::iterator col_iterator = stiff_matrix.begin(k);   
                SparseMatrix<double>::iterator col_end = stiff_matrix.end(k);   
    	    	for (++col_iterator; col_iterator != col_end; ++col_iterator)
    			if (col_iterator->column() == i)
    			    break;
    		if (col_iterator == col_end)
    		{
    			std::cerr << "Error!" << std::endl;
    			exit(-1);
    		}
    		right_hand_side(k) -= col_iterator->value() * bnd_value; 
    		col_iterator->value() = 0.0;	
            }  
            
    	}	
    }		


//  BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
//  BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
//  boundary_admin.add(boundary);
//  boundary_admin.apply(stiff_matrix, solution, right_hand_side);

    // 明显弄个CG法.
    SolverControl           solver_control (100000, 1e-15);
    // 测试求解精度.
    SolverCG<>              solver (solver_control);

    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner;
    Amg_preconditioner = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG>(new TrilinosWrappers::PreconditionAMG());

    Amg_preconditioner->initialize(stiff_matrix);
    solver.solve (stiff_matrix, solution, right_hand_side,
		  *Amg_preconditioner);	
//  AMGSolver solver(stiff_matrix);
//  solver.solve(solution, right_hand_side, 1.0e-08, 2000);	

  solution.writeOpenDXData("u.dx");
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 4);
  std::cerr << "\nL2 error = " << error << std::endl;

  return 0;
};

double u(const double * p)
{
  return sin(PI*p[0]) * sin(2*PI*p[1]);
};

double f(const double * p)
{
  return 5*PI*PI*u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
