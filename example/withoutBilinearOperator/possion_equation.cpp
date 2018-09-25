/**
 * @file   possion_equation.cpp
 * @author Heyu Wang <scshw@cslin107.csunix.comp.leeds.ac.uk>
 * @date   Mon May 19 13:19:19 2014
 * 
 * @brief 一个例子, 如何脱离 AFEPack 的 BilinearOperator 结构, 自己构建
 * 一个刚度矩阵. 目的是进一步在 NS 方程混合元求解等较为复杂的问题中直接
 * 构建大刚度矩阵. 甚至可以直接替换 deal.II 的矩阵结构. 建议对比
 * AFEPack 中 Possion 方程的例子看.
 * 
 * 
 */

#include <iostream>
#include <fstream>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_ilu.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_mic.h>
#include <lac/sparse_decomposition.h>

#define PI (4.0*atan(1.0))

#define DIM 2

double u(const double *);
double f(const double *);

int main(int argc, char * argv[])
{
    HGeometryTree<DIM> h_tree;

    h_tree.readEasyMesh(argv[1]);

    IrregularMesh<DIM> *irregular_mesh0;
    IrregularMesh<DIM> *irregular_mesh1;

    irregular_mesh0 = new IrregularMesh<DIM>;
    irregular_mesh0->reinit(h_tree);
    irregular_mesh0->semiregularize();
    irregular_mesh0->regularize(false);
    RegularMesh<DIM> &mesh0 = irregular_mesh0->regularMesh();

    irregular_mesh1 = new IrregularMesh<DIM>(*irregular_mesh0);
    irregular_mesh1->globalRefine(3);
    irregular_mesh1->semiregularize();
    irregular_mesh1->regularize(false);
    RegularMesh<DIM> &mesh1 = irregular_mesh1->regularMesh();
    
    mesh0.writeOpenDXData("D0.dx");
    mesh1.writeOpenDXData("D1.dx");

    TemplateGeometry<DIM> triangle_template_geometry;
    triangle_template_geometry.readData("triangle.tmp_geo");
    CoordTransform<DIM, DIM> triangle_coord_transform;
    triangle_coord_transform.readData("triangle.crd_trs");
    TemplateDOF<DIM> triangle_template_dof(triangle_template_geometry);
    triangle_template_dof.readData("triangle.3.tmp_dof");
    BasisFunctionAdmin<double, DIM, DIM> triangle_basis_function(triangle_template_dof);
    triangle_basis_function.readData("triangle.3.bas_fun");

    std::vector<TemplateElement<double, DIM, DIM> > template_element(1);
    template_element[0].reinit(triangle_template_geometry,
			       triangle_template_dof,
			       triangle_coord_transform,
			       triangle_basis_function);


    FEMSpace<double, DIM> fem_space0(mesh0, template_element);
	
    int n_element = mesh0.n_geometry(DIM);
    fem_space0.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
	fem_space0.element(i).reinit(fem_space0, i, 0);

    fem_space0.buildElement();
    fem_space0.buildDof();
    fem_space0.buildDofBoundaryMark();


    FEMSpace<double, DIM> fem_space1(mesh1, template_element);
	
    n_element = mesh1.n_geometry(DIM);
    fem_space1.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
    	fem_space1.element(i).reinit(fem_space1, i, 0);

    fem_space1.buildElement();
    fem_space1.buildDof();
    fem_space1.buildDofBoundaryMark();
    /// 到此为止, 和之前的做法毫无不同.

    /// 准备统计系数矩阵的每一行有多少个非零元.
    int n_total_dof = fem_space0.n_dof();
    std::vector<unsigned int> n_non_zero_per_row(n_total_dof);

    /// 准备一个遍历全部单元的迭代器.
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space0.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space0.endElement();

    /// 第一次循环遍历全部单元, 只是为了统计每一行的非零元个数.
    for (; the_element != end_element; ++the_element) 
    {
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
	for (int i = 0; i < n_element_dof; ++i)
	    n_non_zero_per_row[element_dof[i]] += n_element_dof;
    }

    /// 根据每一行的非零元创建矩阵模板.
    SparsityPattern sp_stiff_matrix(fem_space0.n_dof(), n_non_zero_per_row);

    /// 第二次遍历, 指定每个非零元的坐标.
    for (the_element = fem_space0.beginElement(); 
	 the_element != end_element; ++the_element) 
    {
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
	for (int i = 0; i < n_element_dof; ++i)
	    for (int j = 0; j < n_element_dof; ++j)
		sp_stiff_matrix.add(element_dof[i], element_dof[j]);
    }

    /// 矩阵模板压缩. 创建矩阵.
    sp_stiff_matrix.compress();
    SparseMatrix<double> stiff_matrix(sp_stiff_matrix);

    Vector<double> rhs(fem_space0.n_dof());
    /// 第三次遍历, 给每个矩阵元素正确的赋值. 这一块和
    /// BilinearOperator 中的 getElementMatrix 中做的事情一致.
    for (the_element = fem_space0.beginElement(); 
	 the_element != end_element; ++the_element) 
    {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(6);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<2> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<std::vector<double> > > basis_gradient = the_element->basis_function_gradient(q_point);
        std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
	for (int l = 0; l < n_quadrature_point; ++l)
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	    for (int i = 0; i < n_element_dof; ++i)
            {
		for (int j = 0; j < n_element_dof; ++j)
		{
		    double cont = Jxw * innerProduct(basis_gradient[i][l], basis_gradient[j][l]);
		    stiff_matrix.add(element_dof[i], element_dof[j], cont);
		}
                double cont_rhs = Jxw * f(q_point[l]) * basis_value[i][l];
		rhs(element_dof[i]) += cont_rhs;
            }
	}
    }

    /// 接下去做的事情和之前一样.
    FEMFunction<double, DIM> solution0(fem_space0);

    // for (int i = 0; i < stiff_matrix.m(); i++)
    // {
    // 	SparseMatrix<double>::iterator row_iterator = stiff_matrix.begin(i);
    // 	SparseMatrix<double>::iterator row_end = stiff_matrix.end(i);
    
    // 	for ( ;row_iterator != row_end; ++row_iterator)
    // 		std::cout << "(" << i << ", " << row_iterator->column() << ") = " << row_iterator->value() << ",";
    //     std::cout << std::endl;
    // }	

    for (int i = 0; i < fem_space0.n_dof(); i++)
    {
    	FEMSpace<double, DIM>::dof_info_t dof = fem_space0.dofInfo(i);
    	if (dof.boundary_mark == 1)
    	{
    	    SparseMatrix<double>::iterator row_iterator = stiff_matrix.begin(i);
    	    SparseMatrix<double>::iterator row_end = stiff_matrix.end(i);
    	    double diag = row_iterator->value();
    	    double bnd_value = u(dof.interp_point); 
            rhs(i) = diag * bnd_value;
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
    		rhs(k) -= col_iterator->value() * bnd_value; 
    		col_iterator->value() = 0.0;	
            }  
            
    	}	
    }		

   // BoundaryFunction<double, DIM> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
   // BoundaryConditionAdmin<double, DIM> boundary_admin(fem_space0);
   // boundary_admin.add(boundary);
   // boundary_admin.apply(stiff_matrix, solution0, rhs);

    // for (int i = 0; i < stiff_matrix.m(); i++)
    // {
    // 	SparseMatrix<double>::iterator row_iterator = stiff_matrix.begin(i);
    // 	SparseMatrix<double>::iterator row_end = stiff_matrix.end(i);
    
    // 	for ( ;row_iterator != row_end; ++row_iterator)
    // 		std::cout << "(" << i << ", " << row_iterator->column() << ") = " << row_iterator->value() << ",";
    //     std::cout << std::endl;
    // }	


//    for (int i = 0; i < fem_space0.n_dof(); i++)
//        std::cout << "rhs[" << i << "]=" << right_hand_side(i) << "=" << rhs(i) << std::endl;
	




//    SparseILU<double> preconditioner(stiff_matrix);
    SparseMIC<double> mic;
    
    SparseMIC<double>::AdditionalData ad;
    mic.initialize(stiff_matrix, ad);

    PreconditionSSOR<SparseMatrix<double> > ssor;
    ssor.initialize (stiff_matrix, 1.2);

    SolverControl solver_control (200000, 1e-15);
    SolverCG<Vector<double> > cg (solver_control);
//    cg.solve (stiff_matrix,  solution0,  rhs, PreconditionIdentity());
    cg.solve (stiff_matrix,  solution0,  rhs, ssor);
//    cg.solve (stiff_matrix,  solution0,  rhs, mic);

//     AMGSolver solver(stiff_matrix);
//     solver.solve(solution0, rhs);	
    solution0.writeOpenDXData("u0.dx");


    FEMFunction<double, DIM> solution1(fem_space1);
	    

    Operator::L2Interpolate(solution0, solution1);
//    Operator::L2Project(solution0, solution1, 3, 3);
    solution1.writeOpenDXData("u1.dx");

    double error = Functional::L2Error(solution0, FunctionFunction<double>(&u), 10);
    std::cout << "L2 Error = " << error << std::endl;
    return 0;
};

double u(const double * p)
{
    return sin(PI * p[0]) * sin(2.0 * PI * p[1]);
};

double f(const double * p)
{
    return 5 * PI * PI * u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
