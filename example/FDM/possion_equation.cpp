#include <iostream>
#include <cmath>

#include <AFEPack/Geometry.h>

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_ilu.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_mic.h>
#include <lac/sparse_decomposition.h>

#define DIM 2
#define EPS 1e-7
#define PI (atan(1.0) * 4.0)

double u(const double * p);
double f(const double * p);

class FDMDomain
{
public:
    typedef AFEPack::Point<DIM> pnt_t;
private:
    pnt_t lb_conner; 
    pnt_t rt_conner; 
    pnt_t center;
    double width;
    double height;
    double hx;
    double hy;
    u_int Nx;
    u_int Ny;
public:
    ~FDMDomain()
	{};
    FDMDomain()
	{};
    FDMDomain(pnt_t &lb, pnt_t &rt, double h);
    FDMDomain(pnt_t &lb, pnt_t &rt, int n);
    FDMDomain(pnt_t &lb, pnt_t &rt, int nx, int ny);
    int n_dofs();
    int list_pnts();
    int get_Nx();
    int get_Ny();
    double get_hx();
    double get_hy();
    pnt_t interp_point(int k);
    pnt_t interp_point(int i, int j);
    int index_mesh2dof(int i, int j);
    int index_dof2x(int n);
    int index_dof2y(int n);
    int central5(SparsityPattern &spA, SparseMatrix<double> &A);
    int apply_Dirichlet(SparseMatrix<double> &A,
			Vector<double> &b);
};

FDMDomain::FDMDomain(pnt_t &lb, pnt_t &rt, double h)
{
    lb_conner = lb;
    rt_conner = rt;
    center = midpoint(lb, rt);
    width = rt[0] - lb[0];
    height = rt[1] - lb[1];
    hx = h;
    hy = h;
    Nx = u_int(width / hx) + 1;
    Ny = u_int(height / hy) + 1;
    if (std::fabs(width - (Nx - 1) * hx) > EPS)
    {
	std::cerr << "Waring: grid space is not match in x direction."
		  << std::endl;
    }
    if (std::fabs(width - (Ny - 1) * hy) > EPS)
    {
	std::cerr << "Waring: grid space is not match in y direction."
		  << std::endl;
    }
};

FDMDomain::FDMDomain(pnt_t &lb, pnt_t &rt, int n)
{
    lb_conner = lb;
    rt_conner = rt;
    center = midpoint(lb, rt);
    width = rt[0] - lb[0];
    height = rt[1] - lb[1];
    hx = width / (n - 1);
    hy = height / (n - 1);
    Nx = n;
    Ny = n;
    if (std::fabs(hx - hy) > EPS)
    {
	std::cerr << "Waring: grid spaces are not equal to each other between x and y directions."
		  << std::endl;
    }
};

FDMDomain::FDMDomain(pnt_t &lb, pnt_t &rt, int nx, int ny)
{
    lb_conner = lb;
    rt_conner = rt;
    center = midpoint(lb, rt);
    width = rt[0] - lb[0];
    height = rt[1] - lb[1];
    hx = width / (nx - 1);
    hy = height / (ny - 1);
    Nx = nx;
    Ny = ny;
    if (std::fabs(hx - hy) > EPS)
    {
	std::cerr << "Waring: grid spaces are not equal to each other between x and y directions."
		  << std::endl;
    }
};

int FDMDomain::n_dofs()
{
    return Nx * Ny;
};

int FDMDomain::get_Nx()
{
    return Nx;
};

int FDMDomain::get_Ny()
{
    return Ny;
};

double FDMDomain::get_hx()
{
    return hx;
};

double FDMDomain::get_hy()
{
    return hy;
};

int FDMDomain::index_mesh2dof(int i, int j)
{
    return (i * Nx + j);
};

int FDMDomain::index_dof2x(int n)
{
    return n % Nx;
};

int FDMDomain::index_dof2y(int n)
{
    return n / Nx;
};

int FDMDomain::list_pnts()
{
    for (int i = 0; i < Ny; i++)
    {
	for (int j = 0; j < Nx; j++)
	{
	    double x = lb_conner[0] + hx * j;
	    double y = lb_conner[1] + hy * i;
	    std::cout.setf(std::ios::fixed);
	    std::cout.precision(7);
	    std::cout << "(" << x << ", " << y <<")" << "\t";
	}
	std::cout << std::endl;
    }
    return 0;
};

FDMDomain::pnt_t FDMDomain::interp_point(int k)
{
    int i = k / Nx;
    int j = k % Nx;
    pnt_t x(lb_conner[0] + hx * j, lb_conner[1] + hy * i);
    return x;
};

FDMDomain::pnt_t FDMDomain::interp_point(int i, int j)
{
    return FDMDomain::pnt_t(lb_conner[0] + hx * j, lb_conner[1] + hy * i);
};

int FDMDomain::central5(SparsityPattern &spA, SparseMatrix<double> &A)
{
    int n = this->n_dofs();
    spA.reinit(n, n, 5);
    
    for (int k = 0; k < n; k ++)
	spA.add(k, k);
    for (int j = 0; j < Nx - 1; j++)
	for (int i = 0; i < Ny; i++)
	{
	    int k = i * Nx + j;
	    spA.add(k, k + 1);
	}
    for (int j = 1; j < Nx; j++)
	for (int i = 0; i < Ny; i++)
	{
	    int k = i * Nx + j;
	    spA.add(k, k - 1);
	}
    for (int j = 0; j < Nx; j++)
	for (int i = 0; i < Ny - 1; i++)
	{
	    int k = i * Nx + j;
	    spA.add(k, k + Nx);
	}
    for (int j = 0; j < Nx; j++)
	for (int i = 1; i < Ny; i++)
	{
	    int k = i * Nx + j;
	    spA.add(k, k - Nx);
	}

    spA.compress();
    A.reinit(spA);
    for (int k = 0; k < n; k ++)
	A.set(k, k, 4.0);
    for (int j = 0; j < Nx - 1; j++)
	for (int i = 0; i < Ny; i++)
	{
	    int k = i * Nx + j;
	    A.set(k, k + 1, -1.0);
	}
    for (int j = 1; j < Nx; j++)
	for (int i = 0; i < Ny; i++)
	{
	    int k = i * Nx + j;
	    A.set(k, k - 1, -1.0);
	}
    for (int j = 0; j < Nx; j++)
	for (int i = 0; i < Ny - 1; i++)
	{
	    int k = i * Nx + j;
	    A.set(k, k + Nx, -1.0);
	}
    for (int j = 0; j < Nx; j++)
	for (int i = 1; i < Ny; i++)
	{
	    int k = i * Nx + j;
	    A.set(k, k - Nx, -1.0);
	}
    return 0;
};

int FDMDomain::apply_Dirichlet(SparseMatrix<double> &A,
			       Vector<double> &b)
{
    int n = this->n_dofs();
    for (int k = 0; k < n; k++)
    {
	int j = this->index_dof2x(k);
	int i = this->index_dof2y(k);
	if (i == 0 || i == Ny - 1 || j == 0 || j == Nx - 1)
	{
	    SparseMatrix<double>::iterator row_iterator = A.begin(k);
    	    SparseMatrix<double>::iterator row_end = A.end(k);
    	    double diag = row_iterator->value();
    	    double bnd_value = u(this->interp_point(k)); 
            b(i) = diag * bnd_value;
    	    for (++row_iterator; row_iterator != row_end; ++row_iterator)
            {
            	row_iterator->value() = 0.0;
    		int c = row_iterator->column();
                SparseMatrix<double>::iterator col_iterator = A.begin(c);   
                SparseMatrix<double>::iterator col_end = A.end(c);   
    	    	for (++col_iterator; col_iterator != col_end; ++col_iterator)
    			if (col_iterator->column() == k)
    			    break;
    		if (col_iterator == col_end)
    		{
    			std::cerr << "Error!" << std::endl;
    			exit(-1);
    		}
    		b(c) -= col_iterator->value() * bnd_value; 
    		col_iterator->value() = 0.0;	
            }  
	}
    }
};

int main (int argc, char *argv[])
{
    int n = 100;
    FDMDomain::pnt_t x0(1.0, 1.0); 
    FDMDomain::pnt_t x1(3.0, 3.0); 
    FDMDomain D0(x0, x1, n);
    SparsityPattern spA;
    SparseMatrix<double> A;
    D0.central5(spA, A);
    int n_dofs = D0.n_dofs();
    Vector<double> rhs(n_dofs);
    for (int k = 0; k < n_dofs; k++)
    {
	double h = D0.get_hx();
	rhs(k) = f(D0.interp_point(k)) * h * h;
    }
    Vector<double> solution(n_dofs);
    D0.apply_Dirichlet(A, rhs);

    // 明显弄个CG法.
    SolverControl           solver_control (10000, 1e-15);
    // 测试求解精度.
    SolverCG<>              solver (solver_control);
    solver.solve (A, solution, rhs, PreconditionIdentity());

    double error = 0;
    for (int k = 0; k < n_dofs; k++)
    {
	double diff = fabs(u(D0.interp_point(k)) - solution(k));
	error += diff * diff;
    }
    std::cout << "h = " << D0.get_hx() << std::endl;
    std::cerr << "\nL2 error = " << sqrt(error) / n_dofs  << std::endl;

    return 0;
};

double u(const double * p)
{
  return sin(p[0]) * sin(p[1]);
};

double f(const double * p)
{
  return 2.0 * u(p);
};

