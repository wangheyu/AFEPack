/**
 * @file   AMGSolver.h
 * @author Robert Lie
 * @date   Sun Mar 13 16:25:11 2002
 * 
 * @brief  algebraic multigrid solver for positive defined system.
 * 
 * 
 */

#ifndef _AMGSolver_h_
#define _AMGSolver_h_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <list>
#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include "Miscellaneous.h"

AFEPACK_OPEN_NAMESPACE

//! An algebraic multigrid solver.
/**
 * This is an implelementation of the algebraic multigrid solver for a general 
 * sparse matrix. The projection operators include a very simple element-free type 
 * operator and a comparative complex one. The simple one can be used for solving
 * once for a sequence of systems with the same sparsity pattern(often in a time
 * depended problem). It is faster in preparation but slower in solving.
 * The mode can be accesed by the so-called lazy-mode. The complexier one can be used
 * to solve one system for many times, since it is slower in preparation but faster
 * in solving. The smoother is a Gauss-Sidel iteration. For the two projection
 * operator, please refer Brezina etc, SIAM JSC 22(5), pp1570-1592 and Cleary etc, 
 * SIAM JSC 21(5), pp1888-1908.
 *
 * The usage of this solver is quite simple. The following is an example:
 * <pre>
 * SparseMatrix<double> M;
 * Vector<double> x;
 * Vector<double> rhs;
 * ... ...
 * AMGSolver amg_solver(M);
 * amg_solver.solve(x, rhs);
 * </pre>
 * For the lazy mode, we often want to solve a sequence of systems as follows:
 * <pre>
 * SparseMatrix<double> M;
 * Vector<double> x;
 * Vector<double> rhs;
 * AMGSolver amg_solver;
 * ... ...
 * while (1) {
 *   amg_solver.lazyReinit(M);
 *   amg_solver.solve(x, rhs);
 *   ... ... // made some changes to the matrix
 * }
 * </pre>
 *
 */
class AMGSolver {
 public:
  typedef SparseMatrix<double> Matrix;
 private:
  bool is_initialized; /// if the solver is initialized
  u_int n_project; /// number of grid levels
  std::vector<Matrix *> project_matrix; /// the left projection matrix
  std::vector<Matrix *> project_matrix_r; /// the right projection matrix, which is the transpose of the left projection matrix
  std::vector<const Matrix *> projected_matrix; /// the projected matrix
  bool is_most_project_full; /// if the last projected matrix is full
  bool is_solve_most_project_exactly; /// if use gauss-elimination to solve the most projected system exactly
  FullMatrix<double> M_n; /// if the last projected matrix is full, this is the inverse of this matrix
  double toler; /// tolerence
  double alpha; /// the parameter for given C-point
  u_int smooth_step; /// GS smooth step in every 
  u_int n_most_project_order; /// lower bound of the order of the last projected matrix
  double n_most_project_sparse_degree; ///  upper bound of the sparse degree of the last projected matrix
 public:
  AMGSolver(); /// default constructor
  AMGSolver(const Matrix&, 
	    double tol = 1.0e-12, 
	    u_int s = 3, 
	    u_int nmpo = 50,
	    double nmpsd = 0.382,
	    double alp = 0.25);
  virtual ~AMGSolver(); /// destructor
  void clear(); /// clear memory used by the solver

  double tolerence() const {return toler;}; /// 
  double& tolerence() {return toler;};
  u_int smoothStep() const {return smooth_step;};
  u_int& smoothStep() {return smooth_step;};
  bool& isSolveMostProjectExactly() {return is_solve_most_project_exactly;}
  bool isSolveMostProjectExactly() const {return is_solve_most_project_exactly;}

  void reinit(const Matrix&); /// reinitialize using another matrix
  /** 
   * solve the linear system with $x$ as solution and $r$ as right hand side.
   * 
   * @param x the solution vector
   * @param r the right hand side vector
   * @param tol tolerence
   * @param step upper bound of iteration step
   * @param mode 0: solve mode; 1: precondition mode
   *
   */
  void solve(Vector<double>& x, 
	     const Vector<double>& r, 
	     double tol = 0.0, 
	     u_int step = 20,
	     int mode = 0) const; 

  void lazyReinit(const Matrix&); /// reinitialize using another matrix in lazy mode

 private:
  void lazyInit(const Matrix&);
  void lazyProject(const Matrix& M, 
		   Matrix *& P, 
		   Matrix *& PMPt,
		   Matrix *& Pt);

  void Project(const Matrix& M, 
	       Matrix *& P, 
	       Matrix *& PMPt);

  void GaussSidel(const Matrix& M, 
		  Vector<double>& x, 
		  const Vector<double>& r, 
		  const int& s) const;
  Matrix * getPMPt(const Matrix & P, 
		   const Matrix & M,
		   const Matrix & Pt) const;
};



//! An algebraic multigrid preconditioner
/**
 * With this class, we can use AMG solver as the proconditioner of deal.II
 * As a class to be a preconditioner, the only requirement is to have a
 * method as
 * <pre>
 * void vmult(Vector<double>&, const Vector<double>&);
 * </pre>
 * We here modified the solve procedure of the AMGSolver to iterate a given
 * steps and give it the name "vmult".
 *
 * This is an example to use this preconditioner
 * <pre>
 * SparseMatrix<double> M;
 * Vector<double> x;
 * Vector<double> rhs;
 * ... ...
 * AMGPreconditioner amg(M);
 * SolverControl solver_control (50, 1.0e-8);
 * PrimitiveVectorMemory<> vector_memory;
 * SolverBicgstab<> solver(solver_control, vector_memory);
 * solver.solve(M, x, rhs,
 *              PreconditionUseMatrix<AMGPreconditioner, Vector<double> >
 *              (amg, &AMGPreconditioner::vmult));
 * </pre>
 *
 */
class AMGPreconditioner : public AMGSolver
{
 private:
  u_int iterate_step; /// iteration steps in every precondition multiplication
 public:
  AMGPreconditioner() {};
  AMGPreconditioner(const Matrix& M, 
		    u_int is = 5, 
		    u_int s = 3) : 
    iterate_step(is) {
    lazyReinit(M);
    smoothStep() = s;
  }
  ~AMGPreconditioner() {}
 public:
  void reinit(const Matrix& M,
	      u_int is = 5) {
    AMGSolver::lazyReinit(M);
    iterate_step = is;
  }
  u_int& iterateStep() {return iterate_step;}
  const u_int& iterateStep() const {return iterate_step;}
  void vmult(Vector<double>& x, 
	     const Vector<double>& b) const
    {solve(x, b, 0, iterate_step, 1);}
};

AFEPACK_CLOSE_NAMESPACE

#endif

//
// end of file
///////////////////////////////////////////////////////////////////////////////

