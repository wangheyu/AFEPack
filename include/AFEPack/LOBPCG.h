/**
 * @file   LOBPCG.h
 * @author Hu Guanghui <gary@Brin>
 * @date   Fri Nov 18 16:07:26 2016
 * 
 * @brief  A * \psi = \lambda * B * \psi
 * 
 * 
 */
#ifndef _LOBPCG_H_
#define _LOBPCG_H_

#include <AFEPack/AMGSolver.h>
#include <AFEPack/BilinearOperator.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>

//#include <lac/sparse_matrix.h>
#include <AFEPack/SparseMatrix.h>
//#include <lac/lapack_full_matrix.h>
#include <AFEPack/DenseMatrix.h>
//#include <lac/vector.h>
#include <AFEPack/Vector.h>

#include <iomanip>
#include <omp.h>
#include <fstream>

//#define DIM 3

class LOBPCG
{
 public:
  LOBPCG();
  ~LOBPCG();
 private:
  SparseMatrix<double> *A;
  SparseMatrix<double> *B;
  //StiffMatrix<DIM, double> *stiff_matrix;
  //SparseMatrix<double> *stiff_matrix;
  //SparseMatrix<double> D;
  int n_eigenpair;
  std::vector<Vector<double>* >  eigenVector;
  std::vector<Vector<double> > residual;
  std::vector<Vector<double> > precondResidual;  
  std::vector<Vector<double> > searchDirection;
  std::vector<double>* eigenValue;
  std::vector<Vector<double> > AVec, BVec, RRVector;
  std::vector<Vector<double> * > tmpVec;/// for orthogonalization

  ///LAPACKFullMatrix for RR procedure
  //LAPACKFullMatrix<double> * lapackA, * lapackB, * lapackTMP, * lapackTMP1, * lapackTMP2, * lapackTMP3, * lapackTMP4, * lapackTMP5, * lapackTMP6, * lapackA0, * lapackB0;
  DenseMatrix<double> * lapackA, * lapackB, * lapackTMP, * lapackTMP1, * lapackTMP2, * lapackTMP3, * lapackTMP4, * lapackTMP5, * lapackTMP6, * lapackA0, * lapackB0;
  
public:
  //void reinit(SparseMatrix<double>&, MassMatrix<DIM, double>&, std::vector<FEMFunction<double, DIM> >&, std::vector<double>&, AMGSolver&);
  void reinit(SparseMatrix<double>&, SparseMatrix<double>&, std::vector<Vector<double> >&, std::vector<double>&);
 // void reinit(SparseMatrix<double>&, SparseMatrix<double>&, std::vector<Vector<double> >&, std::vector<double>&, AMGSolver&);  
  void Borthonormalize(SparseMatrix<double>&, std::vector<Vector<double>* >&);
  void normalize(std::vector<Vector<double>* >&);

  /*We need the four subfunctions below to cover all 4 possible calculations */
  void getPMQt(DenseMatrix<double>&,
	       const std::vector<Vector<double>* >&,
	       const SparseMatrix<double>&,
	       const std::vector<Vector<double> >&);
  void getPMQt(DenseMatrix<double>&,
	       const std::vector<Vector<double> >&,
	       const SparseMatrix<double>&,
	       const std::vector<Vector<double>* >&);
  void getPMQt(DenseMatrix<double>&,
	       const std::vector<Vector<double>* >&,
	       const SparseMatrix<double>&,
	       const std::vector<Vector<double>* >&);
  void getPMQt(DenseMatrix<double>&,
	       const std::vector<Vector<double> >&,
	       const SparseMatrix<double>&,
	       const std::vector<Vector<double> >&);
  /* RR procedure, first one for the first step, the second one for the rest steps*/
  void RRProcedure(std::vector<Vector<double> >&,
		   const std::vector<Vector<double>* >&,
		   const std::vector<Vector<double> >&,
		   const SparseMatrix<double>&,
		   const SparseMatrix<double>&);
  void RRProcedure(std::vector<Vector<double> >&,
		   const std::vector<Vector<double>* >&,
		   const std::vector<Vector<double> >&,
		   const std::vector<Vector<double> >&,
		   const SparseMatrix<double>&,
		   const SparseMatrix<double>&);

  void catMatrix(DenseMatrix<double>&,
		 const DenseMatrix<double>&,
		 const int &,
		 const int &,
		 const int &);
  void solve(const double& tol_LOBPCG,
	     const int& max_LOBPCG,
	     bool is_ortho = true,
	     const double& tol_AMG_Precond = 1.0e-3,
             const int& max_AMG_Precond = 10,
	     const SparseMatrix<double>& P = SparseMatrix<double>());
  
  void getProduct();


public:
  std::ofstream Record; 
  void open_log(std::string& fileName){Record.open(fileName);};
  void close_log(){Record.close();};

  int iter;
  
  std::ofstream timeinRR;
};

#endif//_LOBPCG_H_
/**
 * end of file
 * 
 */
