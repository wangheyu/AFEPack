/**
 * @file   LOBPCG.cpp
 * @author Hu Guanghui <gary@Brin>
 * @date   Fri Nov 18 23:07:14 2016
 * 
 * @brief  
 * 
 * 
 */
#include "AFEPack/LOBPCG.h"
#include <assert.h>

LOBPCG::LOBPCG()
{
  lapackA = new DenseMatrix<double>();
    lapackB= new DenseMatrix<double>();
    lapackTMP= new DenseMatrix<double>();
    lapackTMP1 = new DenseMatrix<double>();
    lapackTMP2 = new DenseMatrix<double>();
    lapackTMP3 = new DenseMatrix<double>();
    lapackTMP4 = new DenseMatrix<double>();
    lapackTMP5 = new DenseMatrix<double>();
    lapackTMP6 = new DenseMatrix<double>();
    lapackA0 = new DenseMatrix<double>();
    lapackB0 = new DenseMatrix<double>();
}

LOBPCG::~LOBPCG()
{
  delete
    lapackA,
    lapackB,
    lapackTMP,
    lapackTMP1,
    lapackTMP2,
    lapackTMP3,
    lapackTMP4,
    lapackTMP5,
    lapackTMP6,
    lapackA0,
    lapackB0;
}

void LOBPCG::getProduct()
{
  /// suppose that AVec and BVec have been initialized
  const int& n_dof = A->m();
  AVec.resize(n_eigenpair, Vector<double>(n_dof));
  BVec.resize(n_eigenpair, Vector<double>(n_dof));

#pragma omp parallel for
  for(int i = 0;i < n_eigenpair;++ i){
    A->vmult(AVec[i], (*eigenVector[i]));
    B->vmult(BVec[i], (*eigenVector[i]));
  }
}

void LOBPCG::reinit(SparseMatrix<double>& _A, 
		    SparseMatrix<double>& _B, 
		    std::vector<Vector<double> >& _KSOrbital, 
		    std::vector<double>& _KSEigenValue)
{
  timeinRR.open("RRtime.txt", std::ios_base::app);
  
  n_eigenpair = _KSOrbital.size();
  A = &(_A);
  B = &(_B);

  //stiff_matrix = &(_stiff_matrix);
  /// since 0.5 * stiff_matrix is used in the preconditioning process, we do the following here
  //(*stiff_matrix) *= 0.5;
  
  eigenVector.resize(n_eigenpair);
  for(int i = 0;i < n_eigenpair;++ i){
    eigenVector[i] = &(_KSOrbital[i]);    
  }
  eigenValue = &(_KSEigenValue);

  const int& n_dof = A->m();
  searchDirection.clear();
  searchDirection.resize(n_eigenpair, Vector<double>(n_dof));
  residual.clear();
  residual.resize(n_eigenpair, Vector<double>(n_dof));
  precondResidual.clear();
  precondResidual.resize(n_eigenpair, Vector<double>(n_dof));
  AVec.clear();
  AVec.resize(n_eigenpair, Vector<double>(n_dof));
  BVec.clear();
  BVec.resize(n_eigenpair, Vector<double>(n_dof));

  tmpVec.resize(3*n_eigenpair);

  lapackA = new DenseMatrix<double>(3*n_eigenpair, 3*n_eigenpair);
  lapackB = new DenseMatrix<double>(3*n_eigenpair, 3*n_eigenpair);
  lapackTMP = new DenseMatrix<double>(n_eigenpair);
  lapackTMP1 = new DenseMatrix<double>(n_eigenpair);
  lapackTMP2 = new DenseMatrix<double>(n_eigenpair);
  lapackTMP3 = new DenseMatrix<double>(n_eigenpair);
  lapackTMP4 = new DenseMatrix<double>(n_eigenpair);
  lapackTMP5 = new DenseMatrix<double>(n_eigenpair);
  lapackTMP6 = new DenseMatrix<double>(n_eigenpair);
  lapackA0 = new DenseMatrix<double>(2*n_eigenpair, 2*n_eigenpair);
  lapackB0 = new DenseMatrix<double>(2*n_eigenpair, 2*n_eigenpair);

}


void LOBPCG::solve(const double& tol_LOBPCG,
		   const int& max_LOBPCG,
		   bool is_ortho,
		   const double& tol_AMG_Precond,
		   const int& max_AMG_Precond,
		   const SparseMatrix<double>& P)
{ 
  double TimeStamp; // time stamp for recording the computing time.
  TimeStamp = omp_get_wtime();
  
  iter = 0;
  //int iter = 0;
  bool stop = true;
  //SparseMatrix<double> D(stiff_matrix->get_sparsity_pattern());
  
  getProduct();

  std::vector<Vector<double> > tmpVec;
  tmpVec.resize(n_eigenpair);
  for(int i = 0;i < n_eigenpair;++ i){
    tmpVec[i] = *eigenVector[i];
  }
#pragma omp  parallel for 
    for(int i = 0;i < n_eigenpair;++ i){
      residual[i].equ(-(*eigenValue)[i], BVec[i], 1., AVec[i]);
    }
    
#pragma omp parallel for
  for(int i = 0;i < n_eigenpair;++ i){
    (*eigenValue)[i] = (*eigenVector[i]) * AVec[i]/((*eigenVector[i]) * BVec[i]);
  }

  for(int i = 0;i < n_eigenpair;++ i)
    std::cout << std::setprecision(10) << "The " << i << " eigenvalue is " << eigenValue->at(i)
	      << ", and the residual's norm is " << residual[i].l2_norm()
	      << std::endl;

  stop = true;
  double max_residual = 0.;
  for(int i = 0;i < n_eigenpair;++ i){
    if(residual[i].l2_norm() > tol_LOBPCG)
      stop = false;
    if(residual[i].l2_norm() > max_residual)
      max_residual = residual[i].l2_norm();
  }

  double max_orthogonal_error = 0.;

  for(int i = 0;i < n_eigenpair;++ i)
    {
      for(int j = i + 1;j < n_eigenpair;j ++)
	{
	  double tmp_ortho_error; 
	  tmp_ortho_error = (*eigenVector[j]) * BVec[i];
	  //Record << tmp_ortho_error << " ";
	  if(fabs(tmp_ortho_error) > max_orthogonal_error)
	    {
	      max_orthogonal_error =  fabs(tmp_ortho_error);
	    }
	}
    }

 // Record << std::endl;
  Record << 0. << " ";
  Record << 0. << " ";
  Record << 0. << " ";
  Record << iter << " ";//Record the LOBPCG iteration number.
  Record << max_residual << " ";
  Record << max_orthogonal_error <<" ";
  Record << omp_get_wtime() - TimeStamp << std::endl; // Recording the time for calculation one time.
  if(stop) return;
  
  double time4precond;
  if(!P.is_initialize()) // Precondition matrix is not given;
    {
      std::cout<<"P.is_initalize():: " << P.is_initialize() << std::endl;
      std::cout<<"The precondition matrix is not given. The diagonal preconditioner is used as default." << std::endl;

      time4precond = omp_get_wtime();
#pragma omp parallel
      {
	//SparseMatrix<double> D(A->get_sparsity_pattern());
	//D.copy_from(*A);/// be noted that stiff_matrix = 0.5 * stiff_matrix, which has been done in contructor.

#pragma omp for 
	for(int i = 0;i < n_eigenpair;++ i){
	  for(int j = 0;j < A->m();++ j){
	    precondResidual[i](j) = residual[i](j)/(A->diag_element(j) - eigenValue->at(i)*B->diag_element(j));
	  }	  
	}
      }
    }
  else
    {
      std::cout<<"P.is_initalize():: " << P.is_initialize() << std::endl;
      std::cout<<"The Precondition matrix is given..." << std::endl;
      double time4precond;
      time4precond = omp_get_wtime();
#pragma omp parallel
      {
	//SparseMatrix<double> D(A->get_sparsity_pattern());
	//D.copy_from(*A);/// be noted that stiff_matrix = 0.5 * stiff_matrix, which has been done in contructor.
	AMGSolver mmy_amg;
	mmy_amg.lazyReinit(P);//initialize the AMGSolver for the precondition matrix P.
	    
#pragma omp for 
	for(int i = 0;i < n_eigenpair;++ i){
	  if(eigenValue->at(i) > 0.){/// when eigenvalue is positive,
				   /// the positive definite of the
				   /// preconditioner can not
				   /// guaranteed. Hence we use a
				   /// simple one here.
	    for(int j = 0;j < A->m();++ j){
	      precondResidual[i](j) = residual[i](j)/(A->diag_element(j) - eigenValue->at(i)*B->diag_element(j));
	    }
	  }
	  else{
	    //mmy_amg.lazyReinit(P);
	    mmy_amg.solve(precondResidual[i], residual[i], tol_AMG_Precond, max_AMG_Precond);
	  }	  
	}
      }
    }
  Record << omp_get_wtime() - time4precond << " ";
  std::cout << "Begin the first step RRProcedure" << std::endl;
  double time4RR;
  time4RR = omp_get_wtime();
  RRProcedure(RRVector, eigenVector, precondResidual, *A, *B);
  Record << omp_get_wtime() - time4RR << " ";
  std::cout << "The first step RRProcedure is done..." << std::endl;

#pragma omp parallel for
  for(int i = 0;i < n_eigenpair;++ i){
    eigenVector[i]->equ(1., (*eigenVector[i]),  RRVector[i](n_eigenpair + i)/RRVector[i](i), precondResidual[i]);    
  }
  double time4ortho;
  time4ortho = omp_get_wtime();
  if(is_ortho){
    Borthonormalize(*B, eigenVector);
  }
  Record << omp_get_wtime() - time4ortho << " ";
  
  //normalize(eigenVector); // avoid this part for testing.
  
#pragma omp parallel for
  for(int i = 0;i < n_eigenpair;++ i){
    searchDirection[i].equ(RRVector[i](n_eigenpair + i)/RRVector[i](i), precondResidual[i]);
  }

  /////// Start the iterations of LOBPCG ////////
  do{
    getProduct();
    std::cout << "update Product done..." << std::endl;
#pragma omp parallel for
    for(int i = 0;i < n_eigenpair;++ i){
      tmpVec[i] = *eigenVector[i];
    }
 
#pragma omp  parallel for 
    for(int i = 0;i < n_eigenpair;++ i){
      residual[i].equ(-(*eigenValue)[i], BVec[i], 1., AVec[i]);
    }

    double max_orthogonal_error = 0.;

    for(int i = 0;i < n_eigenpair;++ i)
      {
	for(int j = i+1;j < n_eigenpair;j ++)
	  {
	    double tmp_ortho_error; 
	    tmp_ortho_error = (*eigenVector[j]) * BVec[i];
	    //Record << tmp_ortho_error << " ";
	    if(fabs(tmp_ortho_error) > max_orthogonal_error)
	      {
		max_orthogonal_error =  fabs(tmp_ortho_error) ;
	      }
	  }
      }

    //Record << iter + 1 << " ";//Record the LOBPCG iteration number.
    //Record << max_orthogonal_error <<" ";
      
    stop = true;
    max_residual = 0.;
#pragma omp parallel for
    for(int i = 0;i < n_eigenpair;++ i){
      if(residual[i].l2_norm() > tol_LOBPCG)
	{
	  stop = false;
	  //Record << i << " " << residual[i].l2_norm() << " ";
	}
    }
    for(int i = 0;i < n_eigenpair;++ i){
      if(residual[i].l2_norm() > max_residual)
        max_residual = residual[i].l2_norm();
      std::cout << "The " << iter << " iteration: ";
      std::cout << std::setprecision(10) << "The " << i << " eigenvalue is " << eigenValue->at(i)
		<< ", and the residual's norm is " << residual[i].l2_norm()
		<< std::endl;
    }

//    Record << std::endl;
    Record << iter + 1 << " ";//Record the LOBPCG iteration number.
    Record << max_residual << " ";
    Record << max_orthogonal_error <<" ";
    Record << omp_get_wtime() - TimeStamp << std::endl; // Recording the time for calculation one time.
    if(stop) return;

    if(!P.is_initialize()) // Precondition matrix is not given;
      {
	std::cout<<"The precondition matrix is not given. The diagonal preconditioner is used as default." << std::endl;
	time4precond = omp_get_wtime();
#pragma omp parallel
	{      
	  ///otherwise, introduce the preconditioner
	  //SparseMatrix<double> D(A->get_sparsity_pattern());
	  //D.copy_from(*A);/// be noted that stiff_matrix = 0.5 * stiff_matrix, which has been done in contructor.
#pragma omp for 
	  for(int i = 0;i < n_eigenpair;++ i){
	    for(int j = 0;j < A->m();++ j){
	      precondResidual[i](j) = residual[i](j)/(A->diag_element(j) - eigenValue->at(i)*B->diag_element(j));
	    }
	  }
	}
      }
    else
      {
	std::cout<<"The Precondition matrix is given..." << std::endl;
	double time4precond;
	time4precond = omp_get_wtime();
#pragma omp parallel
	{
	  //SparseMatrix<double> D(A->get_sparsity_pattern());
	  //D.copy_from(*A);/// be noted that stiff_matrix = 0.5 * stiff_matrix, which has been done in contructor.
	  AMGSolver mmy_amg;
	  mmy_amg.lazyReinit(P);//initialize the AMGSolver for the precondition matrix P.
	    
#pragma omp for 
	  for(int i = 0;i < n_eigenpair;++ i){
	    if(eigenValue->at(i) > 0.){/// when eigenvalue is positive,
	      /// the positive definite of the
	      /// preconditioner can not
	      /// guaranteed. Hence we use a
	      /// simple one here.
	      for(int j = 0;j < A->m();++ j){
		precondResidual[i](j) = residual[i](j)/(A->diag_element(j) - eigenValue->at(i)*B->diag_element(j));
	      }
	    }
	    else{
	      //mmy_amg.lazyReinit(P);
	      mmy_amg.solve(precondResidual[i], residual[i], tol_AMG_Precond, max_AMG_Precond);
	    }	  
	  }
	}
      }
    if(iter < max_LOBPCG){
    Record <<  omp_get_wtime() - time4precond << " "; // Recording the time for calculation one time.
    }
    std::cout << "preconditioning is done..." << std::endl;
    ///Rayleigh-Ritz procedure
    time4RR = omp_get_wtime();
    RRProcedure(RRVector, eigenVector, precondResidual, searchDirection, *A, *B);
    if(iter < max_LOBPCG){
    Record <<  omp_get_wtime() - time4RR << " "; // Recording the time for calculation one time.
    }
    std::cout << "RRProcedure is done..." << std::endl;

    ///update eigenvector
#pragma omp parallel for
    for(int i = 0;i < n_eigenpair;++ i){
      *eigenVector[i] = 0.;
      for(int j = 0;j < n_eigenpair;++ j){
	eigenVector[i]->add(RRVector[i](0 * n_eigenpair + j), tmpVec[j]);
	eigenVector[i]->add(RRVector[i](1 * n_eigenpair + j), precondResidual[j]);
	eigenVector[i]->add(RRVector[i](2 * n_eigenpair + j), searchDirection[j]);
      }
    }
    /// Borthonormalize eigenVector    
    time4ortho = omp_get_wtime();
    if(is_ortho){
      Borthonormalize(*B, eigenVector);
    }
    if(iter < max_LOBPCG){
      Record << omp_get_wtime() -time4ortho << " ";
    }
    //std::cout << "Borthonomarlization is done..." << std::endl;

    
    //normalize(eigenVector); // avoid the normalize for testing 
        
    std::cout << "nomarlization is done..." << std::endl;


    
    ///update searchDirction
#pragma omp parallel for
    for(int i = 0;i < n_eigenpair;++ i){
      searchDirection[i].equ(RRVector[i](n_eigenpair + i)/RRVector[i](i), precondResidual[i]);
    }
        
  }while(iter++ < max_LOBPCG && !stop);
}

void LOBPCG::normalize(std::vector<Vector<double>* > & v_h)
{
  const int n_size = v_h.size();
  double l2norm = 0.;

#pragma omp parallel for
  for(int i = 0;i < n_size;++ i){
    l2norm =v_h[i]->l2_norm();
    (*v_h[i]) *= 1./l2norm;
  }
}

/* // Original version.
void LOBPCG::Borthonormalize(SparseMatrix<double>& B, 
			     std::vector<Vector<double>* > & v_h)
{
  const int n_size = v_h.size();
  double l2norm = 0.;
  l2norm = v_h[0]->l2_norm();
  //  l2norm = Functional::L2Norm(v_h[0], algebraicAccuracy);
  (*v_h[0]) *= 1./l2norm;
  Vector<double> tmp(v_h[0]->size());

  for(int i = 1;i < n_size;++ i){
    /// for orthogonalizaion
    for(int j = 0;j < i;++ j){
      B.vmult(tmp, (*v_h[i]));
      l2norm = (*v_h[j]) * tmp;
      v_h[i]->add(-l2norm, (*v_h[j]));
    }
    
    l2norm =v_h[i]->l2_norm();
    //    l2norm = Functional::L2Norm(v_h[i], algebraicAccuracy);
    //l2norm = v_h[0].l2_norm();
    (*v_h[i]) *= 1./l2norm;
  }
}
*/
void LOBPCG::Borthonormalize(SparseMatrix<double>& B,
                             std::vector<Vector<double>* > & v_h)
{
  const int n_size = v_h.size();
  double l2norm = 0.;
  // l2norm = v_h[0]->l2_norm();
  //  l2norm = Functional::L2Norm(v_h[0], algebraicAccuracy);
  //(*v_h[0]) *= 1./l2norm;
  Vector<double> tmp(v_h[0]->size());
  B.vmult(tmp, (*v_h[0]));
  l2norm = (*v_h[0]) * tmp;

  (*v_h[0]) *= 1./l2norm;

  for(int i = 1;i < n_size;++ i){
    /// for orthogonalizaion
    for(int j = 0;j < i;++ j){
      B.vmult(tmp, (*v_h[i]));
      l2norm = (*v_h[j]) * tmp;
      v_h[i]->add(-l2norm, (*v_h[j]));
    }

    B.vmult(tmp, (*v_h[i]));
    l2norm = (*v_h[i]) * tmp;

    (*v_h[i]) *= 1./l2norm;
  }
}


void LOBPCG::getPMQt(DenseMatrix<double>& fm,
		     const std::vector<Vector<double>* >& vec,
		     const SparseMatrix<double>& sm,
		     const std::vector<Vector<double> >& vect)
{/// fm = vec * sm * vect
  const int& vec_size = vec.size();
  const int& vect_size = vect.size();
  assert(vec_size == vect_size);
  assert(vec[0]->size() == vect[0].size());
  assert(vec[0]->size() == sm.m());

  fm.reinit(vec_size, vec_size);
#pragma omp parallel for
  for(int i = 0;i < vect_size;++ i){
      Vector<double> tmp(vec[0]->size());
    sm.vmult(tmp, vect[i]);
    for(int j = 0;j < vec_size;++ j){
      fm(j, i) = (*vec[j]) * tmp;
    }    
  }
}

void LOBPCG::getPMQt(DenseMatrix<double>& fm,
		     const std::vector<Vector<double> >& vec,
		     const SparseMatrix<double>& sm,
		     const std::vector<Vector<double>* >& vect)
{/// fm = vec * sm * vect
  const int& vec_size = vec.size();
  const int& vect_size = vect.size();
  assert(vec_size == vect_size);
  assert(vec[0].size() == vect[0]->size());
  assert(vec[0].size() == sm.m());

  fm.reinit(vec_size, vec_size);
#pragma omp parallel for
  for(int i = 0;i < vect_size;++ i){
  Vector<double> tmp(vec[0].size());
    sm.vmult(tmp, (*vect[i]));
    for(int j = 0;j < vec_size;++ j){
      fm(j, i) = vec[j] * tmp;
    }    
  }
}

void LOBPCG::getPMQt(DenseMatrix<double>& fm,
		     const std::vector<Vector<double>* >& vec,
		     const SparseMatrix<double>& sm,
		     const std::vector<Vector<double>* >& vect)
{/// fm = vec * sm * vect
  const int& vec_size = vec.size();
  const int& vect_size = vect.size();
  assert(vec_size == vect_size);
  assert(vec[0]->size() == vect[0]->size());
  assert(vec[0]->size() == sm.m());

  fm.reinit(vec_size, vec_size);
#pragma omp parallel for
  for(int i = 0;i < vect_size;++ i){
  Vector<double> tmp(vec[0]->size());
    sm.vmult(tmp, (*vect[i]));
    for(int j = 0;j < vec_size;++ j){
      fm(j, i) = (*vec[j]) * tmp;
    }    
  }
}


void LOBPCG::getPMQt(DenseMatrix<double>& fm,
		     const std::vector<Vector<double> >& vec,
		     const SparseMatrix<double>& sm,
		     const std::vector<Vector<double> >& vect)
{/// fm = vec * sm * vect
  const int& vec_size = vec.size();
  const int& vect_size = vect.size();
  assert(vec_size == vect_size);
  assert(vec[0].size() == vect[0].size());
  assert(vec[0].size() == sm.m());

  fm.reinit(vec_size, vec_size);
#pragma omp parallel for
  for(int i = 0;i < vect_size;++ i){
  Vector<double> tmp(vec[0].size());
    sm.vmult(tmp, vect[i]);
    for(int j = 0;j < vec_size;++ j){
      fm(j, i) = vec[j] * tmp;
    }    
  }
}

void getvmult(DenseMatrix<double>& A,
	      const Vector<double>& vec,
	      Vector<double>& Avec)
{
  const int& vec_size = vec.size();
  Avec.reinit(vec.size());
  // Vector<double> tmp(vec_size);
  for(int i = 0;i < vec_size;i ++)
    {
      Avec[i] = 0.;
      for(int j = 0;j < vec_size;j ++)
	{
	  Avec[i] += A(i,j) * vec[j];
	}
    }
}

void LOBPCG::catMatrix(DenseMatrix<double>& dst,
		       const DenseMatrix<double>& src,
		       const int& blockSize,/// 3
		       const int& row_idx,/// 1
		       const int& col_idx)///2
{
  //  assert(dst.n_rows() == blockSize * src.n_rows());
  const int& n_size = src.get_n();

#pragma omp parallel for
  for(int i = 0;i < n_size;++ i){
    for(int j = 0;j < n_size;++ j){
      dst(row_idx * blockSize + i, col_idx * blockSize + j) = src(i, j);
    }
  }
}



/** 
 * 1. form the small eigenvalue system
 * 2. solve it
 * 3. Ax = \lambda Bx
 */
void LOBPCG::RRProcedure(std::vector<Vector<double> >& RRVector,
			 const std::vector<Vector<double>* >& eigenVector,
			 const std::vector<Vector<double> >& residual,
			 const SparseMatrix<double>& A,
			 const SparseMatrix<double>& B)
{
  assert(A.m() == A.n());
  assert(A.m() == B.m());
  assert(B.m() == B.n());
  assert(A.m() == eigenVector[0]->size());
  assert(eigenVector.size() == residual.size());
  assert(eigenVector[0]->size() == residual[0].size());

  ///reset entries
  (*lapackA0) = 0.;
  (*lapackB0) = 0.;
  (*lapackTMP) = 0.;
  
  // LAPACKFullMatrix<double>
  //   lapackA(2*n_eigenpair, 2*n_eigenpair),
  //   lapackB(2*n_eigenpair, 2*n_eigenpair),
  //   lapackTMP(n_eigenpair);
  
  /// puzzling
  double time4puzzling;
  time4puzzling = omp_get_wtime();
  
  getPMQt(*lapackTMP, eigenVector, A, eigenVector);
  catMatrix(*lapackA0, *lapackTMP, n_eigenpair, 0, 0);
  getPMQt(*lapackTMP, eigenVector, A, residual);
  catMatrix(*lapackA0, *lapackTMP, n_eigenpair, 0, 1);
  getPMQt(*lapackTMP, residual, A, eigenVector);
  catMatrix(*lapackA0, *lapackTMP, n_eigenpair, 1, 0);
  getPMQt(*lapackTMP, residual, A, residual);
  catMatrix(*lapackA0, *lapackTMP, n_eigenpair, 1, 1);

  getPMQt(*lapackTMP, eigenVector, B, eigenVector);
  catMatrix(*lapackB0, *lapackTMP, n_eigenpair, 0, 0);
  getPMQt(*lapackTMP, eigenVector, B, residual);
  catMatrix(*lapackB0, *lapackTMP, n_eigenpair, 0, 1);
  getPMQt(*lapackTMP, residual, B, eigenVector);
  catMatrix(*lapackB0, *lapackTMP, n_eigenpair, 1, 0);
  getPMQt(*lapackTMP, residual, B, residual);
  catMatrix(*lapackB0, *lapackTMP, n_eigenpair, 1, 1);
  
  /////////////////////////////////////////////////////////
  time4puzzling = omp_get_wtime() - time4puzzling;
  
  // lapackA.print_formatted(std::cout);
  // std::cout << "---------------------------------------------------" << std::endl;
  // lapackB.print_formatted(std::cout);

  //  getchar();
  /// now we have lapackA x = \lambda lapackB x, solve it
  // Here we can compute the condition number of matrix lapackB.
  // With our designed preconditioning, the condition number of lapackB can be controlled.
  double time4RRVecResize;
  time4RRVecResize = omp_get_wtime();
  RRVector.clear();
  RRVector.resize(2*n_eigenpair, Vector<double>(2*n_eigenpair));
  time4RRVecResize = omp_get_wtime() - time4RRVecResize;
  double time4DSYGV;
  time4DSYGV = omp_get_wtime();
  lapackA0->compute_generalized_eigenvalues_symmetric(*lapackB0, RRVector, 1);
  time4DSYGV = omp_get_wtime() - time4DSYGV;
  for(int i = 0;i < n_eigenpair;++ i){
    (*eigenValue)[i] = lapackA0->eigenvalue(i).real();
  }
  /// we are interested in the first n_eigenpair eigenvectors
  
  timeinRR << iter << " " << time4puzzling << " " << time4RRVecResize << " " << time4DSYGV << std::endl;
  
}

/** 
 * 1. form the small eigenvalue system
 * 2. solve it
 * 3. Ax = \lambda Bx
 */
void LOBPCG::RRProcedure(std::vector<Vector<double> >& RRVector,
			 const std::vector<Vector<double>* >& eigenVector,
			 const std::vector<Vector<double> >& residual,
			 const std::vector<Vector<double> >& searchDirection,
			 const SparseMatrix<double>& A,
			 const SparseMatrix<double>& B)
{
  assert(A.m() == A.n());
  assert(A.m() == B.m());
  assert(B.m() == B.n());
  assert(A.m() == eigenVector[0]->size());
  assert(eigenVector.size() == residual.size());
  assert(eigenVector.size() == searchDirection.size());
  assert(eigenVector[0]->size() == residual[0].size());


  (*lapackA) = 0.;
  (*lapackB) = 0.;
  (*lapackTMP1) = 0.;
  (*lapackTMP2) = 0.;
  (*lapackTMP3) = 0.;
  (*lapackTMP4) = 0.;
  (*lapackTMP5) = 0.;
  (*lapackTMP6) = 0.;
  // LAPACKFullMatrix<double>
  //   lapackA(3*n_eigenpair, 3*n_eigenpair),
  //   lapackB(3*n_eigenpair, 3*n_eigenpair),
  //   lapackTMP1(n_eigenpair),
  //   lapackTMP2(n_eigenpair),
  //   lapackTMP3(n_eigenpair),
  //   lapackTMP4(n_eigenpair),
  //   lapackTMP5(n_eigenpair),
  //   lapackTMP6(n_eigenpair);
    //lapacktA(3*n_eigenpair, 3*n_eigenpair),
    //lapacktB(3*n_eigenpair, 3*n_eigenpair);


  /// puzzling
  double time4puzzling;
  time4puzzling = omp_get_wtime();
  
#pragma omp parallel sections
  {
#pragma omp section
    {
      getPMQt(*lapackTMP1, eigenVector, A, eigenVector);
      catMatrix(*lapackA, *lapackTMP1, n_eigenpair, 0, 0);
      getPMQt(*lapackTMP1, eigenVector, A, residual);
      catMatrix(*lapackA, *lapackTMP1, n_eigenpair, 0, 1);
      getPMQt(*lapackTMP1, eigenVector, A, searchDirection);
      catMatrix(*lapackA, *lapackTMP1, n_eigenpair, 0, 2);
    }
#pragma omp section    
    {
      getPMQt(*lapackTMP2, residual, A, eigenVector);
      catMatrix(*lapackA, *lapackTMP2, n_eigenpair, 1, 0);
      getPMQt(*lapackTMP2, residual, A, residual);
      catMatrix(*lapackA, *lapackTMP2, n_eigenpair, 1, 1);
      getPMQt(*lapackTMP2, residual, A, searchDirection);
      catMatrix(*lapackA, *lapackTMP2, n_eigenpair, 1, 2);
    }
#pragma omp section
    {
      getPMQt(*lapackTMP3, searchDirection, A, eigenVector);
      catMatrix(*lapackA, *lapackTMP3, n_eigenpair, 2, 0);
      getPMQt(*lapackTMP3, searchDirection, A, residual);
      catMatrix(*lapackA, *lapackTMP3, n_eigenpair, 2, 1);
      getPMQt(*lapackTMP3, searchDirection, A, searchDirection);
      catMatrix(*lapackA, *lapackTMP3, n_eigenpair, 2, 2);
    }
#pragma omp section  
    {
      getPMQt(*lapackTMP4, eigenVector, B, eigenVector);
      catMatrix(*lapackB, *lapackTMP4, n_eigenpair, 0, 0);
      getPMQt(*lapackTMP4, eigenVector, B, residual);
      catMatrix(*lapackB, *lapackTMP4, n_eigenpair, 0, 1);
      getPMQt(*lapackTMP4, eigenVector, B, searchDirection);
      catMatrix(*lapackB, *lapackTMP4, n_eigenpair, 0, 2);
    }
#pragma omp section
    {
      getPMQt(*lapackTMP5, residual, B, eigenVector);
      catMatrix(*lapackB, *lapackTMP5, n_eigenpair, 1, 0);
      getPMQt(*lapackTMP5, residual, B, residual);
      catMatrix(*lapackB, *lapackTMP5, n_eigenpair, 1, 1);
      getPMQt(*lapackTMP5, residual, B, searchDirection);
      catMatrix(*lapackB, *lapackTMP5, n_eigenpair, 1, 2);
    }
#pragma omp section
    {
      getPMQt(*lapackTMP6, searchDirection, B, eigenVector);
      catMatrix(*lapackB, *lapackTMP6, n_eigenpair, 2, 0);
      getPMQt(*lapackTMP6, searchDirection, B, residual);
      catMatrix(*lapackB, *lapackTMP6, n_eigenpair, 2, 1);
      getPMQt(*lapackTMP6, searchDirection, B, searchDirection);
      catMatrix(*lapackB, *lapackTMP6, n_eigenpair, 2, 2);
    }
  }
  ///////////////////////////////////////////////
  time4puzzling = omp_get_wtime() - time4puzzling;
  
 /* 
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout<<"lapackA and B before the generalized eigenvalues function:\n";
  lapackA->print_formatted(std::cout);
  std::cout << "-------------------------------------------------" << std::endl;
  lapackB->print_formatted(std::cout);
  std::cout << "-------------------------------------------------" << std::endl;

  LAPACKFullMatrix<double> lapackAt, lapackBt;
  lapackAt = *lapackA;
  lapackBt = *lapackB;
*/
  //  lapackA.print_formatted(std::cout);
 
  /// now we have lapackA x = \lambda lapackB x, solve it
  double time4RRVecResize;
  time4RRVecResize = omp_get_wtime();
  
  RRVector.clear();
  RRVector.resize(3*n_eigenpair, Vector<double>(3*n_eigenpair));

  time4RRVecResize = omp_get_wtime() - time4RRVecResize;
  double time4DSYGV;
  time4DSYGV = omp_get_wtime();

  lapackA->compute_generalized_eigenvalues_symmetric(*lapackB, RRVector, 1); // This dealii function uses the subroutine DSYGVX().
  
  time4DSYGV = omp_get_wtime() - time4DSYGV;
#pragma omp parallel for
  for(int i = 0;i < n_eigenpair;++ i){
    (*eigenValue)[i] = lapackA->eigenvalue(i).real();
  }

  /*
   * Here is the computation of the condition number of the Gram matrix B.
   */
  timeinRR << iter << " " << time4puzzling << " " << time4RRVecResize << " " << time4DSYGV << std::endl;
}


/**
 * end of file
 * 
 */
