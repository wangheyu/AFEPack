#include <AFEPack/DenseMatrix.h>

AFEPACK_OPEN_NAMESPACE

template <typename T>
DenseMatrix<T>::DenseMatrix(const unsigned int& _m,
			 const unsigned int& _n,
			 const T& _entry)
{
  m = _m;
  n = _n;

  std::vector<T>::resize(m * n, _entry);
  if(m == n) is_square = true;
}

template <typename T>
DenseMatrix<T>::DenseMatrix(const unsigned int& _m)
{
  m = _m;
  n = _m;

  std::vector<T>::resize(m * n, 0.);
  if(m == n) is_square = true;
}


template <typename T>
void DenseMatrix<T>::reinit_impl(const unsigned int& _m,
			 const unsigned int& _n)/*,
			 //const T& _entry)*/
{
  m = _m;
  n = _n;

  //std::vector<T>::resize(m * n, _entry);
  std::vector<T>::clear();
  std::vector<T>::resize(m * n); // updated Vecma;
  if(m == n) is_square = true;  
}

template <typename T>
void DenseMatrix<T>::set_entry_impl(const unsigned int& row_idx,
		      const unsigned int& col_idx,
		      const T& entry)
{
  std::vector<T>::operator[](row_idx * n + col_idx) = entry;
}

template <typename T>
T * DenseMatrix<T>::get_data_impl()
{
  return std::vector<T>::data();
}

template <typename T>
const T * DenseMatrix<T>::get_data_impl() const
{
  return std::vector<T>::data();
}

template <typename T>
void DenseMatrix<T>::add_entry_impl(const unsigned int& row_idx,
		      const unsigned int& col_idx,
		      const T& entry)
{
  std::vector<T>::operator[](row_idx * n + col_idx) += entry;
}


/// des = this + des
template <typename T>
void DenseMatrix<T>::add_m_impl(DenseMatrix<T>& des)
{
  assert(m == des.get_m_impl());
  assert(n == des.get_n_impl());

  unsigned int sizeOfData = this->size();

  cblas_daxpy(sizeOfData, 1., get_data_impl(), 1, des.get_data_impl(), 1);
}



/// SparseMatrix(des) = DenseMatrix(this) + SparseMatrix(des).
/// The output is a sparse matrix. However, since
/// (this) is a dense matrix, it is quite reasonable to output a dense matrix too.
/// So 1). we transfer des from sparse to dense, and 2). call above dense + dense subroutine.
template <typename T>
void DenseMatrix<T>::add_m_impl(SparseMatrix<T>& des)
{
  assert(m == des.get_m_impl());
  assert(n == des.get_n_impl());

  unsigned int sizeOfData = this->size();

  Vecma<SparseMatrix<T>, T> * tmp_ptr = static_cast<Vecma<SparseMatrix<T>, T> *> (&des);
  //Vecma<DenseMatrix<T>, T> * dense_ptr = static_cast<Vecma<DenseMatrix<T>, T> *> (tmp_ptr);

  cblas_daxpy(sizeOfData, 1., get_data_impl(), 1, des.get_data_impl(), 1);
}



/// des := this * src
template <typename T>
//void DenseMatrix<T>::mult_v_impl(Vector<T>& src, Vector<T>& des)
void DenseMatrix<T>::mult_v_impl(Vector<T>& des, Vector<T>& src) const
{
  assert(n == src.size());
  assert(des.size() == src.size());
  des.reinit(n);
  
  cblas_dgemv(CblasRowMajor,
	      CblasNoTrans,
	      m,
	      n,
	      1.,/// alpha
	      get_data_impl(),
	      n, //orgin is n
	      src.get_data_impl(),
	      1,
	      0.,/// beta
	      des.get_data_impl(),
	      1);
}

/// des := this^T * src
template <typename T>
//void DenseMatrix<T>::Tmult_v_impl(Vector<T>& src, Vector<T>& des)
void DenseMatrix<T>::Tmult_v_impl(Vector<T>& des, Vector<T>& src)
{
  assert(n == src.size());
  assert(des.size() == src.size());
  
  cblas_dgemv(CblasRowMajor,
	      CblasTrans,
	      m,
	      n,
	      1.,/// alpha
	      get_data_impl(),
	      n,
	      src.get_data_impl(),
	      1,
	      0.,/// beta
	      des.get_data_impl(),
	      1);
}

/// des := alpha * A * src + beta * y,
template <typename T>
void DenseMatrix<T>::mult_v_impl(const T& alpha, Vector<T>& src, const T& beta, Vector<T>& des) const
{
  assert(n == src.size());
  assert(des.size() == src.size());
  
  cblas_dgemv(CblasRowMajor,
	      CblasNoTrans,
	      m,
	      n,
	      alpha,
	      get_data_impl(),
	      n,
	      src.get_data_impl(),
	      1,
	      beta,
	      des.get_data_impl(),
	      1);
}

/// des := alpha * A * src + beta * y,
template <typename T>
void DenseMatrix<T>::Tmult_v_impl(const T& alpha, Vector<T>& src, const T& beta, Vector<T>& des) 
{
  T * des_data = des.get_data_impl();
  T * src_data = src.get_data_impl();
  T * mat_data = std::vector<T>::data();

  cblas_dgemv(CblasRowMajor,
	      CblasTrans,
	      m,
	      n,
	      alpha,
	      mat_data,
	      m,/// leading dimension of srcB
	      src_data,
	      1,
	      beta,
	      des_data,
	      1);
}

/// des := This *src
template <typename T>
void DenseMatrix<T>::mult_m_impl(DenseMatrix<T>& src,
			    DenseMatrix<T>& des)
{
  assert(get_n_impl() == src.get_m_impl());
  assert(src.get_m_impl() == des.get_m_impl());
  assert(get_m_impl() == des.get_m_impl());
  
  cblas_dgemm(CblasRowMajor,
  	      CblasNoTrans,
  	      CblasNoTrans,
  	      get_m_impl(),
  	      src.get_n_impl(),
  	      get_n_impl(),
  	      1.,
  	      get_data_impl(),
	      get_n_impl(),/// leading dimension of srcB
  	      src.get_data_impl(),
  	      src.get_n_impl(),/// leading dimension of srcB
  	      0.,
  	      des.get_data_impl(),
  	      des.get_n_impl()/// leading dimension of srcB
  	      );
}

/// des := This^T *src
template <typename T>
void DenseMatrix<T>::Tmult_m_impl(DenseMatrix<T>& src,
			     DenseMatrix<T>& des)
{
  assert(get_n_impl() == src.get_m_impl());
  assert(src.get_m_impl() == des.get_m_impl());
  assert(get_m_impl() == des.get_m_impl());
  
  cblas_dgemm(CblasRowMajor,
  	      CblasTrans,
  	      CblasNoTrans,
  	      get_m_impl(),
  	      src.get_n_impl(),
  	      get_n_impl(),
  	      1.,
  	      get_data_impl(),
	      get_n_impl(),/// leading dimension of srcB
  	      src.get_data_impl(),
  	      src.get_n_impl(),/// leading dimension of srcB
  	      0.,
  	      des.get_data_impl(),
  	      des.get_n_impl()/// leading dimension of srcB
  	      );
}

/// des := alpha * this * src + beta * des
template <typename T>
void DenseMatrix<T>::mult_m_impl(const T& alpha,
			    DenseMatrix<T>& src,
			    const T& beta,
			    DenseMatrix<T>& des)
{
  assert(get_n_impl() == src.get_m_impl());
  assert(get_m_impl() == des.get_m_impl());
  assert(src.get_n_impl() == des.get_n_impl());

  cblas_dgemm(CblasRowMajor,
  	      CblasNoTrans,
  	      CblasNoTrans,
  	      get_m_impl(),
  	      src.get_n_impl(),
  	      get_n_impl(),
  	      alpha,
  	      get_data_impl(),
  	      get_n_impl(),/// leading dimension of srcB
  	      src.get_data_impl(),
  	      src.get_n_impl(),/// leading dimension of srcB
  	      beta,
  	      des.get_data_impl(),
  	      des.get_n_impl()/// leading dimension of srcB
  	      );
}


/// des := alpha * this^T * src + beta * des
template <typename T>
void DenseMatrix<T>::Tmult_m_impl(const T& alpha,
			     DenseMatrix<T>& src,
			     const T& beta,
			     DenseMatrix<T>& des)
{
  assert(get_n_impl() == src.get_m_impl());
  assert(get_m_impl() == des.get_m_impl());
  assert(src.get_n_impl() == des.get_n_impl());

  cblas_dgemm(CblasRowMajor,
  	      CblasTrans,
  	      CblasNoTrans,
  	      get_m_impl(),
  	      src.get_n_impl(),
  	      get_n_impl(),
  	      alpha,
  	      get_data_impl(),
  	      get_n_impl(),/// leading dimension of srcB
  	      src.get_data_impl(),
  	      src.get_n_impl(),/// leading dimension of srcB
  	      beta,
  	      des.get_data_impl(),
  	      des.get_n_impl()/// leading dimension of srcB
  	      );
}


template <typename T>
void DenseMatrix<T>::print_impl()
{
  for(unsigned int i = 0;i < m;++ i){
    std::cout << "| ";
    for(unsigned int j = 0;j < n;++ j){
      std::cout << std::vector<T>::operator[](i * n + j) << " ";
    }
    std::cout << " |" << std::endl;
  }
}

template <typename T>
const unsigned int& DenseMatrix<T>::get_m_impl() const
{
  return m;
}

template <typename T>
const unsigned int& DenseMatrix<T>::get_n_impl() const
{
  return n;
}

template <typename T>
void DenseMatrix<T>::gauss_jordan()
{
  assert(get_m_impl() == get_n_impl());
  if (get_n_impl() == 2)
    {
      double a = (*this)(0,0);
      double b = (*this)(0,1);
      double c = (*this)(1,0);
      double d = (*this)(1,1);
      assert(a*d != b*c);
      double inv = a*d-b*c;
      (*this)(0,0) = 1.0/inv*d;
      (*this)(0,1) = 1.0/inv*(-b);
      (*this)(1,0) = 1.0/inv*(-c);
      (*this)(1,1) = 1.0/inv*a;
      return;
    }
  // if (get_n_impl() == 3)
  //  {
  //    return;
  //  }
  //if (get_n_impl() > 15)
  else  
  {
      //const int nn = get_n_impl();
      int nn = get_n_impl();
      
      // workspace for permutations
      std::vector<int> ipiv(nn);
      int info;
      
      // Use the LAPACK function getrf for
      // calculating the LU factorization.
      //getrf(&nn, &nn, &this->values[0], &nn, &ipiv[0], &info);
      /********************************************************************/
      // Use the cBLAS function cblas_getrf for
      // calculating the LU factorization.
      //dgetrf(nn, nn, get_data_impl(), get_n_impl(), &ipiv[0], &info);
      dgetrf_(&nn, &nn, get_data_impl(), &nn, &ipiv[0], &info);
      
	// Assert(info >= 0, ExcInternalError());
	//Assert(info == 0, LACExceptions::ExcSingular());
      
      // scratch array
      std::vector<T> inv_work (nn);
      
      // Use the LAPACK function getri for
      // calculating the actual inverse using
      // the LU factorization.
      //getri(&nn, &this->values[0], &nn, &ipiv[0], &inv_work[0], &nn, &info);
      /********************************************************************/
      // Use the cBLAS function cblas_getrf for
      // calculating the actual inverse using
      // the LU factorization.
      dgetri_(&nn, get_data_impl(), &nn, &ipiv[0], &inv_work[0], &nn, &info);
	
	//Assert(info >= 0, ExcInternalError());
	//Assert(info == 0, LACExceptions::ExcSingular());
      return;
    }
}

extern "C"
{
  void dsygv_ (const int *itype, const char *jobz, const char *uplo,
               const int *n, double *A, const int *lda, double *B,
               const int *ldb, double *w, double *work,
               const int *lwork, int *info);
}
template <typename T>
void DenseMatrix<T>::compute_generalized_eigenvalues_symmetric(DenseMatrix<T> &B, std::vector<Vector<T>> &eigenvectors, int itype)
{
  const int nn = this->get_n();

  wr.resize(nn);
  wi.resize(nn); //This is set purley for consistency reasons with the
		  //eigenvalues() function.
		  //
  T *values_A = const_cast<T *> (this->get_data());
  T *values_B = const_cast<T *> (B.get_data());

  //  int info  = 0;
  int lwork = 1;
  const char jobz = ((eigenvectors.size() > 0) ? ('V') : ('N'));
  const char uplo = 'U';
  int info;
  lwork = -1;
  work.resize(1);
//  std::cout<<"The Matrix A is: " << std::endl;
//  this->print();
//  std::cout<<"The Matrix B is: " << std::endl;
//  B.print();
  dsygv_ (&itype, &jobz, &uplo, &nn, values_A, &nn, values_B, &nn, &wr[0], &work[0], &lwork, &info);
  //int info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 1, jobz, uplo, nn, values_A, nn, values_B, nn, &wr[0]);

//				    std::cout<<"After DSYGV()..." << std::endl;
//  std::cout<<"The Matrix A is: " << std::endl;
//  this->print();
//  std::cout<<"The Matrix B is: " << std::endl;
//    B.print();
  lwork = (int) (work[0]+.1);
  work.resize((size_t) lwork);
  dsygv_ (&itype, &jobz, &uplo, &nn, values_A, &nn, values_B, &nn, &wr[0], &work[0], &lwork, &info);
  //dsygv (itype, jobz, uplo, &nn, values_A, &nn, values_B, &nn, &wr[0], &work[0], lwork, info);

  if(info == 0){
    for (size_t i=0; i < eigenvectors.size(); ++i)
      {
	size_t col_begin(i*nn);
	//eigenvectors[i].reinit(nn, true);
	eigenvectors[i].reinit(nn);
	for (size_t j=0; j < static_cast<size_t>(nn); ++j)	
	  {
	    eigenvectors[i](j) = values_A[col_begin+j];
	  }
      }
  }
  else{
    std::cout << "Failed to compute eigenvalues and eigenvectors." << std::endl;
  }
}
template<typename T> 
std::complex<T> DenseMatrix<T>::eigenvalue (const size_t i) const
{
  return std::complex<T>(wr[i], wi[i]);
}

template<typename T>
DenseMatrix<T> &DenseMatrix<T>::operator = (double val)
{
  std::vector<T>::clear();
  std::vector<T>::resize(m * n, val);
  return *this;
}

template class DenseMatrix<double>;

AFEPACK_CLOSE_NAMESPACE
