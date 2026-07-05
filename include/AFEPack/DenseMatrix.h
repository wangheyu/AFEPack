#ifndef _DENSEMATRIX_H_
#define _DENSEMATRIX_H_

#include <AFEPack/Vector.h>
#include <AFEPack/MatrixBase.h>
#include <cassert>

// DenseMatrix class
template <FloatingPoint T>
class DenseMatrix : public MatrixBase<DenseMatrix<T>, T> {
public://///////////////////////////////////////// constructors and Initilization
  DenseMatrix(size_t rows, size_t cols) : n_rows_(rows), n_cols_(cols), data_(rows * cols, T{}) {}
  DenseMatrix(){};
  void reinit(size_t rows, size_t cols){///all entries become zero
    n_rows_ = rows;
    n_cols_ = cols;
    data_.resize(rows * cols);
    std::fill(data_.begin(), data_.end(), 0.);    
  }

  DenseMatrix& operator=(const DenseMatrix& other)
  {
    if(this == &other) return *this;

    n_rows_ = other.n_rows_;
    n_cols_ = other.n_cols_;
    data_ = other.data_;

    tau_ = other.tau_;
    is_QR_factorized_ = other.is_QR_factorized_;

    return *this;
  }
  


public://///////////////////////////////////////// Element Access and Modification
  T operator()(size_t i, size_t j) const {
    if (i >= n_rows_ || j >= n_cols_) throw std::out_of_range("Invalid indices");
    return data_[i * n_cols_ + j];
  }

  T& operator()(size_t i, size_t j) {
    if (i >= n_rows_ || j >= n_cols_) throw std::out_of_range("Invalid indices");
    return data_[i * n_cols_ + j];
  }
  
public://///////////////////////////////////////// Arithmetic Operators
  DenseMatrix add(const DenseMatrix& other) const {
    if (n_rows() != other.n_rows() || n_cols() != other.n_cols()) {
      throw std::invalid_argument("Matrix dimensions mismatch");
    }
    DenseMatrix result(n_rows(), n_cols());
    if constexpr (std::is_same_v<T, float>) {
      cblas_saxpy(n_rows() * n_cols(), 1.0, data_.data(), 1, other.data_.data(), 1, result.data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_daxpy(n_rows() * n_cols(), 1.0, data_.data(), 1, other.data_.data(), 1, result.data_.data(), 1);
    }
    return result;
  }

  // Method to perform matA = matA + a * matB
  void add(T a, const DenseMatrix& B) {
    if (n_rows_ != B.n_rows_ || n_cols_ != B.n_cols_) {
      throw std::invalid_argument("Matrix dimensions mismatch");
    }

    size_t n = n_rows_ * n_cols_;
    if constexpr (std::is_same_v<T, float>) {
      cblas_saxpy(n, a, B.data_.data(), 1, data_.data(), 1);
    } else if constexpr (std::is_same_v<T, double>) {
      cblas_daxpy(n, a, B.data_.data(), 1, data_.data(), 1);
    } else {
      throw std::invalid_argument("Unsupported data type");
    }

    // Reset QR factorization flag since the matrix has changed
    is_QR_factorized_ = false;
  }

  DenseMatrix& operator*=(const int val)
  {
    if(val){/// if val is not zero
      std::transform(data_.begin(), data_.end(), data_.begin(), [val](T x) { return x * val; });      
    }
    else{/// otherwise, val is zero, so we reset all entries by 0.
      std::fill(data_.begin(), data_.end(), 0.);
    }

    is_QR_factorized_ = false;
    return *this;
  }

  DenseMatrix& operator*=(const T val)
  {
    std::transform(data_.begin(), data_.end(), data_.begin(), [val](T x) { return x * val; });      
    
    is_QR_factorized_ = false;/// since the value has changed

    return *this;
  }

  DenseMatrix multiply(const DenseMatrix& other) const {
    if (n_cols() != other.n_rows()) throw std::invalid_argument("Matrix dimensions mismatch");
    DenseMatrix result(n_rows(), other.n_cols());
    if constexpr (std::is_same_v<T, float>){
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  n_rows(), other.n_cols(), n_cols(), 1.0, data_.data(), n_cols(),
		  other.data_.data(), other.n_cols(), 0.0, result.data_.data(), result.n_cols());
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		  n_rows(), other.n_cols(), n_cols(), 1.0, data_.data(), n_cols(),
		  other.data_.data(), other.n_cols(), 0.0, result.data_.data(), result.n_cols());
    }
    return result;
  }  

  void solve_linear_system(std::vector<T>& b) const {
    if (n_rows() != n_cols() || b.size() != n_rows()) throw std::invalid_argument("Invalid dimensions");
    std::vector<int> ipiv(n_rows());
    if constexpr (std::is_same_v<T, float>){
      LAPACKE_sgesv(LAPACK_ROW_MAJOR, n_rows(), 1, data_.data(), n_cols(), ipiv.data(), b.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      LAPACKE_dgesv(LAPACK_ROW_MAJOR, n_rows(), 1, data_.data(), n_cols(), ipiv.data(), b.data(), 1);
    }
  }  

  void gauss_jordan(){
    /// to create this gauss_jordan just for the compatibility.
    inverse();
  }

  void inverse() {
    if (n_rows() != n_cols()) throw std::invalid_argument("Matrix must be square to compute the inverse");

    int n = static_cast<int>(n_rows_);
    std::vector<int> ipiv(n);

    if constexpr (std::is_same_v<T, float>){
      // Perform LU decomposition
      lapack_int info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, data_.data(), n, ipiv.data());
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgetrf failed with error code: " +
				 std::to_string(info));
      }

      // Compute the inverse using the LU decomposition
      info = LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, data_.data(), n, ipiv.data());
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgetri failed with error code: " +
				 std::to_string(info));
      }
    }
    else if constexpr (std::is_same_v<T, double>){
      // Perform LU decomposition
      lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, data_.data(), n, ipiv.data());
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgetrf failed with error code: " +
				 std::to_string(info));
      }

      // Compute the inverse using the LU decomposition
      info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, data_.data(), n, ipiv.data());
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgetri failed with error code: " +
				 std::to_string(info));
      }
    }
  }

  /*
   * find the inverse of mat, and assign it to *this
   */
  void invert(const DenseMatrix<T>& other)
  {
    *this = other;

    inverse();
  }

  void QRFactorize()
  {/// after this QR factorization is implemented, the tau_ will also be filled.


    lapack_int info;

    if constexpr (std::is_same_v<T, float>){
      // Query optimal workspace size for LAPACKE_dgeqrf
      float work_query;
      lapack_int lwork_query = -1;
      info = LAPACKE_sgeqrf_work(LAPACK_ROW_MAJOR, n_rows_, n_cols_, data_.data(), n_cols_, 
				 tau_.data(), &work_query, lwork_query);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgeqrf workspace query failed with info = " + std::to_string(info));
      }

      // Resize work array
      lwork_ = static_cast<lapack_int>(work_query);
      work_.resize(lwork_);

      // Perform QR factorization
      info = LAPACKE_sgeqrf_work(LAPACK_ROW_MAJOR, n_rows_, n_cols_, data_.data(), n_cols_, 
				 tau_.data(), work_.data(), lwork_);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgeqrf failed with info = " + std::to_string(info));
      }
    }
    else if constexpr (std::is_same_v<T, double>){
      // Query optimal workspace size for LAPACKE_dgeqrf
      double work_query;
      lapack_int lwork_query = -1;
      info = LAPACKE_dgeqrf_work(LAPACK_ROW_MAJOR, n_rows_, n_cols_, data_.data(), n_cols_, 
				 tau_.data(), &work_query, lwork_query);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgeqrf workspace query failed with info = " + std::to_string(info));
      }

      // Resize work array
      lwork_ = static_cast<lapack_int>(work_query);
      work_.resize(lwork_);

      // Perform QR factorization
      info = LAPACKE_dgeqrf_work(LAPACK_ROW_MAJOR, n_rows_, n_cols_, data_.data(), n_cols_, 
				 tau_.data(), work_.data(), lwork_);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dgeqrf failed with info = " + std::to_string(info));
      }
    }
      
    is_QR_factorized_ = true;       
  }

  void lssolve(Vector<T>& x, Vector<T>& b){
    if (b.size() != static_cast<size_t>(n_rows_)) {
      throw std::invalid_argument("Right-hand-side vector size does not match matrix rows");
    }

    // Perform QR factorization if not already done
    if (!is_QR_factorized_) {
      QRFactorize();
    }

    if constexpr(std::is_same_v<T, float>){
      lapack_int info;
      lapack_int nrhs = 1; // One right-hand-side
      Vector<float> b_copy = b; // Avoid modifying input b

      // Apply Q^T to b using LAPACKE_dormqr
      lapack_int k = std::min(n_rows_, n_cols_);
      info = LAPACKE_sormqr_work(LAPACK_ROW_MAJOR, 'L', 'T', n_rows_, nrhs, k, 
				 data_.data(), n_cols_, tau_.data(), 
				 b_copy.data(), n_rows_, work_.data(), lwork_);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dormqr failed with info = " + std::to_string(info));
      }

      // Solve triangular system Rx = c_1 (first n elements of Q^T b)
      // x.resize(n_cols_);
      info = LAPACKE_strsv(LAPACK_ROW_MAJOR, 'U', 'N', 'N', n_cols_, 
			   data_.data(), n_cols_, b_copy.data(), 1);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dtrsv failed with info = " + std::to_string(info));
      }

      // Copy solution to x
      for (lapack_int i = 0; i < n_cols_; ++i) {
	x(i) = b_copy(i);
      }
    }
    else if constexpr(std::is_same_v<T, double>){
      lapack_int info;
      lapack_int nrhs = 1; // One right-hand-side
      Vector<double> b_copy = b; // Avoid modifying input b

      // Apply Q^T to b using LAPACKE_dormqr
      lapack_int k = std::min(n_rows_, n_cols_);
      info = LAPACKE_dormqr_work(LAPACK_ROW_MAJOR, 'L', 'T', n_rows_, nrhs, k, 
				 data_.data(), n_cols_, tau_.data(), 
				 b_copy.data(), n_rows_, work_.data(), lwork_);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dormqr failed with info = " + std::to_string(info));
      }

      // Solve triangular system Rx = c_1 (first n elements of Q^T b)
      // x.resize(n_cols_);
      info = LAPACKE_dtrsv(LAPACK_ROW_MAJOR, 'U', 'N', 'N', n_cols_, 
			   data_.data(), n_cols_, b_copy.data(), 1);
      if (info != 0) {
	throw std::runtime_error("LAPACKE_dtrsv failed with info = " + std::to_string(info));
      }

      // Copy solution to x
      for (lapack_int i = 0; i < n_cols_; ++i) {
	x(i) = b_copy(i);
      }
    }
  }
  
  void vmult(Vector<T>& res, const Vector<T>& inp, const bool adding=false) const{
    /*
      Matrix-vector-multiplication.
   
      The optional parameter
      <tt>adding</tt> determines, whether the
      result is stored in <tt>w</tt> or added
      to <tt>w</tt>.
   
      if (adding)
      <i>w += Av</i>
   
      if (!adding)
      <i>w = Av</i>
   
      Source and destination must
      not be the same vector.
    */

    if (n_cols() != inp.size()) throw std::invalid_argument("Matrix and vector dimensions mismatch");
    if (res.size() != n_rows()) throw std::invalid_argument("Result vector size mismatch");

    T alpha = 1.0;
    T beta = adding ? 1.0 : 0.0;
    

    if constexpr (std::is_same_v<T, float>) {
      cblas_sgemv(CblasRowMajor, CblasNoTrans,
		  n_rows_, n_cols_,
		  alpha, data_.data(), n_cols_,
		  inp.begin(), 1,
		  beta, res.begin(), 1);
      
    }
    else if constexpr (std::is_same_v<T, double>) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans,
		  n_rows_, n_cols_,
		  alpha, data_.data(), n_cols_,
		  inp.begin(), 1,
		  beta, res.begin(), 1);	
    }
  }

  /// w += (*this) * v
  void vmult_add (Vector<T> &w,
                  const Vector<T> &v) const
  {
    vmult(w, v, true);
  }
  
  
  void diagadd (const T src)
  {
    //Assert (!this->empty(), ExcEmptyMatrix());
    assert (n_rows_ == n_cols_);

    for (size_t i=0; i<n_rows_; ++i)
      (*this)(i,i) += src;
  }
  
  
public://///////////////////////////////////////// Memory Management
public://///////////////////////////////////////// Iterators and Range-based Access
public://///////////////////////////////////////// Utility of Query functions
  size_t n_rows() const { return n_rows_; }
  size_t n_cols() const { return n_cols_; }

  T* data(){return data_.data();}
  
public://///////////////////////////////////////// Interoperability and IO

private://///////////////////////////////////////// for priviate data members
  size_t n_rows_;
  size_t n_cols_;
  std::vector<T> data_;


  bool is_QR_factorized_ = false;
  std::vector<T> work_;/// workspace for LAPACK routines
  int lwork_;///Optimal workspace size
  std::vector<T> tau_;/// for QR factorization using DGEQRF  
  
};

#endif//_DENSEMATRIX_H_
