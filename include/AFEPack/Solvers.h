#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <AFEPack/Vector.h>
#include <AFEPack/DenseMatrix.h>
#include <AFEPack/SparseMatrix.h>

#include <cassert>

/// for solving an overdetermined linear system Ax=b, where A is m x n with m > n.

template <typename T>
class LSSolver
{
public:
  LSSolver(){};
  ~LSSolver(){
    delete[] tau_;
  };

public:
  void reinit(DenseMatrix<T>& A)
  {
    assert(A.n_rows() > A.n_cols());

    A_ = &A; ///it also indicates that QR factorization will happen IN-PLACE
    m_ = A_->n_rows();
    n_ = A_->n_cols();
    LDA_ = n_;
    tau_ = new double[LDA_];

    QRFactorization();    
  }

  void solve(Vector<T>& x, const Vector<T>& b)
  {
    x = b;
    
    /// Apply the orthogonal transformation to vector x using LAPACKE_dormqr
    int info = LAPACKE_dormqr(LAPACK_ROW_MAJOR, 'L', 'T', m_, 1, n_, A_->data(), LDA_, tau_, x.data(), 1);

    if (info != 0) {
      std::cerr << "Error in dormqr: " << info << std::endl;
      throw std::runtime_error("Orthogonal transformation failed");
    }

    // Solve the triangular system Rx = Q^Tx using LAPACKE_dtrtrs
    info = LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 'U', 'N', 'N', n_, 1, A_->data(), n_, x.data(), 1);

    if (info != 0) {
        std::cerr << "Error in dtrtrs: " << info << std::endl;
        throw std::runtime_error("Solving triangular system failed");
    }
  }

private:
  double* tau_;

  void QRFactorization()
  {/// after this function, A_ will be changed, and tau_ is filled.

    // Perform the actual QR factorization with the determined workspace
    int info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m_, n_, A_->data(), LDA_, tau_);

    if (info != 0) {
      std::cerr << "Error in dgeqrf: " << info << std::endl;
      throw std::runtime_error("QR factorization failed");
    }    
  }

private:
  DenseMatrix<T> * A_;
  int m_, n_, LDA_;
  
};

#endif//_SOLVER_H_
