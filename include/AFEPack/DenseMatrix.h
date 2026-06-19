#ifndef _DENSEMATRIX_H_
#define _DENSEMATRIX_H_

#include <AFEPack/Vector.h>
#include <AFEPack/SparseMatrix.h>
#include <iostream>

#include <AFEPack/Miscellaneous.h>

AFEPACK_OPEN_NAMESPACE

template <typename T> class DenseMatrix;

template <typename T>
class DenseMatrix : public Vecma<DenseMatrix<T>, T>
{
public:
  DenseMatrix(){};
  ~DenseMatrix(){};

  DenseMatrix(const unsigned int& _m);
  
  DenseMatrix(const unsigned int& _m, const unsigned int& _n, const T& _entry = 0.);

public:
  //void reinit_impl(const unsigned int&, const unsigned int&, const T&);
  void reinit_impl(const unsigned int&, const unsigned int&);
  //void reinit_impl(const unsigned int&, const T&);
public:
  void set_entry_impl(const unsigned int&, const unsigned int&, const T&);
  void add_entry_impl(const unsigned int&, const unsigned int&, const T&);

public:
  const unsigned int& get_m_impl() const;
  const unsigned int& n_rows() const{return get_m_impl();};
  const unsigned int& get_n_impl() const;      
  const unsigned int& n_cols() const{return get_n_impl();};
  T* get_data_impl();
  const T* get_data_impl() const;
  
public:
  /// des = this + src
  void add_m_impl(DenseMatrix<T>&);
  void add_m_impl(SparseMatrix<T>&);

  /// des = this * src
  void mult_v_impl(Vector<T>&, Vector<T>&) const;
  /// des = this * src
  void mult_v_impl(Vector<T>&, Vector<T>&);
  /// des = this^T * src
  void Tmult_v_impl(Vector<T>&, Vector<T>&);
  /// alpha * A * x + beta * y
  void mult_v_impl(const T&, Vector<T>&, const T&, Vector<T>&) const;
  /// alpha * A * x + beta * y
  void mult_v_impl(const T&, Vector<T>&, const T&, Vector<T>&);
  /// alpha * A^T * x + beta * y
  void Tmult_v_impl(const T&, Vector<T>&, const T&, Vector<T>&);  
    
  /// des = this * src
  void mult_m_impl(DenseMatrix<T>&, DenseMatrix<T>&);
  void mult_m_impl(DenseMatrix<T>&, SparseMatrix<T>&);  
  void mult_m_impl(SparseMatrix<T>&, DenseMatrix<T>&);
  void mult_m_impl(SparseMatrix<T>&, SparseMatrix<T>&);  
  /// des = this^T * src
  void Tmult_m_impl(DenseMatrix<T>&, DenseMatrix<T>&);
  void Tmult_m_impl(DenseMatrix<T>&, SparseMatrix<T>&);  
  void Tmult_m_impl(SparseMatrix<T>&, DenseMatrix<T>&);
  void Tmult_m_impl(SparseMatrix<T>&, SparseMatrix<T>&);  
  
  /// alpha * this * A + beta * B  
  void mult_m_impl(const T&, DenseMatrix<T>&, const T&, DenseMatrix<T>&);
  void mult_m_impl(const T&, DenseMatrix<T>&, const T&, SparseMatrix<T>&);
  void mult_m_impl(const T&, SparseMatrix<T>&, const T&, DenseMatrix<T>&);
  void mult_m_impl(const T&, SparseMatrix<T>&, const T&, SparseMatrix<T>&);  
  /// alpha * this^T * A + beta * B  
  void Tmult_m_impl(const T&, DenseMatrix<T>&, const T&, DenseMatrix<T>&);
  void Tmult_m_impl(const T&, DenseMatrix<T>&, const T&, SparseMatrix<T>&);
  void Tmult_m_impl(const T&, SparseMatrix<T>&, const T&, DenseMatrix<T>&);
  void Tmult_m_impl(const T&, SparseMatrix<T>&, const T&, SparseMatrix<T>&);

  void gauss_jordan();
  //void compute_generalized_eigenvalues_symmetric(DenseMatrix<number> &B, const number lower_bound, const number upper_bound, const number abs_accuracy, Vector<number> &eigenvalues, std::vector<Vector<number> > &eigenvectors, const int itype = 1);
  void compute_generalized_eigenvalues_symmetric(DenseMatrix<T> &, std::vector<Vector<T> > &, int itype = 1);


  std::complex<T> eigenvalue (const size_t i) const;

public:
  void output();
  void output_denseMatrix();

public:
  void print_impl();

private:
  unsigned int m = 0;
  unsigned int n = 0;

    /**
     *    * Real parts of eigenvalues or
     *       * the singular values. Filled by
     *          * compute_eigenvalues() or compute_svd().
     *             */
  std::vector<T> wr;

      /**
       *    * Imaginary parts of
       *       * eigenvalues. Filled by
       *          * compute_eigenvalues.
       *             */
  std::vector<T> wi;

  mutable std::vector<T> work;

  bool is_square = false;

  const char N = 'N';
  const char U = 'U';
  const char V = 'V';

public:
  T operator() (const unsigned int& row_idx, const unsigned int& col_idx) const
  {
    return std::vector<T>::operator[](row_idx * n + col_idx);
  };
  T &operator() (const unsigned int& row_idx, const unsigned int& col_idx)
  {
    return std::vector<T>::operator[](row_idx * n + col_idx);
  };
  
  DenseMatrix<T> &operator = (double val);
};

AFEPACK_CLOSE_NAMESPACE
#endif//_DENSEMATRIX_H_
