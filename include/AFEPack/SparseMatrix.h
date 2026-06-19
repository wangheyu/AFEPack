#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <AFEPack/Vector.h>
#include <AFEPack/DenseMatrix.h>
#include <AFEPack/SparsityPattern.h>
#include <iostream>
#include <set>

#include <AFEPack/Miscellaneous.h>

AFEPACK_OPEN_NAMESPACE

template <typename T> class SparseMatrix;

template <typename T>
class SparseMatrix : public Vecma<SparseMatrix<T>, T>
{

 public:
  SparseMatrix(){};
  SparseMatrix(SparsityPattern&);
  SparseMatrix(const SparsityPattern&);
  ~SparseMatrix(){};

  SparseMatrix<T> & copy_from (const SparseMatrix<T> &src);

  SparseMatrix<T> &operator = (const T s);
  SparseMatrix<T> &operator = (const SparseMatrix<T>& src);
  
 public:
  void reinit(SparsityPattern&);
  void reinit(const SparsityPattern&);
  void clear_impl();
      
 public:
  const SparsityPattern& get_sparsity_pattern() const;//original deallii version.
  //SparsityPattern* get_sparsity_pattern();//original vecma version, not;
  const SparsityPattern* get_sparsity_pattern_ptr() const;//original vecma version, not;

  void set_entry_impl(const unsigned int&, const unsigned int&, const T&);
  void add_entry_impl(const unsigned int&, const unsigned int&, const T&);      

 public:
  void print_impl();
  void print_sparseMatrix();

 public:
  const unsigned int& m() const;
  const unsigned int& get_m_impl() const;
  const unsigned int& n() const;
  const unsigned int& get_n_impl() const;      
  T* get_data_impl();
  const T* get_data_impl() const;
  T& global_entry(const unsigned int&);
  const T& global_entry(const unsigned int&) const;  

public:
  T& diag_element_impl(const unsigned int&);
  const T& diag_element_impl(const unsigned int&) const;

  T el(const unsigned int&, const unsigned int&) const;
  
 public:
  /// old: des = this * src
  void mult_v_impl(Vector<T>&, Vector<T>&) const;
  
  void mult_v_impl(Vector<T>&, const Vector<T>&) const;
  /// des = this^T * src
  void Tmult_v_impl(Vector<T>&, Vector<T>&);
  /// alpha * this * x + beta * y      
  void mult_v_impl(const T&, Vector<T>&, const T&, Vector<T>&) const;
  /// alpha * this^T * x + beta * y      
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

  /// des = alpha * this * srcB + beta * des      
  void mult_m_impl(const T&, DenseMatrix<T>&, const T&, DenseMatrix<T>&);
  void mult_m_impl(const T&, DenseMatrix<T>&, const T&, SparseMatrix<T>&);
  void mult_m_impl(const T&, SparseMatrix<T>&, const T&, DenseMatrix<T>&);
  void mult_m_impl(const T&, SparseMatrix<T>&, const T&, SparseMatrix<T>&);
  /// des = alpha * this^T * srcB + beta * des      
  void Tmult_m_impl(const T&, DenseMatrix<T>&, const T&, DenseMatrix<T>&);
  void Tmult_m_impl(const T&, DenseMatrix<T>&, const T&, SparseMatrix<T>&);
  void Tmult_m_impl(const T&, SparseMatrix<T>&, const T&, DenseMatrix<T>&);
  void Tmult_m_impl(const T&, SparseMatrix<T>&, const T&, SparseMatrix<T>&);

 private:
  /// we use row compressed storage to store a sparse matrix.  
  /// std::vector<T> entries; /// those "nonzero" entries

  //origin version:
  const SparsityPattern * sparsity_pattern; /// need to attach this to an external sparsity pattern. Q: how to prevent deleting the object?
  //const SparsityPattern * sparsity_pattern; /// need to attach this to an external sparsity pattern. Q: how to prevent deleting the object?

  bool is_square = false;

  int initialized = 0;
    
 public:
  //const unsigned int& n_nonzero_elements() const;
  const unsigned int n_nonzero_elements() const;

  const int is_initialize() const {return initialized;};
};

AFEPACK_CLOSE_NAMESPACE
#endif//_SPARSEMATRIX_H_
