#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <cstddef>
#include <cassert>

#include <AFEPack/MatrixBase.h>
#include <AFEPack/Vector.h>
#include <AFEPack/SparsityPattern.h>



// SparseMatrix class in CRS format
template <FloatingPoint T, Integral IndexType = size_t>
class SparseMatrix : public MatrixBase<SparseMatrix<T, IndexType>, T>
{

public://///////////////////////////////////////// constructors and Initilization
  SparseMatrix(){};
  SparseMatrix(SparsityPattern& sp){
    reinit(sp);
  };
  void reinit(SparsityPattern& sp){
    pattern_ = &sp;
    n_rows_ = pattern_->n_rows();
    n_cols_ = pattern_->n_cols();

    values_.resize(pattern_->n_nnz());    
  }
  ~SparseMatrix(){};
 

  
public://///////////////////////////////////////// Element Access and Modification
  void add (const size_t i,
	    const size_t j,
	    const T value){
    const size_t* row_ptr = pattern_->get_rowstart_indices();
    const size_t* col_indices = pattern_->get_column_numbers();

    for(size_t k = row_ptr[i];k < row_ptr[i + 1];++ k){
      if(col_indices[k] == j){
	values_[k] += value;
	return;
      }
    }

    /// Need to test if there is entry at the position (i, j).
    
  }
  
  T diag_element (const size_t i) const{/// suppose the firsrt entry in each row is the diagnoal element
    const size_t* row_ptr = pattern_->get_rowstart_indices();
    return values_[row_ptr[i]];
  };


  T el (const size_t i, const size_t j) const{/// also return entry at (i, j).
    const size_t* row_ptr = pattern_->get_rowstart_indices();
    const size_t* col_indices = pattern_->get_column_numbers();

    for(size_t k = row_ptr[i];k < row_ptr[i + 1];++ k){
      if(col_indices[k] == j){
	return values_[k];
      }
    }

    return 0;/// if there is such a position (i,j) in the pattern, it means that entry is 0.
    
  }
public://///////////////////////////////////////// Memory Management

  /// Back to the default construction
  void clear(){
    n_rows_ = 0;
    n_cols_ = 0;
    pattern_ = nullptr;
    values_.clear();
  };  
public://///////////////////////////////////////// Iterators and Range-based Access
public://///////////////////////////////////////// Utility of Query Access
  size_t m() const {return n_rows_;}
  size_t n() const {return n_cols_;}

  T  global_entry (const size_t i) const{
    return values_[i];
  }
  T& global_entry (const size_t i){
    return values_[i];
  }
  size_t n_nonzero_elements() const{
    return pattern_->n_nnz();
  }
  const SparsityPattern& get_sparsity_pattern () const{
    return *pattern_;
  }
public://///////////////////////////////////////// Interoperability and IO

  
  void vmult (Vector<T> &dst, const Vector<T> &src) const{
    assert(dst.size() == n_rows_);
    assert(n_cols_ == src.size());

    const size_t* row_ptr = pattern_->get_rowstart_indices();
    const size_t* col_indices = pattern_->get_column_numbers();

    dst = 0;

    for(size_t i = 0;i < n_rows_;++ i){
      for(size_t j = row_ptr[i];j < row_ptr[i + 1];++ j){
	dst(i) += values_[j] * src(col_indices[j]);
      }
    }    
  }

  void Tvmult (Vector<T> &dst, const Vector<T> &src) const{
    assert(n_rows_ == src.size());
    assert(n_cols_ == dst.size());

    const size_t* row_ptr = pattern_->get_rowstart_indices();
    const size_t* col_indices = pattern_->get_column_numbers();

    dst = 0;
    

    for(size_t i = 0;i < n_rows_;++ i){
      for(size_t j = row_ptr[i];j < row_ptr[i + 1];++ j){
	dst(col_indices[j]) += values_[j] * src(i);
      }
    }

  }

private://///////////////////////////////////////// for private data member

  size_t n_rows_;
  size_t n_cols_;

  SparsityPattern * pattern_;
  std::vector<T> values_;


 
};
#endif//_SPARSEMATRIX_H_
