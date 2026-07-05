#ifndef _SPARSITYPATTERN_H
#define _SPARSITYPATTERN_H
#pragma once
#include <limits>
#include <algorithm>
#include <utility>
#include <cassert>

const size_t invalid_size_t = std::numeric_limits<size_t>::max();

// SparsityPattern class
class SparsityPattern {
public://///////////////////////////////////////// constructors and Initilization
  SparsityPattern(){};

  void reinit (const size_t m,
	       const size_t n,
	       const std::vector<unsigned int>& row_lengths)
  {
    // Validate input dimensions
    if (m == 0 || n == 0) {
      n_rows_ = 0;
      n_cols_ = 0;
      nnz_ = 0;
      row_ptr_.clear();
      col_indices_.clear();
      is_compressed_ = false;
      max_row_length_ = 0;
      capacity_ = 0;
      is_valid_ = false;
      is_sorted_ = false;
      return;
    }

    // Set basic dimensions
    n_rows_ = m;
    n_cols_ = n;

    is_diag_element_optimized_ = (m == n);/// if it is a square matrix, put the diagonal element the first one in each row    

    // Calculate total number of non-zeros and max row length
    nnz_ = 0;
    max_row_length_ = 0;
    for (size_t i = 0; i < n_rows_; ++i) {
      nnz_ += row_lengths[i];
      max_row_length_ = std::max(max_row_length_, static_cast<size_t>(row_lengths[i]));
    }

    // Reserve and resize space for column indices
    capacity_ = nnz_ > 0 ? nnz_ : n_rows_ * n_cols_ / 4; // Heuristic for initial capacity if nnz_ is 0
    col_indices_.reserve(capacity_);
    col_indices_.resize(nnz_, invalid_size_t);

    // Initialize row pointers
    row_ptr_.resize(n_rows_ + 1);
    row_ptr_[0] = 0;
    for (size_t i = 0; i < n_rows_; ++i) {
      row_ptr_[i + 1] = row_ptr_[i] + row_lengths[i];
    }

    if(is_diag_element_optimized_){
      for (size_t i = 0; i < n_rows_; ++i) {
	col_indices_[row_ptr_[i]] = i;/// set the first entry be the diagonal one
      }      
    }

    

    // Initialize other members
    is_compressed_ = false;
    is_valid_ = true;
    is_sorted_ = true;

  }
    
  void reinit (const size_t m,
	       const size_t n,
	       const size_t max_row_lengths)
  {
    // Validate input dimensions
    if (m == 0 || n == 0) {
      n_rows_ = 0;
      n_cols_ = 0;
      nnz_ = 0;
      row_ptr_.clear();
      col_indices_.clear();
      is_compressed_ = false;
      max_row_length_ = 0;
      capacity_ = 0;
      is_valid_ = false;
      is_sorted_ = false;
      is_diag_element_optimized_ = false;
      return;
    }

    // Set basic dimensions
    n_rows_ = m;
    n_cols_ = n;
    nnz_ = n_rows_ * max_row_lengths;
    max_row_length_ = max_row_lengths;

    is_diag_element_optimized_ = (m == n);/// if it is a square matrix, put the diagonal element the first one in each row

    // Reserve and resize space for column indices
    capacity_ = nnz_;
    col_indices_.reserve(capacity_);
    col_indices_.resize(nnz_, invalid_size_t);///    

    // Initialize row pointers
    row_ptr_.resize(n_rows_ + 1);
    row_ptr_[0] = 0;
    for (size_t i = 0; i < n_rows_; ++i) {
      row_ptr_[i + 1] = row_ptr_[i] + max_row_length_;
    }

    if(is_diag_element_optimized_){
      for (size_t i = 0; i < n_rows_; ++i) {
	col_indices_[row_ptr_[i]] = i;/// set the first entry be the diagonal one      
      }
    }

    // Initialize other members
    is_compressed_ = false;
    is_valid_ = true;
    is_sorted_ = true;
  }


  
  SparsityPattern (const size_t m,
                   const size_t n,
                   const std::vector<unsigned int> &row_lengths)
  {
    reinit(m, n, row_lengths);
  }
  
  SparsityPattern (const size_t m,
                   const size_t n,
                   const size_t max_row_lengths)
  {
    reinit(m, n, max_row_lengths);
  }

  SparsityPattern (const size_t m,
                   const std::vector<unsigned int> &row_lengths)
  {
    reinit(m, m, row_lengths);
  }
  SparsityPattern (const size_t m,
                   const size_t max_row_lengths)
  {
    reinit(m, m, max_row_lengths);
  }
  
public://///////////////////////////////////////// Element Access and Modification
  // void add_entry(size_t row, size_t col) {
  //   if (row >= n_rows_ || col >= n_cols_) throw std::out_of_range("Invalid indices");
  //   row_ptr_[row + 1]++;
  //   col_indices_.push_back(col);
  // }

  // void finalize() {
  //   for (size_t i = 1; i <= n_rows_; ++i) {
  //     row_ptr_[i] += row_ptr_[i - 1];
  //   }
  // }

  /*
    Step 1. Check if the matrix has been compressed already. If yes,
    just return, otherwise

    Step 2. Introduce two auxiliary vectors, tmp_col_indices_(nnz) for
    storing the final result (will use move to assign the results to
    col_indices finally), and tmp_col_per_row_(max_row_length_) to
    store the temporary entries from each row

    Step 3. Loop all rows

      Step 3.1. Copy entries from current row to tmp_col_per_row_, and
      optimize and sort it

      Step 3.2. Copy the optimized and sorted tmp_col_per_row_ to
      tmp_col_indices_, and revise the corresponding row_ptr_

    Step 4. Move tmp_col_indices_ to col_indices (now col_indices will
    be the final results, and tmp_col_indices will be empty)

    Step 5. Check if dimensions of row_ptr_ and col_indices_ are
    consistent, and update all associated variables such as
    is_compressed_.
    
   */
  void compress()
  {
    if(is_compressed_) return; /// do nothing when is_compressed_ is true

    /// calculate the total number for non-zeros.
    nnz_ -= std::count(col_indices_.begin(), col_indices_.end(), invalid_size_t);/// the entry with invalid_size_t is zero-entry
    std::vector<size_t> tmp_col_indices_(nnz_);
    std::vector<size_t> tmp_col_per_row_(max_row_length_);
    auto it = tmp_col_indices_.begin();
    size_t n_current_row_ = 0;
    size_t new_row_ptr_ = 0;

    for(int i = 0;i < n_rows_;++ i){
      for(size_t j = row_ptr_[i];j < row_ptr_[i + 1];++ j, n_current_row_++){
	if(col_indices_[j] == invalid_size_t) break;

	tmp_col_per_row_[n_current_row_] = col_indices_[j];
      }/// when jump out, tmp_col_pre_row_ stores all non-zero entries
       /// from current row, and n_current_row_ is the total number
       /// for those non-zero entries.

      std::sort(is_diag_element_optimized_ ? tmp_col_per_row_.begin() + 1 : tmp_col_per_row_.begin(),
		tmp_col_per_row_.begin() + n_current_row_);/// non-descring order

      it = std::copy(tmp_col_per_row_.begin(), tmp_col_per_row_.begin() + n_current_row_, it);

      //row_ptr_[i + 1] = row_ptr_[i] + n_current_row_;/// update row_ptr_


      row_ptr_[i] = new_row_ptr_;
      new_row_ptr_ += n_current_row_;

      max_row_length_ = max_row_length_ > n_current_row_ ? max_row_length_ : n_current_row_;/// to update max_row_length_

      n_current_row_ = 0;/// for next row processcing      
    }/// at this moment, both row_ptr_ and tmp_col_indices_ are fullfilled.
    row_ptr_[n_rows_] = new_row_ptr_;/// to modify the last entry of row_ptr_.
    
    col_indices_ = std::move(tmp_col_indices_);/// tmp_col_indices_ should be empty now.

    is_compressed_ = true;

    /// done
    
  }
  
  size_t max_entries_per_row () const
  {
    if (!is_compressed_)
      return max_row_length_;

    size_t m = 0;
    for (size_t i=1; i<n_rows_; ++i)
      m = std::max (m, row_ptr_[i]-row_ptr_[i-1]);

    return m;
  }

  void add (const size_t i, const size_t j){
    
    for (size_t k = row_ptr_[i]; k < row_ptr_[i+1]; k++){
      // entry already exists
      if (col_indices_[k] == j) return;
      // empty entry found, put new entry here
      if (col_indices_[k] == invalid_size_t){
	col_indices_[k] = j;
	return;
      };
    };
  }

  size_t operator() (const size_t i,
		     const size_t j) const
  {
    //assert ((row_ptr_!=nullptr) && (col_indices_!=nullptr));
    assert (i<n_rows_);
    assert (j<n_cols_);
    assert (is_compressed_);

    // let's see whether there is
    // something in this line
    if (row_ptr_[i] == row_ptr_[i+1])
      return invalid_size_t;

    // If special storage of diagonals
    // was requested, we can get the
    // diagonal element faster by this
    // query    if (is_diag_element_optimized_ && (i==j))
      return row_ptr_[i];

    // all other entries are sorted, so
    // we can use a binary search
    // algorithm
    //
    // note that the entries are only
    // sorted upon compression, so this
    // would fail for non-compressed
    // sparsity patterns; however, that
    // is why the Assertion is at the
    // top of this function, so it may
    // not be called for noncompressed
    // structures.
    const size_t *sorted_region_start = (is_diag_element_optimized_ ?
					    &col_indices_[row_ptr_[i]+1] :
					    &col_indices_[row_ptr_[i]]);
    // const size_type *const p
    //   = Utilities::lower_bound<const size_type *> (sorted_region_start,
    // 						   &colnums[rowstart[i+1]],
    // 						   j);

    const size_t *const p
      = std::lower_bound(sorted_region_start,
			 &col_indices_[row_ptr_[i+1]],
			 j);    
    if ((p != &col_indices_[row_ptr_[i+1]])  &&  (*p == j))
      return (p - &col_indices_[0]);
    else
      return invalid_size_t;
    
  }
  

public://///////////////////////////////////////// Memory Management
  void reserve(size_t nnz) { col_indices_.reserve(nnz); }  
public://///////////////////////////////////////// Iterators and Range-based Access
public://///////////////////////////////////////// Utility of Query functions
  size_t n_rows() const { return n_rows_; }
  size_t n_cols() const { return n_cols_; }
  size_t n_nnz() const {return nnz_;}
  size_t n_nonzero_elements() const {return nnz_;}
  const std::vector<size_t>& row_ptr() const { return row_ptr_; }
  const std::vector<size_t>& col_indices() const { return col_indices_; }  
  size_t row_length (const size_t row) const{ return row_ptr_[row + 1] - row_ptr_[row];}
  bool is_compressed() const {return is_compressed_;}
  const size_t* get_rowstart_indices() const {return row_ptr_.data();}
  const size_t* get_column_numbers() const {return col_indices_.data();}
  
public://///////////////////////////////////////// Interoperability and IO

private://///////////////////////////////////////// for private data member

  size_t n_rows_;                     // Number of rows
  size_t n_cols_;                     // Number of columns
  size_t nnz_;                        // Number of non-zero elements
  bool is_compressed_;                   // Compression state
  size_t max_row_length_;             // Max non-zeros in any row
  size_t capacity_;                   // Allocated capacity for vectors
  bool is_valid_;                     // Matrix validity flag
  bool is_sorted_;                    // Are column indices sorted?
  bool is_diag_element_optimized_;       // Is the diagnol element stored first?

  std::vector<size_t> row_ptr_;       // Row pointers for CRS
  std::vector<size_t> col_indices_;   // Column indices of non-zeros
};

#endif//_SPARSITYPATTERN_H
