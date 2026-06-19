#ifndef _SPARSITYPATTERN_H_
#define _SPARSITYPATTERN_H_

#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>

#include <AFEPack/Miscellaneous.h>

AFEPACK_OPEN_NAMESPACE

class SparsityPattern
{
 public:
  SparsityPattern();
  SparsityPattern(const u_int&, const u_int&, const u_int&);
  SparsityPattern(const u_int&, const u_int&, const std::vector<u_int>&);
  SparsityPattern(const u_int&, const std::vector<u_int>&);
  ~SparsityPattern();
  

 public:
  /// initialization
  void reinit(const u_int& _m, const u_int& _n);
  void reinit(const u_int& _m, const u_int& _n, const u_int& _max_entries_per_row);
  void reinit(const u_int& _m, const u_int& _n, const std::vector<u_int>& _entries_per_row);

 public:
  const u_int& get_m() const;
  //u_int& n_rows(){return m;};
  const u_int& n_rows()const{return m;};
  const u_int& get_n() const;
  //u_int& n_cols(){return n;};
  const u_int& n_cols()const{return n;};
  const u_int get_n_nonzeros() const;

  void compress();

  void add(const u_int&, const u_int&);

  std::size_t * get_row();
  const std::size_t * get_row() const;
  std::size_t * get_rowstart_indices();  
  const std::size_t * get_rowstart_indices() const;
  u_int * get_col();
  const u_int * get_col() const;
  u_int * get_column_numbers();
  const u_int * get_column_numbers() const;

  //const u_int get_max_entries_per_row() const{
  //    return max_entries_per_row;  // Add a new function for obtaining the max_entries_per_row. Used in AMGSolver.cpp
  //}
  u_int max_entries_per_row() const; //Add a new function for obtaining the max_entries_per_row. Used in AMGSolver.cpp

  u_int row_length(int i)const {return row[i+1] - row[i];};

  bool is_compressed(); //
  const bool is_compressed() const;
 public:
  
  void set_max_entries_per_row(const u_int&);
  void put_diagEle_first_in_each_row();
  
 public:
  const void print() const;
 private:
  u_int m = 0;
  u_int n = 0;
  u_int n_nonzeros=0; /// including those entries with "0".
  
  bool shape_square = true;
  bool compressed = false;

 private:
  //u_int max_entries_per_row;
  u_int max_row_length;/// initial rough estimate for the maximum length for each row
  
 private:
  /// we use row compressed storage to store a sparse matrix.    
  //std::vector<u_int> row; /// row information of the matrix
  std::vector<size_t> row; /// updated version by Liu C.
  std::vector<u_int> col; /// column informaton of the matrix

 public:
  class ExcNotCompressed : public std::exception{
  public:
    ExcNotCompressed(){}
    ~ExcNotCompressed(){}
    void print_info (std::ostream &out) const {                   
      out << "Not Compressed" << std::endl;                           
    }                        
  };
};

AFEPACK_CLOSE_NAMESPACE

#endif//_SPARSITYPATTERN_H_
