#ifndef _VECMA_H_
#define _VECMA_H_

#include <vector>
#include <algorithm>
#include <cblas.h>
#include <cmath>
#include <cassert>

#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>

#include <AFEPack/Miscellaneous.h>

AFEPACK_OPEN_NAMESPACE

template<typename T> class Vector;
template<typename T> class DenseMatrix;
template<typename T> class SparseMatrix;

template <typename cT, typename T>
class Vecma : public std::vector<T>
{
public:
  Vecma(){};
  ~Vecma(){};

public:
  /// for dense and sparse matrix
//  void reinit(const unsigned int& _m, const unsigned int& _n, const T& _entry = 0.){
//    static_cast<cT*>(this)->reinit_impl(_m, _n, _entry);
//  };
   
  /// for dense and possible sparse vector
//  void reinit(const unsigned int& _n, const T& _entry = 0.){
//    static_cast<cT*>(this)->reinit_impl(_n, _entry);
//  };

  // Updated version for avoiding ambiguous error
  void reinit(const unsigned int& _m, const unsigned int& _n){
    static_cast<cT*>(this)->reinit_impl(_m, _n);
  };
  void reinit(const unsigned int& _n){
    static_cast<cT*>(this)->reinit_impl(_n);
  };

  // resize function similar to the vector::resize()
  void resize(const unsigned int& _m, const T& _entry = 0.){
    static_cast<cT*>(this)->resize_impl(_m, _entry);
  }
  void clear(){
    static_cast<cT*>(this)->clear_impl();
  };


  /// accessor
public:
  const unsigned int& get_m() const{
    return static_cast<const cT*>(this)->get_m_impl();
  };  

  const unsigned int& get_n() const{
    return static_cast<const cT*>(this)->get_n_impl();
  };

  /// for dense and sparse matrix
  const T& get_entry(const unsigned int& _m, const unsigned int& _n){
    return static_cast<cT*>(this)->get_entry_impl(_m, _n);
  }

  T* get_data(){
    return static_cast<cT*>(this)->get_data_impl();
  };

  const T* get_data() const {
    return static_cast<const cT*>(this)->get_data_impl();
  };
  
  T& diag_element(const unsigned int& _n){
    return static_cast<cT*>(this)->diag_element_impl(_n);
  };
  //const T& diag_element(const unsigned int& _n) const{
  //  return static_cast<cT*>(this)->diag_element_impl(_n);
  //};
  const T& diag_element(const unsigned int& _n) const {
    return static_cast<const cT*>(this)->diag_element_impl(_n);
  };

  /// modifier
public:
  /// for vectors
  void set_entry(const unsigned int& _idx, const T & _entry){
    static_cast<cT*>(this)->set_impl(_idx, _entry);
  };
  void add_entry(const unsigned int& _idx, const T & _entry){
    static_cast<cT*>(this)->add_impl(_idx, _entry);
  };

  /// for matrix
  void set_entry(const unsigned int& row_idx, const unsigned int& col_idx, const T & entry){
    static_cast<cT*>(this)->set_entry_impl(row_idx, col_idx, entry);
  };
  void add_entry(const unsigned int& row_idx, const unsigned int& col_idx, const T & entry){
    static_cast<cT*>(this)->add_entry_impl(row_idx, col_idx, entry);
  };
  void add(const unsigned int& row_idx, const unsigned int& col_idx, const T & entry){
    static_cast<cT*>(this)->add_entry_impl(row_idx, col_idx, entry);
  };

  /// for vectors, could be for matrix
  //void rescale(const T& r){
  void scale(const T& r){
    static_cast<cT*>(this)->rescale_impl(r);    
  };

public:

  /// for Vector(_des) = DenseMatrix(this) * Vector(_src), or Vector(_des) = SparseMatrix(this) * Vector(_src)
  //void mult_v(Vector<T>& _src, Vector<T>& _des){
  // void mult_v(Vector<T>& _des, Vector<T>& _src){
  void vmult(Vector<T>& _des, Vector<T>& _src)const {
    static_cast<const cT*>(this)->mult_v_impl(_des, _src);
  };

  void vmult(Vector<T>& _des, const Vector<T>& _src)const {
    static_cast<const cT*>(this)->mult_v_impl(_des, _src);
  };
  
  /*
  void vmult(Vector<T>& _des, Vector<T>& _src){
    static_cast<cT*>(this)->mult_v_impl(_des, _src);
  };
*/
  /// for Vector(_des) = DenseMatrix(this)^T * Vector(_src), or Vector(_des) = SparseMatrix(this)^T * Vector(_src)
  //void Tmult_v(Vector<T>& _src, Vector<T>& _des){
  //void Tmult_v(Vector<T>& _des, Vector<T>& _src){
  void Tvmult(Vector<T>& _des, Vector<T>& _src) {
    static_cast<cT*>(this)->Tmult_v_impl(_des, _src);
  };

  /// alpha * A * x + beta * y
  //void mult_v(const T& _alpha, Vector<T>& _src, const T& _beta, Vector<T>& _des)
  void vmult(const T& _alpha, Vector<T>& _src, const T& _beta, Vector<T>& _des) const{
    static_cast<const cT*>(this)->mult_v_impl(_alpha, _src, _beta, _des);
  }
  /*
  void vmult(const T& _alpha, Vector<T>& _src, const T& _beta, Vector<T>& _des){
    static_cast<cT*>(this)->mult_v_impl(_alpha, _src, _beta, _des);
  }*/
  /// alpha * A^T * x + beta * y
  //void Tmult_v(const T& _alpha, Vector<T>& _src, const T& _beta, Vector<T>& _des)
  void Tvmult(const T& _alpha, Vector<T>& _src, const T& _beta, Vector<T>& _des) {
    static_cast<cT*>(this)->Tmult_v_impl(_alpha, _src, _beta, _des);
  }

  

  /// 
  void add_v(const T& _alpha, const Vector<T>& _des){
    static_cast<cT*>(this)->add_v_impl(_alpha, _des);
  };
  /// 
  //void add_v(DenseMatrix<T>& _src, SparseMatrix<T>& _des){
  //    static_cast<cT*>(this)->add_impl(_des);
  //};

  /*
   * Function saa_v() is used for scaling vector and adding a vector
   * *this = s*(*this) + _src;
   */
  void sadd(const T& s, const Vector<T>& _src){
    static_cast<cT*>(this)->sadd_v_impl(s, _src);
  }
  
  /// for DenseMatrix(_des) = DenseMatrix(this) + DenseMatrix(_src), or DenseMatrix(_des) = SparseMatrix(this) + DenseMatrix(_src)
  void add_m(DenseMatrix<T>& _des){
    static_cast<cT*>(this)->add_m_impl(_des);
  };
  /// for SparseMatrix(_des) = DenseMatrix(this) + DenseMatrix(_src), or SparseMatrix(_des) = SparseMatrix(this) + DenseMatrix(_src)
  void add_m(SparseMatrix<T>& _des){
    static_cast<cT*>(this)->add_m_impl(_des);
  };
  
  //////////////////// C = alpha * A * B + beta * C
  /// for DenseMatrix(_des) = DenseMatrix(this) * DenseMatrix(_src), or DenseMatrix(_des) = SparseMatrix(this) * DenseMatrix(_src)
  void mult_m(DenseMatrix<T>& _src, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_src, _des);
  };
  /// for SparseMatrix(_des) = DenseMatrix(this) * DenseMatrix(_src), or SparseMatrix(_des) = SparseMatrix(this) * DenseMatrix(_src)
  void mult_m(DenseMatrix<T>& _src, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_src, _des);
  };
  /// for DenseMatrix(_des) = DenseMatrix(this) * SparseMatrix(_src), or DenseMatrix(_des) = SparseMatrix(this) * SparseMatrix(_src)
  void mult_m(SparseMatrix<T>& _src, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_src, _des);
  };
  /// for SparseMatrix(_des) = DenseMatrix(this) * SparseMatrix(_src), or SparseMatrix(_des) = SparseMatrix(this) * SparseMatrix(_src)
  void mult_m(SparseMatrix<T>& _src, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_src, _des);
  };

  //////////////// C = alpha * A^T * B + beta * C
  /// for DenseMatrix(_des) = DenseMatrix(this) * DenseMatrix(_src), or DenseMatrix(_des) = SparseMatrix(this) * DenseMatrix(_src)
  void Tmult_m(DenseMatrix<T>& _src, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_src, _des);
  };
  /// for SparseMatrix(_des) = DenseMatrix(this) * DenseMatrix(_src), or SparseMatrix(_des) = SparseMatrix(this) * DenseMatrix(_src)
  void Tmult_m(DenseMatrix<T>& _src, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_src, _des);
  };
  /// for DenseMatrix(_des) = DenseMatrix(this) * SparseMatrix(_src), or DenseMatrix(_des) = SparseMatrix(this) * SparseMatrix(_src)
  void Tmult_m(SparseMatrix<T>& _src, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_src, _des);
  };
  /// for SparseMatrix(_des) = DenseMatrix(this) * SparseMatrix(_src), or SparseMatrix(_des) = SparseMatrix(this) * SparseMatrix(_src)
  void Tmult_m(SparseMatrix<T>& _src, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_src, _des);
  };

  /// alpha * this * A + beta * B  
  void mult_m(const T& _alpha, DenseMatrix<T>& _src, const T& _beta, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_alpha, _src, _beta, _des);
  };
  void mult_m(const T& _alpha, DenseMatrix<T>& _src, const T& _beta, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_alpha, _src, _beta, _des);
  };
  void mult_m(const T& _alpha, SparseMatrix<T>& _src, const T& _beta, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_alpha, _src, _beta, _des);
  };
  void mult_m(const T& _alpha, SparseMatrix<T>& _src, const T& _beta, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->mult_m_impl(_alpha, _src, _beta, _des);
  };

  /// alpha * this^T * A + beta * B  
  void Tmult_m(const T& _alpha, DenseMatrix<T>& _src, const T& _beta, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_alpha, _src, _beta, _des);    
  };
  void Tmult_m(const T& _alpha, DenseMatrix<T>& _src, const T& _beta, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_alpha, _src, _beta, _des);    
  };
  void Tmult_m(const T& _alpha, SparseMatrix<T>& _src, const T& _beta, DenseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_alpha, _src, _beta, _des);    
  };
  void Tmult_m(const T& _alpha, SparseMatrix<T>& _src, const T& _beta, SparseMatrix<T>& _des){
    static_cast<cT*>(this)->Tmult_m_impl(_alpha, _src, _beta, _des);    
  };

public:
  T l1_norm(){return static_cast<cT*>(this)->l1_norm_impl();}
  T l2_norm(){return static_cast<cT*>(this)->l2_norm_impl();}
  T l8_norm(){return static_cast<cT*>(this)->l8_norm_impl();}
  T lf_norm(){return static_cast<cT*>(this)->lf_norm_impl();}

public:
  void print(){static_cast<cT*>(this)->print_impl();};

};

AFEPACK_CLOSE_NAMESPACE
#endif//_VECMA_H_
