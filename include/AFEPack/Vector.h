#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <AFEPack/Vecma.h>
#include <iostream>

#include <AFEPack/Miscellaneous.h>

AFEPACK_OPEN_NAMESPACE

template <typename T>
class Vector : public Vecma<Vector<T>, T>
{
public:
  Vector(){};
  ~Vector(){};
  Vector(const unsigned int&, const T& = 0.);

public:
  //void reinit_impl(const unsigned int&, const T&);
  void reinit_impl(const unsigned int&);
  void reinit_impl(const unsigned int&, bool);
  void resize_impl(const unsigned int&, const T& = 0.);
  void clear_impl();

  /**
   * Assignment <tt>*this = a*u </tt>.
   */
  void equ (const T a, const Vector<T> &u);

  
  /**
   * Assignment <tt>*this = a*u + b*v</tt>.
   */
  void equ (const T a, const Vector<T> &u,
            const T b, const Vector<T> &v);

  /**
   * Add the given vector to the 
   * present one.
   */
  Vector<T> &operator += (const Vector<T> &V);
  /**
   * Subtract the given vector from the
   * present one.
   */
  Vector<T> &operator -= (const Vector<T> &V);

  /**
   * Set all entries of the vector to s;
   *
   */
  Vector<T> &operator = (const T s);

  /**
   * Multiply all entries by s
   */
  Vector<T> &operator *= (const T s);

  /**
   * return the scalar product of the two vectors
   */
  T operator * (const Vector<T> &V) const;

public:
  const unsigned int& get_n_impl() const;
  const unsigned int& size() const;

public:
  T& operator()(const unsigned int&);
  const T& operator()(const unsigned int&) const;  
public:
  T* get_data_impl();
  const T* get_data_impl() const;
  void rescale_impl(const T&);

public:
  void block_write (std::ostream &out) const{};
  void block_read (std::istream &in){};  
public:
  void set_entry(const unsigned int&, const T&){};
  void add_entry(const unsigned int&, const T&){};  
  
  /// uniqely for Vector class
public:
  T innerProduct(Vector<T>&);

  
  /// des = alpha * this + des
  void add_v_impl(const T&, const Vector<T>&);

  /// this = s*this + src;
  void sadd_v_impl(const T&, const Vector<T>&);

  /// A new function, this += alpha* des,
  //void add_av_impl(const T&, const Vector<T>&);

  /// A new function, equals to the operator +=, i.e., this += des;
  void add(const Vector<T>&);

  void add(const T& a, const Vector<T>& src)
  {
    add_v_impl(a, src);
  };
public:
  T l1_norm_impl();
  T l2_norm_impl();
  T l8_norm_impl();

public:
  void print_impl();

private:
  unsigned int n; 
};

AFEPACK_CLOSE_NAMESPACE

#endif//_VECTOR_H_
