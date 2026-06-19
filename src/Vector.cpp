#include <AFEPack/Vector.h>

AFEPACK_OPEN_NAMESPACE

template <typename T>
Vector<T>::Vector(const unsigned int& _n, const T& _entry)
{
  std::vector<T>::resize(_n, _entry);
  n = _n;
}

template <typename T>
void Vector<T>::reinit_impl(const unsigned int& _n)
{
  std::vector<T>::clear();
  std::vector<T>::resize(_n);
  n = _n;
}

template <typename T>
void Vector<T>::resize_impl(const unsigned int& _m, const T& _entry)
{
  std::vector<T>::clear();
  std::vector<T>::resize(_m, _entry);
  n = _m;
}

template <typename T>
void Vector<T>::clear_impl()
{
  std::vector<T>::clear();
}

template <typename T>
const unsigned int& Vector<T>::get_n_impl() const
{
  return n;
}

template <typename T>
const unsigned int& Vector<T>::size() const
{
  return n;
}

template <typename T>
T& Vector<T>::operator()(const unsigned int& _n)
{
  return (*this)[_n];
}

template <typename T>
const T& Vector<T>::operator()(const unsigned int& _n) const
{
  return (*this)[_n];
}

template <typename T>
T* Vector<T>::get_data_impl()
{
  return std::vector<T>::data();
}

template <typename T>
const T* Vector<T>::get_data_impl() const
{
  return std::vector<T>::data();
}

template <typename T>
void Vector<T>::rescale_impl(const T& r)
{
  cblas_dscal(n, r, get_data_impl(), 1);
}


template <typename T>
T Vector<T>::innerProduct(Vector<T>& src)
{
  return cblas_ddot(n, get_data_impl(), 1, src.get_data_impl(), 1);
}

template <typename T>
void Vector<T>::add_v_impl(const T& alpha,const Vector<T>& des)
{
  //cblas_daxpy(n, alpha, get_data_impl(), 1, des.get_data_impl(), 1);
  cblas_daxpy(n, alpha, des.get_data_impl(), 1, get_data_impl(), 1);
}

template <typename T>
void Vector<T>::sadd_v_impl(const T& s,const Vector<T>& src)
{
  cblas_dscal(n, s, get_data_impl(), 1);
  cblas_daxpy(n, 1, src.get_data_impl(), 1, get_data_impl(), 1);
}

template <typename T>
void Vector<T>::add(const Vector<T>& des)
{
  cblas_daxpy(n, 1.0, des.get_data_impl(), 1, get_data_impl(), 1);
}

template <typename T>
Vector<T> &Vector<T>::operator += (const Vector<T> &v)
{
  assert(n == v.size());
  add (v);
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator -= (const Vector<T> &v)
{
  //Assert (vec_size!=0, ExcEmptyObject());
  //Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  cblas_daxpy(n, -1.0, v.get_data_impl(), 1, get_data_impl(), 1);
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator = (const T s)
{
  std::vector<T>::clear();
  this->resize_impl(n, s);
  //std::vector<T>::resize(n, s);
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator *= (const T s)
{
  this->rescale_impl(s);
  return *this;
}

template <typename T>
T Vector<T>::operator * (const Vector<T> &v) const
{
  assert(n==v.size());
  double sum = 0.;
  for(size_t i = 0;i < n;i ++)
    {
      sum += (*this)(i) * v(i);
    }
  return sum;
}

template <typename T>
T Vector<T>::l1_norm_impl()
{
  return cblas_dasum(n, get_data_impl(), 1);/// for double by default!!!!!!!!!!!!!!!
}


template <typename T>
T Vector<T>::l2_norm_impl()
{
  return cblas_dnrm2(n, get_data_impl(), 1);/// for double by default!!!!!!!!!!!!!!
}

template <typename T>
T Vector<T>::l8_norm_impl()
{
  return fabs(get_data_impl()[cblas_idamax(n, get_data_impl(), 1)]);/// for double by default!!!!!!!!!!!!!!!!!!

}


template <typename T>
void Vector<T>::print_impl()
{
  T* data = get_data_impl();

  std::cout << "[ ";
  for(unsigned int i = 0;i < n;++ i){
    std::cout << data[i] << " ";
  }
  std::cout << " ]'" << std::endl;
}

template <typename T>
void Vector<T>::equ(const T a, const Vector<T> &u)
{
  assert(n == u.size());
  for (size_t i=0; i < n; ++i)
    { //val[i] = a * u.val[i] + b * v.val[i];
      get_data_impl()[i] = a * u.get_data_impl()[i];
    }
}

template <typename T>
void Vector<T>::equ(const T a, const Vector<T> &u,
		    const T b, const Vector<T> &v)
{
  assert(n == u.size());
  assert(n == v.size());
  for (size_t i = 0;i < n; ++i)
    { //val[i] = a * u.val[i] + b * v.val[i];
      get_data_impl()[i] = a * u.get_data_impl()[i] + b * v.get_data_impl()[i];
    }
}

template class Vector<double>;

AFEPACK_CLOSE_NAMESPACE
