/**
 * @file   Miscellaneous.h
 * @author Ruo Li <rli@aztec>
 * @date   Fri Dec 14 09:47:04 2007
 * 
 * @brief  
 * 
 * 
 */

#ifndef _Miscellaneous_h_
#define _Miscellaneous_h_

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <wordexp.h>
#include <dlfcn.h>

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <vector>

#include <boost/iostreams/filtering_stream.hpp>
#include <lac/vector.h>

namespace dealii {};
using namespace dealii;

#ifndef __SERIALIZATION__
#define __SERIALIZATION__
#endif

#ifndef MULTITHREAD
#define MULTITHREAD
#endif

#define AFEPACK_OPEN_NAMESPACE namespace AFEPack {
#define AFEPACK_CLOSE_NAMESPACE }

AFEPACK_OPEN_NAMESPACE

typedef void * dlhandle_t;
typedef boost::iostreams::filtering_istream filtering_istream;

void ExpandString(std::string& str);
void StringToWord(const std::string& str, const char& c, std::vector<std::string>& result);
void CombineString(const std::vector<std::string>& prefix, 
                   const std::vector<std::string>& suffix, 
                   std::vector<std::string>& result);

/** 
 * 搜索指定的路径，根据文件名找到包含库文件的完整路径的文件名
 * 
 * @param filename 文件名
 */
std::string FindAFEPackLibraryFilePath(const std::string& filename);

/**
 * 打开共享库并返回一个共享库的句柄
 * 
 */
dlhandle_t AFEPackDLOpen(const std::string& filename);

/** 
 * 搜索指定的路径，找到库文件，并打开为一个流
 * 
 * @param filename 库文件名
 * @param is 打开的流
 */
void OpenAFEPackLibraryFile(const std::string& filename, 
                            filtering_istream& is);

/** 
 * 从共享库中载入一个函数
 * 
 * @param handle 共享库句柄
 * @param sym 函数名
 * @param fun_ptr 返回的函数指针
 */
void LoadLibraryFunction(dlhandle_t& handle,
                         const std::string& sym,
                         dlhandle_t& fun_ptr);

/** 
 * 根据文件名，将其打开为一个过滤的流，流中的 shell 脚本型的注释将被
 * 过滤掉。RealHP 的打开的各种文本型的库文件都事实上使用本函数打开为
 * 一个流，从而其中支持使用 shell 脚本型的注释。
 * 
 * @param filename 文件名
 * @param is 打开的过滤流
 */
void OpenFilteredStream(const std::string& filename,
                        filtering_istream& is);

/**
 * 二维的希尔伯特空间曲线重排。
 */
void hsfc_renumerate(int, double *, double *, int *);
/**
 * 三维的希尔伯特空间曲线重排。
 */
void hsfc_renumerate(int, double *, double *, double *, int *);

template <class V, template <class T> class C>
  V innerProduct(const C<V>& c0, const C<V>& c1) {
  V v = 0;
  typename C<V>::const_iterator 
    the_v0 = c0.begin(), 
    end_v0 = c0.end(),
    the_v1 = c1.begin();
  for (;the_v0 != end_v0; ++ the_v0, ++ the_v1) {
    v += (*the_v0) * (*the_v1);
  }
  return v;
}

/**
 * 此模板函数可以兼容 gcc 4.2.x 版本的情况，这个版本的编译器似乎不允许
 * 对模板参数的缺省情况作模板匹配，从而导致上面的模板 template <class
 * T> class C 不能和 std::vector<value_type, allocator = default
 * allocator> 进行匹配。蔡振宁提供。
 * 
 */
template <class C>
typename C::value_type innerProduct(const C& c0, const C& c1) {
  typename C::value_type v = 0;
  typename C::const_iterator
    the_v0 = c0.begin(),
    end_v0 = c0.end(),
    the_v1 = c1.begin();
  for (;the_v0 != end_v0; ++ the_v0, ++ the_v1) {
    v += (*the_v0) * (*the_v1);
  }
  return v;
}

/**
 * template valued function class. 
 * This class provided the virtual base of a mathematical function with its
 * value and gradient expression provided. Then it can be used to make those 
 * finite element coding more convenient.
 */
template <class value_type>
class Function
{
 public:
  /** Constructor */
  Function() {};
  /** Destructor */
  virtual ~Function() {};
 public:
  /** the value of the function */
  virtual value_type value(const double * p) const {return value_type();};
  /** the gradient of the function */
  virtual std::vector<value_type> gradient(const double * p) const
    {return std::vector<value_type>();};
};

/**
 * function of function. The function adopt pointers to funtion to evaluate 
 * its value and gradient. This make it convenient to use an old C and C++
 * function in the new formation.
 */
template <class value_type>
class FunctionFunction : public Function<value_type>
{
 public:
  typedef value_type (*ValuePrototype)(const double *);
  typedef std::vector<value_type> (*GradientPrototype)(const double *);
 private:
  /**< the pointer to the value function. */
  ValuePrototype			vf;
  /**< the pointer to the gradient function. */
  GradientPrototype			gf;
 public:
  /** Constructor */
 FunctionFunction(value_type (*v)(const double *) = NULL, 
                  std::vector<value_type> (*g)(const double *) = NULL) :
  vf(v), gf(g) {};
 FunctionFunction(const FunctionFunction<value_type>& f) :
  vf(f.vf), gf(f.gf) {};
  /** Destructor */
  virtual ~FunctionFunction() {};
 public:
  /** the pointer of the value function. */
  ValuePrototype valueFunction() const {return vf;};
  /** the pointer of the value function. */
  ValuePrototype& valueFunction() {return vf;};
  /** the pointer of the gradient function. */
  GradientPrototype gradientFunction() const {return gf;};
  /** the pointer of the gradient function. */
  GradientPrototype& gradientFunction() {return gf;};
  operator ValuePrototype() const {return vf;};
  operator ValuePrototype() {return vf;};
  operator GradientPrototype() const {return gf;};
  operator GradientPrototype() {return gf;};
 public:
  /** the value of the function. */
  virtual value_type value(const double * p) const {return (*vf)(p);};
  /** the gradient of the function. */
  virtual std::vector<value_type> gradient(const double * p) const {return (*gf)(p);};
};

/**
 * The vector with length to be "n". Used for vector value basis function
 * finite element space.
 * 
 */
template <int n, class _Tp>
  class nVector : public std::vector<_Tp>
{
 public:
 nVector() : std::vector<_Tp>(n, _Tp()) {};
  explicit nVector(const _Tp& t) : std::vector<_Tp>(n, t) {};
 nVector(const nVector<n,_Tp>& v) : std::vector<_Tp>(n) {
    for (int i = 0;i < n;i ++)
      (*this)[i] = v[i];
  };
  nVector<n,_Tp>& operator=(const nVector<n,_Tp>& v) {
    for (int i = 0;i < n;i ++)
      (*this)[i] = v[i];
    return *this;
  };
};

AFEPACK_CLOSE_NAMESPACE

#ifndef ENABLE_AFEPACK_NAMESPACE
using namespace AFEPack;
#endif

#endif

/**
 * end of file
 * 
 */

