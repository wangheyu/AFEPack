////////////////////////////////////////////////////////////////////////////
// Operator.h :
//

#ifndef _Operator_h_
#define _Operator_h_

#include <typeinfo>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>

#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include "AMGSolver.h"
#include "Geometry.h"
#include "TemplateElement.h"
#include "FEMSpace.h"
#include "HGeometry.h"
#include "BilinearOperator.h"

AFEPACK_OPEN_NAMESPACE

/**
 * This namespace packs a list of operator operating on the finite element
 * functions and analytic functions. Such operations include $L^2$ interpolation,
 * projection and discretization. Generally, an numerical intergration will
 * be involved in the implentation of those operations. A parameter \p{Method}
 * is sumetimes in the parameter list of sume operators to require certain
 * operation style.
 */
namespace Operator {
  typedef int				Method;
  static const Method		MASS_ACCUMULATION = 1;
  static const Method		LEAST_SQUARE = 2;
  static const Method		LOCAL_LEAST_SQUARE = 3;

  /**
   * Interpolate a finite element function \p{src} in one finite element space
   * to a finite element funtion \p{des} in another finite element space.
   */
  template <class value_type, int DIM>
    void L2Interpolate(const FEMFunction<value_type, DIM>& src, FEMFunction<value_type, DIM>& des);
  /**
   * Interpolate an analytic function \p{src} to a finite element funtion \p{des}
   * in another finite element space.
   */
  template <class value_type, int DIM>
    void L2Interpolate(value_type (*)(const double *), FEMFunction<value_type, DIM>& fun);
  /**
   * Interpolate an analytic function \p{src} to a finite element funtion \p{des}
   * in another finite element space.
   */
  template <class value_type, int DIM>
    void L2Interpolate(value_type (*)(const Point<DIM>&), FEMFunction<value_type, DIM>& fun);
  /**
   * Interpolate an analytic function \p{src} to a finite element funtion \p{des}
   * in another finite element space.
   */
  template <class value_type, int DIM>
    void L2Interpolate(const Function<value_type>&, FEMFunction<value_type, DIM>& fun);
	
  /**
   * Project a finite element function \p{src} in one finite element space to
   * a finite element function \p{des} in another finite element space, with algebriac
   * integrate accuracy order \p{algebric_accuracy}
   */
  template <class value_type, int DIM> 
    void L2Project(const FEMFunction<value_type, DIM>& src, FEMFunction<value_type, DIM>& des, 
                   Method method, int algebric_accuracy);
  /**
   * Project an analytic \p{src} in one finite element space to a finite element
   * function \p{des} in another finite element space, with algebriac integrate
   * accuracy order \p{algebric_accuracy}
   */
  template <class value_type, int DIM>
    void L2Project(value_type (*)(const double *), FEMFunction<value_type, DIM>& fun, 
                   Method method, int algebric_accuracy);
  /**
   * Project an analytic \p{src} in one finite element space to a finite element
   * function \p{des} in another finite element space, with algebriac integrate
   * accuracy order \p{algebric_accuracy}
   */
  template <class value_type, int DIM>
    void L2Project(value_type (*)(const Point<DIM>&), FEMFunction<value_type, DIM>& fun, 
                   Method method, int algebric_accuracy);
  /**
   * Project an analytic \p{src} in one finite element space to a finite element
   * function \p{des} in another finite element space, with algebriac integrate
   * accuracy order \p{algebric_accuracy}
   */
  template <class value_type, int DIM>
    void L2Project(const Function<value_type>&, FEMFunction<value_type, DIM>& fun, 
                   Method method, int algebric_accuracy);

  /**
   * Project an analytic function \p{f} of a finite element function \p{src} to a
   * finite element function \p{des} in another finite element space, which means
   * \p{des = Project( f(src) )}.
   */
  template <class value_type, int DIM>
    void L2Project(value_type (*f)(const value_type&), const FEMFunction<value_type,DIM>& src,
                   FEMFunction<value_type,DIM>& des, Method method, int algebric_accuracy);
  /**
   * Project an analytic function \p{f} of two finite element functions \p{src0} and
   * \p{src1} to a finite element function \p{des} in another finite element space,
   * which means \p{des = Project( f(src0, src1) )}. Because we can now only cope with
   * those operation with two different meshes, it's required that the three finite
   * element functions \p{src0}, \p{src1} and \p{des} should be on only two different
   * meshes.
   */
  template <class value_type, int DIM>
    void L2Project(value_type (*f)(const value_type&,const value_type&),
                   const FEMFunction<value_type,DIM>& src0, const FEMFunction<value_type,DIM>& src1,
                   FEMFunction<value_type,DIM>& des, Method method, int algebric_accuracy);

  /**
   * Discrete formation of a finite element function in the finite element space it's in.
   */
  template <class value_type, int DIM>
    void L2Discretize(const FEMFunction<value_type, DIM>& src, Vector<double>& des, int algebric_accuracy);
  /**
   * Discrete formation of a finite element function in finite element space \p{space}.
   */
  template <class value_type, int DIM>
    void L2Discretize(const FEMFunction<value_type, DIM>& src, const FEMSpace<value_type, DIM>& space,
                      Vector<double>& des, int algebric_accuracy);
  /**
   * Discrete formation of an analytic function in finite element space \p{space}.
   */
  template <class value_type, int DIM> 
    void L2Discretize(value_type (*)(const double *), const FEMSpace<value_type, DIM>& space, 
                      Vector<double>& des, int algebric_accuracy);
  /**
   * Discrete formation of an analytic function in finite element space \p{space}.
   */
  template <class value_type, int DIM>
    void L2Discretize(value_type (*)(const Point<DIM>&), const FEMSpace<value_type, DIM>& space, 
                      Vector<double>& des, int algebric_accuracy);
  /**
   * Discrete formation of an analytic function in finite element space \p{space}.
   */
  template <class value_type, int DIM>
    void L2Discretize(const Function<value_type>&, const FEMSpace<value_type, DIM>& space, 
                      Vector<double>& des, int algebric_accuracy);

  /**
   * Discrete formation of an analytic function \p{f} of a finite element function \p{src}
   * in the finite element space \p{src} in.
   */
  template <class value_type, int DIM>
    void L2Discretize(value_type (*f)(const value_type&), const FEMFunction<value_type, DIM>& src,
                      Vector<double>& des, int algebric_accuracy);
  /**
   * Discrete formation of an analytic function \p{f} of a finite element function \p{src}
   * in another finite element space \p{space}.
   */
  template <class value_type, int DIM>
    void L2Discretize(value_type (*f)(const value_type&), const FEMFunction<value_type, DIM>& src,
                      const FEMSpace<value_type, DIM>& space, Vector<double>& des, int algebric_accuracy);
  /**
   * Discrete formation of an analytic function \p{f} of two finite element functions \p{src0} and
   * \p{src1} to a finite element function \p{des} in another finite element space,
   * which means \p{des = Discrtize( f(src0, src1) )}. Because we can now only cope with
   * those operation with two different meshes, it's required that the three finite
   * element functions \p{src0}, \p{src1} and \p{des} should be on only two different
   * meshes.
   */
  template <class value_type, int DIM>
    void L2Discretize(value_type (*f)(const value_type&, const value_type&),
                      const FEMFunction<value_type, DIM>& src0, const FEMFunction<value_type, DIM>& src1,
                      const FEMSpace<value_type, DIM>& space, Vector<double>& des, int algebric_accuracy);
};

AFEPACK_CLOSE_NAMESPACE

#endif

//
// end of file
////////////////////////////////////////////////////////////////////////////
