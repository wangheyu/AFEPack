/**
 * @file   Operator.discretize.cpp
 * @author Robert Lie
 * @date   Thu Mar 30 11:47:27 2006
 * 
 * @brief  discretizations
 * 
 * 
 */

#include "Operator.discretize.templates.h"

AFEPACK_OPEN_NAMESPACE

#define value_type double
#define DIM 2
namespace Operator {
  template void L2Discretize(const FEMFunction<value_type, DIM>&, 
                             Vector<double>&, 
                             int);
  template void L2Discretize(const FEMFunction<value_type, DIM>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const double *), 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const Point<DIM>&), 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(const Function<value_type>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);

  template void L2Discretize(value_type (*)(const value_type&), 
                             const FEMFunction<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const value_type&), 
                             const FEMFunction<value_type, DIM>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const value_type&, const value_type&), 
                             const FEMFunction<value_type, DIM>&, 
                             const FEMFunction<value_type, DIM>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
};

#undef DIM

#define DIM 3
namespace Operator {
  template void L2Discretize(const FEMFunction<value_type, DIM>&, 
                             Vector<double>&, 
                             int);
  template void L2Discretize(const FEMFunction<value_type, DIM>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const double *), 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const Point<DIM>&), 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(const Function<value_type>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);

  template void L2Discretize(value_type (*)(const value_type&), 
                             const FEMFunction<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const value_type&), 
                             const FEMFunction<value_type, DIM>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
  template void L2Discretize(value_type (*)(const value_type&, const value_type&), 
                             const FEMFunction<value_type, DIM>&, 
                             const FEMFunction<value_type, DIM>&, 
                             const FEMSpace<value_type, DIM>&, 
                             Vector<double>&, int);
}
#undef DIM
#undef value_type

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

