////////////////////////////////////////////////////////////////////////////
//  BilinearOperator.h : by R.Lie
//

#ifndef _BilinearOperator_h_
#define _BilinearOperator_h_

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

#include <base/exceptions.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include "Geometry.h"
#include "TemplateElement.h"
#include "FEMSpace.h"
#include "HGeometry.h"

AFEPACK_OPEN_NAMESPACE

/**
 * class BilinearOperator is the discretized version of the bilinear
 * form on two linear space
 * \f[
 *   a(u,v), \qquad u \in U, v \in V
 * \f]
 * thus for \f$ U \f$ and \f$ V \f$ to be finite dimensional space
 * with basis function as \f$ \phi^i, 1\leq i\leq M \f$ and \f$ \psi^j,
 * 1\leq j \leq N \f$, respectively, the bilinear form can be given
 * as a matrix with its entry as
 * \f[
 *   a^{ij} = a(\phi^i, \psi^j)
 * \f]
 * The linear finite dimensional spaces here are two finite element
 * space and the resulted matrix is a sparse matrix. The method <i>
 * getElementMatrix</i> is a pure virtual function that you must specify
 * this method in derived class to contruct an object. Typical code
 * to use this class is as:
 * <pre>
 * // definition of the derived class
 * class a_bilinear_operator : 
 *   public BilinearOperator<DIM, double>
 * {
 *   ... ...
 *   virtual void getElementMatrix(const Element<double,DIM>&,
 *                                 const Element<double,DIM>&);
 *   ... ...
 * };
 *
 * a_bilinear_operator::getElementMatrix(const Element<double,DIM>&,
 *                                       const Element<double,DIM>&)
 * {
 *   // implementation code of this bilinear operator
 *   ... ...
 * };
 *
 * // use this bilinear operator
 * void some_function()
 * {
 *   ... ...
 *   a_bilinear_operator op(fem_space0, fem_space1);
 *   a_bilinear_operator.algebricAccuracy() = 3;
 *   a_bilinear_operator.build();
 *   ... ...
 * };
 * </pre>
 *
 */
template <int DIM, class value_type0, class value_type1=value_type0,
  int DOW=DIM, int TDIM0=DIM, int TDIM1 = DIM>
  class BilinearOperator : public SparseMatrix<double>
  {
  private:
  SparsityPattern					sparsity_pattern;
  FEMSpace<value_type0,DIM,DOW,TDIM0> *		fem_space0;
  FEMSpace<value_type1,DIM,DOW,TDIM1> *		fem_space1;
  int						n_dof0;
  int						n_dof1;
  int						n_max_coupling_dof;
  const std::vector<int> *			element_dof0;
  const std::vector<int> *			element_dof1;
  FullMatrix<double>				element_matrix;
  int						algebric_accuracy;
  public:
  BilinearOperator();
  BilinearOperator(FEMSpace<value_type0,DIM,DOW,TDIM0>&, FEMSpace<value_type1,DIM,DOW,TDIM1>&);
  virtual ~BilinearOperator();
  public:
  void reinit(FEMSpace<value_type0,DIM,DOW,TDIM0>&, FEMSpace<value_type1,DIM,DOW,TDIM1>&);
  const FEMSpace<value_type0,DIM,DOW,TDIM0>& FEMSpace0() const {return *fem_space0;};
  const FEMSpace<value_type1,DIM,DOW,TDIM1>& FEMSpace1() const {return *fem_space1;};
  FEMSpace<value_type0,DIM,DOW,TDIM0> * FEMSpace0() {return fem_space0;};
  FEMSpace<value_type1,DIM,DOW,TDIM1> * FEMSpace1() {return fem_space1;};
  const SparsityPattern& getSparsityPattern() const {return sparsity_pattern;};
  SparsityPattern& getSparsityPattern() {return sparsity_pattern;};
  public:
  const int& nDof0() const {return n_dof0;};
  const int& nDof1() const {return n_dof1;};
  const int& nMaxCouplingDof() const {return n_max_coupling_dof;};
  int& nDof0() {return n_dof0;};
  int& nDof1() {return n_dof1;};
  int& nMaxCouplingDof() {return n_max_coupling_dof;};

  const std::vector<int>& elementDof0() const {return *element_dof0;};
  const std::vector<int>& elementDof1() const {return *element_dof1;};
  const int& elementDof0(const int& i) const {return (*element_dof0)[i];};
  const int& elementDof1(const int& i) const {return (*element_dof1)[i];};

  const FullMatrix<double>& elementMatrix() const {return element_matrix;};
  const double elementMatrix(const int& i, const int& j) const {return element_matrix(i,j);};
  FullMatrix<double>& elementMatrix() {return element_matrix;};
  double& elementMatrix(const int& i, const int& j) {return element_matrix(i,j);};
  void buildDofInfo();
  void addElementPattern();
  void addElementMatrix();
  void getElementPattern(const Element<value_type0,DIM,DOW,TDIM0>&, const Element<value_type1,DIM,DOW,TDIM1>&);
  public:
  const int& algebricAccuracy() const {return algebric_accuracy;};
  int& algebricAccuracy() {return algebric_accuracy;};

  virtual void getElementMatrix(const Element<value_type0,DIM,DOW,TDIM0>&, 
                                const Element<value_type1,DIM,DOW,TDIM1>&,
                                const typename ActiveElementPairIterator<DIM>::State state = ActiveElementPairIterator<DIM>::EQUAL) = 0;
  virtual void build();
  virtual void buildSparsityPattern();
  virtual void buildSparseMatrix();
  };

/**
 * class L2InnerProduct is a special bilinear operator standing as the \f$ L^2 \f$
 * inner product of two \f$ L^2 \f$ integratable functions.
 * 
 */
template <int DIM, class value_type0, class value_type1=value_type0,
  int DOW=DIM, int TDIM0=DIM, int TDIM1 = DIM>
  class L2InnerProduct : public BilinearOperator<DIM, value_type0, value_type1, DOW, TDIM0, TDIM1>
  {
  public:
  L2InnerProduct() {};
  L2InnerProduct(FEMSpace<value_type0,DIM,DOW,TDIM0>& sp0, FEMSpace<value_type1,DIM,DOW,TDIM1>& sp1) :
  BilinearOperator<DIM, value_type0, value_type1, DOW, TDIM0, TDIM1>(sp0, sp1) {};
  virtual ~L2InnerProduct() {};
  public:
  void reinit(FEMSpace<value_type0,DIM,DOW,TDIM0>&, FEMSpace<value_type1,DIM,DOW,TDIM1>&);
  virtual void getElementMatrix(const Element<value_type0,DIM,DOW,TDIM0>&, 
                                const Element<value_type1,DIM,DOW,TDIM1>&,
                                const typename ActiveElementPairIterator<DIM>::State state = ActiveElementPairIterator<DIM>::EQUAL);
  };

template <int DIM, class value_type, int DOW=DIM, int TDIM=DIM>
  class MassMatrix : public BilinearOperator<DIM, value_type, value_type, DOW, TDIM, TDIM>
  {
  public:
  MassMatrix() {};
  MassMatrix(FEMSpace<typename MassMatrix::value_type,DIM,DOW,TDIM>& sp) :
  BilinearOperator<DIM, typename MassMatrix::value_type, 
  typename MassMatrix::value_type, DOW, TDIM, TDIM>(sp, sp) {};
  virtual ~MassMatrix() {};
  public:
  void reinit(FEMSpace<typename MassMatrix::value_type,DIM,DOW,TDIM>& sp);
  virtual void getElementMatrix(const Element<typename MassMatrix::value_type,DIM,DOW,TDIM>& e0, 
                                const Element<typename MassMatrix::value_type,DIM,DOW,TDIM>& e1,
                                const typename ActiveElementPairIterator<DIM>::State state = ActiveElementPairIterator<DIM>::EQUAL);
  };

template <int DIM, class value_type, int DOW=DIM, int TDIM=DIM>
  class StiffMatrix : public BilinearOperator<DIM, value_type, value_type, DOW, TDIM, TDIM>
  {
  public:
  StiffMatrix() {};
  StiffMatrix(FEMSpace<typename StiffMatrix::value_type,DIM,DOW,TDIM>& sp) :
  BilinearOperator<DIM, typename StiffMatrix::value_type, 
  typename StiffMatrix::value_type,DOW,TDIM,TDIM>(sp, sp) {};
  virtual ~StiffMatrix() {};
  public:
  void reinit(FEMSpace<typename StiffMatrix::value_type,DIM,DOW,TDIM>& sp);
  virtual void getElementMatrix(const Element<typename StiffMatrix::value_type,DIM,DOW,TDIM>&, 
                                const Element<typename StiffMatrix::value_type,DIM,DOW,TDIM>&,
                                const typename ActiveElementPairIterator<DIM>::State state = ActiveElementPairIterator<DIM>::EQUAL);
  };

AFEPACK_CLOSE_NAMESPACE

#endif

//
// end of file
/////////////////////////////////////////////////////////////////////////////
