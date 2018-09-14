/**
 * @file   BilinearOperator.templates.h
 * @author Robert Lie
 * @date   Tue Sep 13 20:12:50 2005
 * 
 * @brief  
 * 
 * 
 */

#ifndef _BilinearOperator_templates_h_
#define _BilinearOperator_templates_h_

#include "BilinearOperator.h"

AFEPACK_OPEN_NAMESPACE

#define TEMPLATE template <int DIM, class value_type0, class value_type1, int DOW, int TDIM0, int TDIM1>
#define THIS BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>

TEMPLATE
THIS::BilinearOperator()
{}



TEMPLATE
THIS::BilinearOperator(FEMSpace<value_type0,DIM,DOW,TDIM0>& sp0, 
                       FEMSpace<value_type1,DIM,DOW,TDIM1>& sp1)
{
  fem_space0 = &sp0;
  fem_space1 = &sp1;
}



TEMPLATE
THIS::~BilinearOperator()
{
  SparseMatrix<double>::clear();
}



TEMPLATE
void THIS::reinit(FEMSpace<value_type0,DIM,DOW,TDIM0>& sp0, 
                  FEMSpace<value_type1,DIM,DOW,TDIM1>& sp1)
{
  fem_space0 = &sp0;
  fem_space1 = &sp1;
}



TEMPLATE
void THIS::build()
{
  buildSparsityPattern();
  buildSparseMatrix();
}



TEMPLATE
void THIS::buildSparsityPattern()
{
  buildDofInfo();
  sparsity_pattern.reinit(nDof0(), nDof1(), nMaxCouplingDof());
  if ((void *)fem_space0 == (void *)fem_space1) {
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      getElementPattern(*the_element0, *the_element1);
      addElementPattern();
    }
  }
  else if (&(fem_space0->mesh()) == &(fem_space1->mesh())) {
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      Assert (the_element0->index() == the_element1->index(), ExcInternalError());
      getElementPattern(*the_element0, *the_element1);
      addElementPattern();
    }
  }
  else {
    Mesh<DIM,DOW>& mesh0 = fem_space0->mesh();
    Mesh<DIM,DOW>& mesh1 = fem_space1->mesh();
    RegularMesh<DIM,DOW>& regular_mesh0 = dynamic_cast<RegularMesh<DIM,DOW>&>(mesh0);
    RegularMesh<DIM,DOW>& regular_mesh1 = dynamic_cast<RegularMesh<DIM,DOW>&>(mesh1);
    IrregularMesh<DIM,DOW>& irregular_mesh0 = regular_mesh0.irregularMesh();
    IrregularMesh<DIM,DOW>& irregular_mesh1 = regular_mesh1.irregularMesh();
    Assert (&(irregular_mesh0.geometryTree()) == &(irregular_mesh1.geometryTree()), ExcInternalError());
    IrregularMeshPair<DIM,DOW> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM,DOW> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM,DOW> end_pair = mesh_pair.endActiveElementPair();
    for (;the_pair != end_pair;++ the_pair) {
      const HElement<DIM,DOW>& h_element0 = the_pair(0);
      const HElement<DIM,DOW>& h_element1 = the_pair(1);
      Element<value_type0,DIM,DOW,TDIM0>& element0 = fem_space0->element(h_element0.index);
      Element<value_type1,DIM,DOW,TDIM1>& element1 = fem_space1->element(h_element1.index);
      Assert (element0.index() == h_element0.index, ExcInternalError());
      Assert (element1.index() == h_element1.index, ExcInternalError());
      getElementPattern(element0, element1);
      addElementPattern();
    }
  }
  sparsity_pattern.compress();
}



TEMPLATE
void THIS::buildDofInfo()
{
  nDof0() = fem_space0->n_dof();
  nDof1() = fem_space1->n_dof();
  std::vector<int> n_coupling_dof(nDof0(), 0);
  if ((void *)fem_space0 == (void *)fem_space1) {
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      getElementPattern(*the_element0, *the_element1);
      int n_element_dof = elementDof0().size();
      for (int i = 0;i < n_element_dof;i ++)
	n_coupling_dof[elementDof0(i)] += n_element_dof;
    }
  }
  else if (&(fem_space0->mesh()) == &(fem_space1->mesh())) {
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      Assert (the_element0->index() == the_element1->index(), ExcInternalError());
      getElementPattern(*the_element0, *the_element1);
      int n_element_dof0 = elementDof0().size();
      int n_element_dof1 = elementDof1().size();
      for (int i = 0;i < n_element_dof0;i ++)
	n_coupling_dof[elementDof0(i)] += n_element_dof1;
    }
  }
  else {
    Mesh<DIM,DOW>& mesh0 = fem_space0->mesh();
    Mesh<DIM,DOW>& mesh1 = fem_space1->mesh();
    RegularMesh<DIM,DOW>& regular_mesh0 = dynamic_cast<RegularMesh<DIM,DOW>&>(mesh0);
    RegularMesh<DIM,DOW>& regular_mesh1 = dynamic_cast<RegularMesh<DIM,DOW>&>(mesh1);
    IrregularMesh<DIM,DOW>& irregular_mesh0 = regular_mesh0.irregularMesh();
    IrregularMesh<DIM,DOW>& irregular_mesh1 = regular_mesh1.irregularMesh();
    Assert (&(irregular_mesh0.geometryTree()) == &(irregular_mesh0.geometryTree()), ExcInternalError());
    IrregularMeshPair<DIM,DOW> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM,DOW> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM,DOW> end_pair = mesh_pair.endActiveElementPair();
    for (;the_pair != end_pair;++ the_pair) {
      const HElement<DIM,DOW>& h_element0 = the_pair(0);
      const HElement<DIM,DOW>& h_element1 = the_pair(1);
      Element<value_type0,DIM,DOW,TDIM0>& element0 = fem_space0->element(h_element0.index);
      Element<value_type1,DIM,DOW,TDIM1>& element1 = fem_space1->element(h_element1.index);
      Assert (element0.index() == h_element0.index, ExcInternalError());
      Assert (element1.index() == h_element1.index, ExcInternalError());
      getElementPattern(element0, element1);
      int n_element_dof0 = elementDof0().size();
      int n_element_dof1 = elementDof1().size();
      for (int i = 0;i < n_element_dof0;i ++)
	n_coupling_dof[elementDof0(i)] += n_element_dof1;
    }
  }
  nMaxCouplingDof() = *max_element(n_coupling_dof.begin(), n_coupling_dof.end());
  if (nMaxCouplingDof() > nDof1()) nMaxCouplingDof() = nDof1();
}



TEMPLATE
void THIS::buildSparseMatrix()
{
  SparseMatrix<double>::reinit(sparsity_pattern);
  if ((void *)fem_space0 == (void *)fem_space1) {
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      getElementPattern(*the_element0, *the_element1);
      elementMatrix().reinit(elementDof0().size(), elementDof1().size());
      getElementMatrix(*the_element0, *the_element1);
      addElementMatrix();
    }
  }
  else if (&(fem_space0->mesh()) == &(fem_space1->mesh())) {
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      Assert (the_element0->index() == the_element1->index(), ExcInternalError());
      getElementPattern(*the_element0, *the_element1);
      elementMatrix().reinit(elementDof0().size(), elementDof1().size());
      getElementMatrix(*the_element0, *the_element1);
      addElementMatrix();
    }
  }
  else {
    Mesh<DIM,DOW>& mesh0 = fem_space0->mesh();
    Mesh<DIM,DOW>& mesh1 = fem_space1->mesh();
    RegularMesh<DIM,DOW>& regular_mesh0 = dynamic_cast<RegularMesh<DIM,DOW>&>(mesh0);
    RegularMesh<DIM,DOW>& regular_mesh1 = dynamic_cast<RegularMesh<DIM,DOW>&>(mesh1);
    IrregularMesh<DIM,DOW>& irregular_mesh0 = regular_mesh0.irregularMesh();
    IrregularMesh<DIM,DOW>& irregular_mesh1 = regular_mesh1.irregularMesh();
    Assert (&(irregular_mesh0.geometryTree()) == &(irregular_mesh1.geometryTree()), ExcInternalError());
    IrregularMeshPair<DIM,DOW> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM,DOW> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM,DOW> end_pair = mesh_pair.endActiveElementPair();
    for (;the_pair != end_pair;++ the_pair) {
      const HElement<DIM,DOW>& h_element0 = the_pair(0);
      const HElement<DIM,DOW>& h_element1 = the_pair(1);
      Element<value_type0,DIM,DOW,TDIM0>& element0 = fem_space0->element(h_element0.index);
      Element<value_type1,DIM,DOW,TDIM1>& element1 = fem_space1->element(h_element1.index);
      Assert (element0.index() == h_element0.index, ExcInternalError());
      Assert (element1.index() == h_element1.index, ExcInternalError());
      getElementPattern(element0, element1);
      elementMatrix().reinit(elementDof0().size(), elementDof1().size());
      getElementMatrix(element0, element1, the_pair.state());
      addElementMatrix();
    }
  }
}



TEMPLATE
void THIS::addElementPattern()
{
  int n_element_dof0 = elementDof0().size();
  int n_element_dof1 = elementDof1().size();
  for (int i = 0;i < n_element_dof0;i ++) {
    for (int j = 0;j < n_element_dof1;j ++) {
      sparsity_pattern.add(elementDof0(i), elementDof1(j));
    }
  }
}



TEMPLATE
void THIS::addElementMatrix()
{
  int n_element_dof0 = elementDof0().size();
  int n_element_dof1 = elementDof1().size();
  for (int i = 0;i < n_element_dof0;i ++) {
    for (int j = 0;j < n_element_dof1;j ++) {
      SparseMatrix<double>::add(elementDof0(i), elementDof1(j), element_matrix(i,j));
    }
  }
}



TEMPLATE
void THIS::getElementPattern(
                             const Element<value_type0,DIM,DOW,TDIM0>& element0, 
                             const Element<value_type1,DIM,DOW,TDIM1>& element1)
{
  element_dof0 = &(element0.dof());
  element_dof1 = &(element1.dof());
}

#undef THIS
#undef TEMPLATE

#define TEMPLATE template <int DIM, class value_type0, class value_type1, int DOW, int TDIM0, int TDIM1>
#define THIS L2InnerProduct<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>

TEMPLATE
void THIS::reinit(FEMSpace<value_type0,DIM,DOW,TDIM0>& sp0, 
                  FEMSpace<value_type1,DIM,DOW,TDIM1>& sp1)
{
  BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>::reinit(sp0, sp1);	
}

TEMPLATE
void THIS::getElementMatrix(const Element<value_type0,DIM,DOW,TDIM0>& element0,
                            const Element<value_type1,DIM,DOW,TDIM1>& element1,
                            const typename ActiveElementPairIterator<DIM>::State state)
{
  const std::vector<int>& ele_dof0 = element0.dof();
  const std::vector<int>& ele_dof1 = element1.dof();
  int n_element_dof0 = ele_dof0.size();
  int n_element_dof1 = ele_dof1.size();
  if (state == ActiveElementPairIterator<DIM,DOW>::GREAT_THAN) {
    double volume = element1.templateElement().volume();
    const int& alg_acc = this->algebricAccuracy();
    const QuadratureInfo<TDIM1>& quad_info = element1.findQuadratureInfo(alg_acc);
    std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DOW> > q_point = element1.local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<value_type0> > basis_value0 = element0.basis_function_value(q_point);
    std::vector<std::vector<value_type1> > basis_value1 = element1.basis_function_value(q_point);
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      for (int j = 0;j < n_element_dof0;j ++) {
	for (int k = 0;k < n_element_dof1;k ++) {
	  this->elementMatrix(j,k) += Jxw*basis_value0[j][l]*basis_value1[k][l];
	}
      }
    }
  }
  else {
    double volume = element0.templateElement().volume();
    const int& alg_acc = this->algebricAccuracy();
    const QuadratureInfo<TDIM0>& quad_info = element0.findQuadratureInfo(alg_acc);
    std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DOW> > q_point = element0.local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<value_type0> > basis_value0 = element0.basis_function_value(q_point);
    std::vector<std::vector<value_type1> > basis_value1 = element1.basis_function_value(q_point);
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      for (int j = 0;j < n_element_dof0;j ++) {
	for (int k = 0;k < n_element_dof1;k ++) {
	  this->elementMatrix(j,k) += Jxw*basis_value0[j][l]*basis_value1[k][l];
	}
      }
    }
  }
}

#undef THIS
#undef TEMPLATE

#define TEMPLATE template <int DIM, class value_type, int DOW, int TDIM>
#define THIS MassMatrix<DIM,value_type,DOW,TDIM>

TEMPLATE
void THIS::reinit(FEMSpace<typename MassMatrix::value_type,DIM,DOW,TDIM>& sp)
{
  BilinearOperator<DIM,typename MassMatrix::value_type,typename MassMatrix::value_type,DOW,TDIM,TDIM>::reinit(sp, sp);
}



TEMPLATE
void THIS::getElementMatrix(const Element<typename MassMatrix::value_type,DIM,DOW,TDIM>& element0,
                            const Element<typename MassMatrix::value_type,DIM,DOW,TDIM>& element1,
                            const typename ActiveElementPairIterator<DIM>::State state)
{
  const std::vector<int>& ele_dof0 = element0.dof();
  const std::vector<int>& ele_dof1 = element1.dof();
  int n_element_dof0 = ele_dof0.size();
  int n_element_dof1 = ele_dof1.size();
  double volume = element0.templateElement().volume();
  const int& alg_acc = this->algebricAccuracy();
  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(alg_acc);
  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
  int n_quadrature_point = quad_info.n_quadraturePoint();
  std::vector<Point<DOW> > q_point = element0.local_to_global(quad_info.quadraturePoint());
  std::vector<std::vector<typename MassMatrix::value_type> > basis_value = element0.basis_function_value(q_point);
  for (int l = 0;l < n_quadrature_point;l ++) {
    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
    for (int j = 0;j < n_element_dof0;j ++) {
      for (int k = 0;k < n_element_dof1;k ++) {
	this->elementMatrix(j,k) += Jxw*basis_value[j][l]*basis_value[k][l];
      }
    }
  }
}

#undef THIS
#undef TEMPLATE

#define TEMPLATE template <int DIM, class value_type, int DOW, int TDIM>
#define THIS StiffMatrix<DIM,value_type,DOW,TDIM>

TEMPLATE
void THIS::reinit(FEMSpace<typename StiffMatrix::value_type,DIM,DOW,TDIM>& sp)
{
  BilinearOperator<DIM,typename StiffMatrix::value_type,typename StiffMatrix::value_type,DOW,TDIM,TDIM>::reinit(sp, sp);
}



TEMPLATE
void THIS::getElementMatrix(const Element<typename StiffMatrix::value_type,DIM,DOW,TDIM>& element0,
                            const Element<typename StiffMatrix::value_type,DIM,DOW,TDIM>& element1,
                            const typename ActiveElementPairIterator<DIM>::State state)
{
  const std::vector<int>& ele_dof0 = element0.dof();
  const std::vector<int>& ele_dof1 = element1.dof();
  int n_element_dof0 = ele_dof0.size();
  int n_element_dof1 = ele_dof1.size();
  double volume = element0.templateElement().volume();
  const int& alg_acc = this->algebricAccuracy();  
  const QuadratureInfo<TDIM>& quad_info = element0.findQuadratureInfo(alg_acc);
  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
  int n_quadrature_point = quad_info.n_quadraturePoint();
  std::vector<Point<DOW> > q_point = element0.local_to_global(quad_info.quadraturePoint());
  std::vector<std::vector<std::vector<typename StiffMatrix::value_type> > > basis_gradient = element0.basis_function_gradient(q_point);
  for (int l = 0;l < n_quadrature_point;l ++) {
    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
    for (int j = 0;j < n_element_dof0;j ++) {
      for (int k = 0;k < n_element_dof1;k ++) {
	this->elementMatrix(j,k) += Jxw*innerProduct(basis_gradient[j][l],basis_gradient[k][l]);
      }
    }
  }
}

#undef THIS
#undef TEMPLATE

AFEPACK_CLOSE_NAMESPACE

#endif

//
// end of file
/////////////////////////////////////////////////////////////////////////////
