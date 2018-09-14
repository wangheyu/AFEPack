/**
 * @file   Operator.interpolate.templates.h
 * @author Robert Lie
 * @date   Thu Mar 30 11:39:19 2003
 * 
 * @brief  interpolations
 * 
 * 
 */

#ifndef _Operator_interpolate_templates_h_
#define _Operator_interpolate_templates_h_

#include "Operator.h"

AFEPACK_OPEN_NAMESPACE

template <class value_type, int DIM> 
void Operator::L2Interpolate(const FEMFunction<value_type, DIM>& f0, FEMFunction<value_type, DIM>& f1)
{
  const FEMSpace<value_type,DIM>& fem_space0 = f0.femSpace();
  FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
  const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
  const RegularMesh<DIM>& regular_mesh1 = static_cast<const RegularMesh<DIM> &>(fem_space1.mesh());
  IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
  IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
  if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
    std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
	      << std::endl;
    Assert(false, ExcInternalError());
  }
  std::vector<bool> flag(f1.size(), false);
  f1.Vector<value_type>::operator=(0.0);
  IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
  ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
  ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
  for (;the_pair != end_pair;++ the_pair) {
    const HElement<DIM>& h_element0 = the_pair(0);
    const HElement<DIM>& h_element1 = the_pair(1);
    const Element<value_type,DIM>& element0 = fem_space0.element(h_element0.index);
    const Element<value_type,DIM>& element1 = fem_space1.element(h_element1.index);
    Assert(element0.index() == h_element0.index, ExcInternalError());
    Assert(element1.index() == h_element1.index, ExcInternalError());
    const std::vector<int>& element_dof1 = element1.dof();
    unsigned int n_element_dof1 = element_dof1.size();
    if (the_pair.state() != ActiveElementPairIterator<DIM>::LESS_THAN) {
      for (unsigned int i = 0;i < n_element_dof1;i ++) {
	unsigned int j = element_dof1[i];
	const Point<DIM>& interp_point = fem_space1.dofInfo(j).interp_point;
	f1(j) = f0.value(interp_point, element0);
	flag[j] = true;
      }
    }
    else {
      for (unsigned int i = 0;i < n_element_dof1;i ++) {
	unsigned int j = element_dof1[i];
	if (flag[j]) continue;
	const Point<DIM>& interp_point = fem_space1.dofInfo(j).interp_point;
	if (h_element0.isIncludePoint(interp_point)) {
	  f1(j) = f0.value(interp_point, element0);
	  flag[j] = true;
	}
      }
    }
  }
}



template <class value_type, int DIM> 
void Operator::L2Interpolate(value_type (*f0)(const double *), FEMFunction<value_type, DIM>& f1)
{
  FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
  typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
  typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      unsigned int j = element_dof[i];
      const Point<DIM>& interp_point = fem_space.dofInfo(j).interp_point;
      f1(j) = f0(interp_point);
    }
  }
}



template <class value_type, int DIM> 
void Operator::L2Interpolate(value_type (*f0)(const Point<DIM>&), FEMFunction<value_type, DIM>& f1)
{
  FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
  typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
  typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      unsigned int j = element_dof[i];
      const Point<DIM>& interp_point = fem_space.dofInfo(j).interp_point;
      f1(j) = f0(interp_point);
    }
  }
}



template <class value_type, int DIM> 
void Operator::L2Interpolate(const Function<value_type>& f0, FEMFunction<value_type, DIM>& f1)
{
  FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
  typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
  typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      unsigned int j = element_dof[i];
      const Point<DIM>& interp_point = fem_space.dofInfo(j).interp_point;
      f1(j) = f0.value(interp_point);
    }
  }
}

AFEPACK_CLOSE_NAMESPACE

#endif

/**
 * end of file
 * 
 */
