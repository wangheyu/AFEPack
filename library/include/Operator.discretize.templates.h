/**
 * @file   Operator.discretize.templates.h
 * @author Robert Lie
 * @date   Thu Mar 30 11:42:44 2003
 * 
 * @brief  discretizations
 * 
 * 
 */

#ifndef _Operator_discretize_templates_h_
#define _Operator_discretize_templates_h_

#include "Operator.h"

AFEPACK_OPEN_NAMESPACE

template <class value_type, int DIM>
  void Operator::L2Discretize(const FEMFunction<value_type, DIM>& f0, 
                              Vector<double>& f1, 
                              int algebric_accuracy)
{
  const FEMSpace<value_type,DIM>& fem_space = f0.femSpace();
  typename FEMSpace<value_type,DIM>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  unsigned int n_dof = fem_space.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
    f1.Vector<value_type>::operator=(0.0);
  for (;the_element != end_element;++ the_element) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    std::vector<double> f0_value = f0.value(q_point, *the_element);
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      for (unsigned int j = 0;j < n_element_dof;j ++) {
	f1(element_dof[j]) += Jxw*f0_value[l]*basis_value[j][l];
      }
    }
  }
}



template <class value_type, int DIM>
void Operator::L2Discretize(const FEMFunction<value_type, DIM>& f0, 
                            const FEMSpace<value_type, DIM>& fem_space1, 
                            Vector<double>& f1, 
                            int algebric_accuracy)
{
  const FEMSpace<value_type,DIM>& fem_space0 = f0.femSpace();
  const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
  const RegularMesh<DIM>& regular_mesh1 = static_cast<const RegularMesh<DIM> &>(fem_space1.mesh());
  IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
  IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
  if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
    std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
	      << std::endl;
    Assert(false, ExcInternalError());
  }
  unsigned int n_dof = fem_space1.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
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
    if (the_pair.state() == ActiveElementPairIterator<DIM>::GREAT_THAN) {
      double volume = element1.templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = element1.findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
      std::vector<double> f0_value = f0.value(q_point, element0);
      std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof1;j ++) {
	  f1(element_dof1[j]) += Jxw*f0_value[l]*basis_value1[j][l];
	}
      }
    }
    else {
      double volume = element0.templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
      std::vector<double> f0_value = f0.value(q_point, element0);
      std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof1;j ++) {
	  f1(element_dof1[j]) += Jxw*f0_value[l]*basis_value1[j][l];
	}
      }
    }
  }
}



template <class value_type, int DIM> 
void Operator::L2Discretize(value_type (*f0)(const double *), 
                            const FEMSpace<value_type, DIM>& fem_space, 
                            Vector<double>& f1, 
                            int algebric_accuracy)
{
  unsigned int n_dof = fem_space.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
    f1.Vector<value_type>::operator=(0.0);
  typename FEMSpace<value_type,DIM>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (int l = 0;l < n_quadrature_point;l ++) {
      value_type f0_value = (*f0)(q_point[l]);
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      for (unsigned int j = 0;j < n_element_dof;j ++) {
	f1(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
      }
    }
  }
}



template <class value_type, int DIM> 
void Operator::L2Discretize(value_type (*f0)(const Point<DIM>&), 
                            const FEMSpace<value_type, DIM>& fem_space, 
                            Vector<double>& f1, 
                            int algebric_accuracy)
{
  unsigned int n_dof = fem_space.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
    f1.Vector<value_type>::operator=(0.0);
  typename FEMSpace<value_type,DIM>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (int l = 0;l < n_quadrature_point;l ++) {
      value_type f0_value = (*f0)(q_point[l]);
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      for (unsigned int j = 0;j < n_element_dof;j ++) {
	f1(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
      }
    }
  }
}




template <class value_type, int DIM> 
void Operator::L2Discretize(const Function<value_type>& f0, 
                            const FEMSpace<value_type, DIM>& fem_space, 
                            Vector<double>& f1, 
                            int algebric_accuracy)
{
  unsigned int n_dof = fem_space.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
    f1.Vector<value_type>::operator=(0.0);
  typename FEMSpace<value_type,DIM>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (int l = 0;l < n_quadrature_point;l ++) {
      value_type f0_value = f0.value(q_point[l]);
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      for (unsigned int j = 0;j < n_element_dof;j ++) {
	f1(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
      }
    }
  }
}




template <class value_type, int DIM>
void Operator::L2Discretize(value_type (*f)(const value_type&), 
                            const FEMFunction<value_type, DIM>& f0,
                            Vector<double>& f1, 
                            int algebric_accuracy)
{
  const FEMSpace<value_type,DIM>& fem_space = f0.femSpace();
  typename FEMSpace<value_type,DIM>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  unsigned int n_dof = fem_space.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
    f1.Vector<value_type>::operator=(0.0);
  for (;the_element != end_element;++ the_element) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
    std::vector<double> f0_value = f0.value(q_point, *the_element);
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      value_type f_value = (*f)(f0_value[l]);
      for (unsigned int j = 0;j < n_element_dof;j ++) {
	f1(element_dof[j]) += Jxw*f_value*basis_value[j][l];
      }
    }
  }
}



template <class value_type, int DIM>
void Operator::L2Discretize(value_type (*f)(const value_type&), 
                            const FEMFunction<value_type, DIM>& f0,
                            const FEMSpace<value_type, DIM>& fem_space1, 
                            Vector<double>& f1, 
                            int algebric_accuracy)
{
  const FEMSpace<value_type,DIM>& fem_space0 = f0.femSpace();
  unsigned int n_dof = fem_space1.n_dof();
  if (f1.size() != n_dof)
    f1.reinit(n_dof);
  else
    f1.Vector<value_type>::operator=(0.0);
  if (&fem_space0 == &fem_space1 || &(fem_space0.mesh()) == &(fem_space1.mesh())) {
    typename FEMSpace<value_type,DIM>::ConstElementIterator 
      the_element = fem_space0.beginElement(),
      end_element = fem_space0.endElement();
    typename FEMSpace<value_type,DIM>::ConstElementIterator 
      the_element1 = fem_space1.beginElement();
    for (;the_element != end_element;++ the_element, ++ the_element1) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<double> > basis_value = the_element1->basis_function_value(q_point);
      std::vector<double> f0_value = f0.value(q_point, *the_element);
      const std::vector<int>& element_dof = the_element1->dof();
      unsigned int n_element_dof = element_dof.size();
      for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	value_type f_value = (*f)(f0_value[l]);
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  f1(element_dof[j]) += Jxw*f_value*basis_value[j][l];
	}
      }
    }
    return;
  }
  const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
  const RegularMesh<DIM>& regular_mesh1 = static_cast<const RegularMesh<DIM> &>(fem_space1.mesh());
  IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
  IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
  if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
    std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
	      << std::endl;
    Assert(false, ExcInternalError());
  }
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
    if (the_pair.state() == ActiveElementPairIterator<DIM>::GREAT_THAN) {
      double volume = element1.templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = element1.findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
      std::vector<double> f0_value = f0.value(q_point, element0);
      std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	value_type f_value = (*f)(f0_value[l]);
	for (unsigned int j = 0;j < n_element_dof1;j ++) {
	  f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	}
      }
    }
    else {
      double volume = element0.templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
      std::vector<double> f0_value = f0.value(q_point, element0);
      std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	value_type f_value = (*f)(f0_value[l]);
	for (unsigned int j = 0;j < n_element_dof1;j ++) {
	  f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	}
      }
    }
  }
}



template <class value_type, int DIM>
void Operator::L2Discretize(value_type (*f)(const value_type&, const value_type&),
                            const FEMFunction<value_type, DIM>& f0, 
                            const FEMFunction<value_type, DIM>& f1,
                            const FEMSpace<value_type, DIM>& fem_space2, 
                            Vector<double>& f2, 
                            int algebric_accuracy)
{
  const FEMSpace<value_type,DIM>& fem_space0 = f0.femSpace();
  const FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
  if ((&fem_space0 != &fem_space1) &&
      (&fem_space0 != &fem_space2) &&
      (&fem_space1 != &fem_space2)) {
    std::cerr << "The three FEM functions are on three different finite element spaces."
	      << std::endl;
    abort();
  }
  unsigned int n_dof = fem_space2.n_dof();
  if (f2.size() != n_dof)
    f2.reinit(n_dof);
  else
    f2.Vector<value_type>::operator=(0.0);
  if ((&fem_space0 == &fem_space1) && (&fem_space0 == &fem_space2)) {
    typename FEMSpace<value_type,DIM>::ConstElementIterator 
      the_element = fem_space0.beginElement(),
      end_element = fem_space0.endElement();
    for (;the_element != end_element;++ the_element) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
      std::vector<double> f0_value = f0.value(q_point, *the_element);
      std::vector<double> f1_value = f1.value(q_point, *the_element);
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      for (int l = 0;l < n_quadrature_point;l ++) {
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	value_type f_value = (*f)(f0_value[l], f1_value[l]);
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  f2(element_dof[j]) += Jxw*f_value*basis_value[j][l];
	}
      }
    }
  }
  else if (&fem_space0 == &fem_space1) {
    const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
    const RegularMesh<DIM>& regular_mesh1 = static_cast<const RegularMesh<DIM> &>(fem_space2.mesh());
    IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
    IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
    if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
      std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		<< std::endl;
      Assert(false, ExcInternalError());
    }
    IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
    for (;the_pair != end_pair;++ the_pair) {
      const HElement<DIM>& h_element0 = the_pair(0);
      const HElement<DIM>& h_element1 = the_pair(1);
      const Element<value_type,DIM>& element0 = fem_space0.element(h_element0.index);
      const Element<value_type,DIM>& element1 = fem_space2.element(h_element1.index);
      Assert(element0.index() == h_element0.index, ExcInternalError());
      Assert(element1.index() == h_element1.index, ExcInternalError());
      const std::vector<int>& element_dof1 = element1.dof();
      unsigned int n_element_dof1 = element_dof1.size();
      if (the_pair.state() == ActiveElementPairIterator<DIM>::GREAT_THAN) {
	double volume = element1.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element1.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	std::vector<double> f0_value = f0.value(q_point, element0);
	std::vector<double> f1_value = f1.value(q_point, element0);
	std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  value_type f_value = (*f)(f0_value[l], f1_value[l]);
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f2(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	  }
	}
      }
      else {
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<double> f0_value = f0.value(q_point, element0);
	std::vector<double> f1_value = f1.value(q_point, element0);
	std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  value_type f_value = (*f)(f0_value[l], f1_value[l]);
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f2(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	  }
	}
      }
    }
  }
  else if (&fem_space0 == &fem_space2) {
    const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space1.mesh());
    const RegularMesh<DIM>& regular_mesh1 = static_cast<const RegularMesh<DIM> &>(fem_space2.mesh());
    IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
    IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
    if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
      std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		<< std::endl;
      Assert(false, ExcInternalError());
    }
    IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
    for (;the_pair != end_pair;++ the_pair) {
      const HElement<DIM>& h_element0 = the_pair(0);
      const HElement<DIM>& h_element1 = the_pair(1);
      const Element<value_type,DIM>& element0 = fem_space1.element(h_element0.index);
      const Element<value_type,DIM>& element1 = fem_space2.element(h_element1.index);
      Assert(element0.index() == h_element0.index, ExcInternalError());
      Assert(element1.index() == h_element1.index, ExcInternalError());
      const std::vector<int>& element_dof1 = element1.dof();
      unsigned int n_element_dof1 = element_dof1.size();
      if (the_pair.state() == ActiveElementPairIterator<DIM>::GREAT_THAN) {
	double volume = element1.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element1.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	std::vector<double> f0_value = f0.value(q_point, element1);
	std::vector<double> f1_value = f1.value(q_point, element0);
	std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  value_type f_value = (*f)(f0_value[l], f1_value[l]);
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f2(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	  }
	}
      }
      else {
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<double> f0_value = f0.value(q_point, element1);
	std::vector<double> f1_value = f1.value(q_point, element0);
	std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  value_type f_value = (*f)(f0_value[l], f1_value[l]);
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f2(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	  }
	}
      }
    }
  }
  else if (&fem_space1 == &fem_space2) {
    const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
    const RegularMesh<DIM>& regular_mesh1 = static_cast<const RegularMesh<DIM> &>(fem_space2.mesh());
    IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
    IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
    if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
      std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		<< std::endl;
      Assert(false, ExcInternalError());
    }
    IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
    for (;the_pair != end_pair;++ the_pair) {
      const HElement<DIM>& h_element0 = the_pair(0);
      const HElement<DIM>& h_element1 = the_pair(1);
      const Element<value_type,DIM>& element0 = fem_space0.element(h_element0.index);
      const Element<value_type,DIM>& element1 = fem_space2.element(h_element1.index);
      Assert(element0.index() == h_element0.index, ExcInternalError());
      Assert(element1.index() == h_element1.index, ExcInternalError());
      const std::vector<int>& element_dof1 = element1.dof();
      unsigned int n_element_dof1 = element_dof1.size();
      if (the_pair.state() == ActiveElementPairIterator<DIM>::GREAT_THAN) {
	double volume = element1.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element1.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	std::vector<double> f0_value = f0.value(q_point, element0);
	std::vector<double> f1_value = f1.value(q_point, element1);
	std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  value_type f_value = (*f)(f0_value[l], f1_value[l]);
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f2(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	  }
	}
      }
      else {
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<double> f0_value = f0.value(q_point, element0);
	std::vector<double> f1_value = f1.value(q_point, element1);
	std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  value_type f_value = (*f)(f0_value[l], f1_value[l]);
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f2(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	  }
	}
      }
    }
  }
  else { // something must be wrong
    Assert(false, ExcInternalError());
  }
}

AFEPACK_CLOSE_NAMESPACE

#endif

/**
 * end of file
 * 
 */
