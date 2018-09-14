/**
 * @file   Operator.project.templates.h
 * @author Robert Lie
 * @date   Thu Mar 30 11:41:05 2003
 * 
 * @brief  projections
 * 
 * 
 */

#ifndef _Operator_project_templates_h_
#define _Operator_project_templates_h_

#include "Operator.h"

AFEPACK_OPEN_NAMESPACE

template <class value_type, int DIM>
void Operator::L2Project(const FEMFunction<value_type, DIM>& f0, 
			 FEMFunction<value_type, DIM>& f1, 
			 Method method, 
                         int algebric_accuracy)
{
  if (&(f0.femSpace().mesh()) == &(f1.femSpace().mesh())) {
    switch (method) {
    case MASS_ACCUMULATION: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      Vector<double> mass_accumulation(n_dof);
      f1.Vector<value_type>::operator=(0.0);
      typename FEMSpace<value_type,DIM>::ConstElementIterator 
        the_element0 = f0.femSpace().beginElement(),
        the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++, the_element0 ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, *the_element0);
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    f1(element_dof[j]) += Jxw*f0_value[l]*basis_value[j][l];
	    mass_accumulation(element_dof[j]) += Jxw*basis_value[j][l];
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) 
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      f1.Vector<value_type>::operator=(0.0);
      unsigned int n_dof = fem_space.n_dof();
      MassMatrix<DIM, value_type> mass_matrix(fem_space);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      Vector<double> rhs(n_dof);
      typename FEMSpace<value_type,DIM>::ConstElementIterator 
        the_element0 = f0.femSpace().beginElement(),
        the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++, the_element0 ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, *the_element0);
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    rhs(element_dof[j]) += Jxw*f0_value[l]*basis_value[j][l];
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      typename FEMSpace<value_type,DIM>::ConstElementIterator 
        the_element0 = f0.femSpace().beginElement(),
        the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++, the_element0 ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	FullMatrix<double> local_mass_matrix(n_element_dof, n_element_dof);
	Vector<double> local_rhs(n_element_dof);
	Vector<double> local_f1(n_element_dof);
	unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, *the_element0);
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    for (unsigned int k = 0;k < n_element_dof;k ++) {
	      local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	    }
	    local_rhs(j) += Jxw*f0_value[l]*basis_value[j][l];
	  }
	}
	local_mass_matrix.gauss_jordan();
	local_mass_matrix.vmult(local_f1, local_rhs);
	for (unsigned int i = 0;i < n_element_dof;i ++) {
	  f1(element_dof[i]) += local_f1(i);
	  counter[element_dof[i]] ++;
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
    }
    }
    return;
  }
  // otherwise, the two FEM functions are on different meshes
  switch (method) {
  case MASS_ACCUMULATION: {
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
    unsigned int n_dof = fem_space1.n_dof();
    Vector<double> mass_accumulation(n_dof);
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
	std::vector<value_type> f0_value = f0.value(q_point, element0);
	std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f1(element_dof1[j]) += Jxw*f0_value[l]*basis_value1[j][l];
	    mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	  }
	}
      }
      else {
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, element0);
	std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    f1(element_dof1[j]) += Jxw*f0_value[l]*basis_value1[j][l];
	    mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	  }
	}
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) 
      f1(i) /= mass_accumulation(i);
    break;
  }
  case LEAST_SQUARE: {
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
    f1.Vector<value_type>::operator=(0.0);
    MassMatrix<DIM, value_type> mass_matrix(fem_space1);
    mass_matrix.algebricAccuracy() = algebric_accuracy;
    mass_matrix.build();
    unsigned int n_dof = fem_space1.n_dof();
    Vector<double> rhs(n_dof);
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
	std::vector<value_type> f0_value = f0.value(q_point, element0);
	std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    rhs(element_dof1[j]) += Jxw*f0_value[l]*basis_value1[j][l];
	  }
	}
      }
      else {
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, element0);
	std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    rhs(element_dof1[j]) += Jxw*f0_value[l]*basis_value1[j][l];
	  }
	}
      }
    }
    AMGSolver solver(mass_matrix);
    solver.solve(f1, rhs);
    break;
  }
  case LOCAL_LEAST_SQUARE: {
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
    unsigned int n_dof = fem_space1.n_dof();
    std::vector<int> counter(n_dof, 0);
    f1.Vector<value_type>::operator=(0.0);
    IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
    ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
    ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
    int the_element_index;
    unsigned int n_the_element_dof;
    { // only used to make the variables local
      const HElement<DIM>& the_h_element = the_pair(1);
      the_element_index = the_h_element.index;
      const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
      Assert(the_element.index() == the_element_index, ExcInternalError());
      n_the_element_dof = the_element.n_dof();
    }
    Vector<double> local_rhs(n_the_element_dof);
    for (;the_pair != end_pair;) {
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
	unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, element0);
	std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    local_rhs(j) += Jxw*f0_value[l]*basis_value1[j][l];
	  }
	}
      }
      else {
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<value_type> f0_value = f0.value(q_point, element0);
	std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof1;j ++) {
	    local_rhs(j) += Jxw*f0_value[l]*basis_value1[j][l];
	  }
	}
      }
      ++ the_pair;
      if (the_pair == end_pair || the_pair(1).index != the_element_index) {
	const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	const std::vector<int>& the_element_dof = the_element.dof();
	FullMatrix<double> local_mass_matrix(n_the_element_dof, n_the_element_dof);
	Vector<value_type> local_f1(n_the_element_dof);
	double volume = the_element.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element.findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element.local_to_global_jacobian(quad_info.quadraturePoint());
	unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element.local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element.basis_function_value(q_point);
	for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_the_element_dof;j ++) {
	    for (unsigned int k = 0;k < n_the_element_dof;k ++) {
	      local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	    }
	  }
	}
	local_mass_matrix.gauss_jordan();
	local_mass_matrix.vmult(local_f1, local_rhs);
	for (unsigned int i = 0;i < n_the_element_dof;i ++) {
	  f1(the_element_dof[i]) += local_f1(i);
	  counter[the_element_dof[i]] ++;
	}
	if (the_pair != end_pair) {
	  the_element_index = the_pair(1).index;
	  Element<value_type,DIM>& the_new_element = fem_space1.element(the_element_index);
	  n_the_element_dof = the_new_element.n_dof();
	  local_rhs.reinit(n_the_element_dof);
	}
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
    break;
  }
  default: {
    Assert(false, ExcNotImplemented());
  }
  }
}



template <class value_type, int DIM> 
  void Operator::L2Project(value_type (*f0)(const double *), 
                           FEMFunction<value_type, DIM>& f1, 
                           Method method, 
                           int algebric_accuracy)
{
  switch (method) {
  case MASS_ACCUMULATION: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    unsigned int n_dof = fem_space.n_dof();
    Vector<double> mass_accumulation(n_dof);
    f1.Vector<value_type>::operator=(0.0);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = (*f0)(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  f1(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
	  mass_accumulation(element_dof[j]) += Jxw*basis_value[j][l];
	}
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) 
      f1(i) /= mass_accumulation(i);
    break;
  }
  case LEAST_SQUARE: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    f1.Vector<value_type>::operator=(0.0);
    unsigned int n_dof = fem_space.n_dof();
    MassMatrix<DIM, value_type> mass_matrix(fem_space);
    mass_matrix.algebricAccuracy() = algebric_accuracy;
    mass_matrix.build();
    Vector<double> rhs(n_dof);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = (*f0)(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  rhs(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
	}
      }
    }
    AMGSolver solver(mass_matrix);
    solver.solve(f1, rhs);
    break;
  }
  case LOCAL_LEAST_SQUARE: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    unsigned int n_dof = fem_space.n_dof();
    std::vector<int> counter(n_dof, 0);
    f1.Vector<value_type>::operator=(0.0);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      FullMatrix<double> local_mass_matrix(n_element_dof, n_element_dof);
      Vector<double> local_rhs(n_element_dof);
      Vector<double> local_f1(n_element_dof);
      unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = (*f0)(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  for (unsigned int k = 0;k < n_element_dof;k ++) {
	    local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	  }
	  local_rhs(j) += Jxw*f0_value*basis_value[j][l];
	}
      }
      local_mass_matrix.gauss_jordan();
      local_mass_matrix.vmult(local_f1, local_rhs);
      for (unsigned int i = 0;i < n_element_dof;i ++) {
	f1(element_dof[i]) += local_f1(i);
	counter[element_dof[i]] ++;
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
    break;
  }
  default: {
    Assert(false, ExcNotImplemented());
  }
  }
}



template <class value_type, int DIM> 
  void Operator::L2Project(value_type (*f0)(const Point<DIM>&), 
                           FEMFunction<value_type, DIM>& f1, 
                           Method method, 
                           int algebric_accuracy)
{
  switch (method) {
  case MASS_ACCUMULATION: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    unsigned int n_dof = fem_space.n_dof();
    Vector<double> mass_accumulation(n_dof);
    f1.Vector<value_type>::operator=(0.0);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = (*f0)(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  f1(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
	  mass_accumulation(element_dof[j]) += Jxw*basis_value[j][l];
	}
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) 
      f1(i) /= mass_accumulation(i);
    break;
  }
  case LEAST_SQUARE: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    f1.Vector<value_type>::operator=(0.0);
    unsigned int n_dof = fem_space.n_dof();
    MassMatrix<DIM, value_type> mass_matrix(fem_space);
    mass_matrix.algebricAccuracy() = algebric_accuracy;
    mass_matrix.build();
    Vector<double> rhs(n_dof);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = (*f0)(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  rhs(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
	}
      }
    }
    AMGSolver solver(mass_matrix);
    solver.solve(f1, rhs);
    break;
  }
  case LOCAL_LEAST_SQUARE: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    unsigned int n_dof = fem_space.n_dof();
    std::vector<int> counter(n_dof, 0);
    f1.Vector<value_type>::operator=(0.0);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      FullMatrix<double> local_mass_matrix(n_element_dof, n_element_dof);
      Vector<double> local_rhs(n_element_dof);
      Vector<double> local_f1(n_element_dof);
      unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = (*f0)(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  for (unsigned int k = 0;k < n_element_dof;k ++) {
	    local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	  }
	  local_rhs(j) += Jxw*f0_value*basis_value[j][l];
	}
      }
      local_mass_matrix.gauss_jordan();
      local_mass_matrix.vmult(local_f1, local_rhs);
      for (unsigned int i = 0;i < n_element_dof;i ++) {
	f1(element_dof[i]) += local_f1(i);
	counter[element_dof[i]] ++;
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
    break;
  }
  default: {
    Assert(false, ExcNotImplemented());
  }
  }
}



template <class value_type, int DIM> 
  void Operator::L2Project(const Function<value_type>& f0, 
                           FEMFunction<value_type, DIM>& f1, 
                           Method method, 
                           int algebric_accuracy)
{
  switch (method) {
  case MASS_ACCUMULATION: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    unsigned int n_dof = fem_space.n_dof();
    Vector<double> mass_accumulation(n_dof);
    f1.Vector<value_type>::operator=(0.0);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = f0.value(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  f1(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
	  mass_accumulation(element_dof[j]) += Jxw*basis_value[j][l];
	}
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) 
      f1(i) /= mass_accumulation(i);
    break;
  }
  case LEAST_SQUARE: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    f1.Vector<value_type>::operator=(0.0);
    unsigned int n_dof = fem_space.n_dof();
    MassMatrix<DIM, value_type> mass_matrix(fem_space);
    mass_matrix.algebricAccuracy() = algebric_accuracy;
    mass_matrix.build();
    Vector<double> rhs(n_dof);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = f0.value(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  rhs(element_dof[j]) += Jxw*f0_value*basis_value[j][l];
	}
      }
    }
    AMGSolver solver(mass_matrix);
    solver.solve(f1, rhs);
    break;
  }
  case LOCAL_LEAST_SQUARE: {
    FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
    unsigned int n_dof = fem_space.n_dof();
    std::vector<int> counter(n_dof, 0);
    f1.Vector<value_type>::operator=(0.0);
    typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
    typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
    for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      const std::vector<int>& element_dof = the_element->dof();
      unsigned int n_element_dof = element_dof.size();
      FullMatrix<double> local_mass_matrix(n_element_dof, n_element_dof);
      Vector<double> local_rhs(n_element_dof);
      Vector<double> local_f1(n_element_dof);
      unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
      for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	value_type f0_value = f0.value(q_point[l]);
	double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	for (unsigned int j = 0;j < n_element_dof;j ++) {
	  for (unsigned int k = 0;k < n_element_dof;k ++) {
	    local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	  }
	  local_rhs(j) += Jxw*f0_value*basis_value[j][l];
	}
      }
      local_mass_matrix.gauss_jordan();
      local_mass_matrix.vmult(local_f1, local_rhs);
      for (unsigned int i = 0;i < n_element_dof;i ++) {
	f1(element_dof[i]) += local_f1(i);
	counter[element_dof[i]] ++;
      }
    }
    for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
    break;
  }
  default: {
    Assert(false, ExcNotImplemented());
  }
  }
}



template <class value_type, int DIM>
  void Operator::L2Project(value_type (*f)(const value_type&), 
                           const FEMFunction<value_type,DIM>& f0,
			   FEMFunction<value_type,DIM>& f1, 
                           Method method, 
                           int algebric_accuracy)
{
  if (&f0.femSpace() == &f1.femSpace()) {
    switch (method) {
    case MASS_ACCUMULATION: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      Vector<double> mass_accumulation(n_dof);
      f1.Vector<value_type>::operator=(0.0);
      typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	std::vector<value_type> f0_value = f0.value(q_point, *the_element);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  value_type f_value = (*f)(f0_value[l]);
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    f1(element_dof[j]) += Jxw*f_value*basis_value[j][l];
	    mass_accumulation(element_dof[j]) += Jxw*basis_value[j][l];
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++)
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      f1.Vector<value_type>::operator=(0.0);
      MassMatrix<DIM, value_type> mass_matrix(fem_space);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      Vector<double> rhs(n_dof);
      typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	std::vector<value_type> f0_value = f0.value(q_point, *the_element);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  value_type f_value = (*f)(f0_value[l]);
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    rhs(element_dof[j]) += Jxw*f_value*basis_value[j][l];
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	FullMatrix<double> local_mass_matrix(n_element_dof, n_element_dof);
	Vector<double> local_rhs(n_element_dof);
	Vector<double> local_f1(n_element_dof);
	unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	std::vector<value_type> f0_value = f0.value(q_point, *the_element);
	for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	  value_type f_value = (*f)(f0_value[l]);
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    for (unsigned int k = 0;k < n_element_dof;k ++) {
	      local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	    }
	    local_rhs(j) += Jxw*f_value*basis_value[j][l];
	  }
	}
	local_mass_matrix.gauss_jordan();
	local_mass_matrix.vmult(local_f1, local_rhs);
	for (unsigned int i = 0;i < n_element_dof;i ++) {
	  f1(element_dof[i]) += local_f1(i);
	  counter[element_dof[i]] ++;
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
    }
    }
  }
  else {
    switch (method) {
    case MASS_ACCUMULATION: {
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
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> mass_accumulation(n_dof);
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
	  std::vector<value_type> f0_value = f0.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f0_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f0_value = f0.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f0_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++)
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
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
      MassMatrix<DIM, value_type> mass_matrix(fem_space1);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> rhs(n_dof);
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
	  std::vector<value_type> f0_value = f0.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f0_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f0_value = f0.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f0_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
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
      unsigned int n_dof = fem_space1.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
      ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
      ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
      int the_element_index;
      unsigned int n_the_element_dof;
      { // only used to make the variables local
	const HElement<DIM>& the_h_element = the_pair(1);
	the_element_index = the_h_element.index;
	const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	Assert(the_element.index() == the_element_index, ExcInternalError());
	n_the_element_dof = the_element.n_dof();
      }
      Vector<double> local_rhs(n_the_element_dof);
      for (;the_pair != end_pair;) {
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
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f0_value = f0.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f0_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f0_value = f0.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f0_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	++ the_pair;
	if (the_pair == end_pair || the_pair(1).index != the_element_index) {
	  const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	  const std::vector<int>& the_element_dof = the_element.dof();
	  FullMatrix<double> local_mass_matrix(n_the_element_dof, n_the_element_dof);
	  Vector<value_type> local_f1(n_the_element_dof);
	  double volume = element1.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = the_element.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = the_element.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = the_element.local_to_global(quad_info.quadraturePoint());
	  std::vector<std::vector<value_type> > basis_value = the_element.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    for (unsigned int j = 0;j < n_the_element_dof;j ++) {
	      for (unsigned int k = 0;k < n_the_element_dof;k ++) {
		local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	      }
	    }
	  }
	  local_mass_matrix.gauss_jordan();
	  local_mass_matrix.vmult(local_f1, local_rhs);
	  for (unsigned int i = 0;i < n_the_element_dof;i ++) {
	    f1(the_element_dof[i]) += local_f1(i);
	    counter[the_element_dof[i]] ++;
	  }
	  if (the_pair != end_pair) {
	    the_element_index = the_pair(1).index;
	    const Element<value_type,DIM>& the_new_element = fem_space1.element(the_element_index);
	    n_the_element_dof = the_new_element.n_dof();
	    local_rhs.reinit(n_the_element_dof);
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
    }
    }
  }
}



template <class value_type, int DIM>
  void Operator::L2Project(value_type (*f)(const value_type&,const value_type&), 
                           const FEMFunction<value_type,DIM>& f00, 
                           const FEMFunction<value_type,DIM>& f01,
			   FEMFunction<value_type,DIM>& f1, 
                           Method method, 
                           int algebric_accuracy)
{
  if ((&f00.femSpace() != &f01.femSpace())
      && (&f00.femSpace() != &f1.femSpace())
      && (&f01.femSpace() != &f1.femSpace())) {
    std::cerr << "The three FEM functions are on three different finite element spaces."
	      << std::endl;
    abort();
  }
  else if ((&f00.femSpace() == &f01.femSpace()) &&
	   (&f00.femSpace() == &f1.femSpace())) {
    switch (method) {
    case MASS_ACCUMULATION: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      Vector<double> mass_accumulation(n_dof);
      f1.Vector<value_type>::operator=(0.0);
      typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	std::vector<value_type> f00_value = f00.value(q_point, *the_element);
	std::vector<value_type> f01_value = f01.value(q_point, *the_element);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  value_type f_value = (*f)(f00_value[l], f01_value[l]);
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    f1(element_dof[j]) += Jxw*f_value*basis_value[j][l];
	    mass_accumulation(element_dof[j]) += Jxw*basis_value[j][l];
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++)
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      f1.Vector<value_type>::operator=(0.0);
      MassMatrix<DIM, value_type> mass_matrix(fem_space);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      Vector<double> rhs(n_dof);
      typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	std::vector<value_type> f00_value = f00.value(q_point, *the_element);
	std::vector<value_type> f01_value = f01.value(q_point, *the_element);
	for (int l = 0;l < n_quadrature_point;l ++) {
	  value_type f_value = (*f)(f00_value[l], f01_value[l]);
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    rhs(element_dof[j]) += Jxw*f_value*basis_value[j][l];
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
      FEMSpace<value_type,DIM>& fem_space = f1.femSpace();
      unsigned int n_dof = fem_space.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
      typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
      for (;the_element != end_element;the_element ++) {
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	const std::vector<int>& element_dof = the_element->dof();
	unsigned int n_element_dof = element_dof.size();
	FullMatrix<double> local_mass_matrix(n_element_dof, n_element_dof);
	Vector<double> local_rhs(n_element_dof);
	Vector<double> local_f1(n_element_dof);
	unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<value_type> > basis_value = the_element->basis_function_value(q_point);
	std::vector<value_type> f00_value = f00.value(q_point, *the_element);
	std::vector<value_type> f01_value = f01.value(q_point, *the_element);
	for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	  value_type f_value = (*f)(f00_value[l], f01_value[l]);
	  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	  for (unsigned int j = 0;j < n_element_dof;j ++) {
	    for (unsigned int k = 0;k < n_element_dof;k ++) {
	      local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	    }
	    local_rhs(j) += Jxw*f_value*basis_value[j][l];
	  }
	}
	local_mass_matrix.gauss_jordan();
	local_mass_matrix.vmult(local_f1, local_rhs);
	for (unsigned int i = 0;i < n_element_dof;i ++) {
	  f1(element_dof[i]) += local_f1(i);
	  counter[element_dof[i]] ++;
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
    }
    }
  }
  else if ((&f00.femSpace() == &f01.femSpace())) {
    switch (method) {
    case MASS_ACCUMULATION: {
      const FEMSpace<value_type,DIM>& fem_space0 = f00.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> mass_accumulation(n_dof);
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
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++)
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
      const FEMSpace<value_type,DIM>& fem_space0 = f00.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      f1.Vector<value_type>::operator=(0.0);
      MassMatrix<DIM, value_type> mass_matrix(fem_space1);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> rhs(n_dof);
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
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
      const FEMSpace<value_type,DIM>& fem_space0 = f00.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      unsigned int n_dof = fem_space1.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
      ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
      ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
      int the_element_index;
      unsigned int n_the_element_dof;
      { // only used to make the variables local
	const HElement<DIM>& the_h_element = the_pair(1);
	the_element_index = the_h_element.index;
	const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	Assert(the_element.index() == the_element_index, ExcInternalError());
	n_the_element_dof = the_element.n_dof();
      }
      Vector<double> local_rhs(n_the_element_dof);
      for (;the_pair != end_pair;) {
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
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	++ the_pair;
	if (the_pair == end_pair || the_pair(1).index != the_element_index) {
	  const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	  const std::vector<int>& the_element_dof = the_element.dof();
	  FullMatrix<double> local_mass_matrix(n_the_element_dof, n_the_element_dof);
	  Vector<value_type> local_f1(n_the_element_dof);
	  double volume = element1.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = the_element.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = the_element.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = the_element.local_to_global(quad_info.quadraturePoint());
	  std::vector<std::vector<value_type> > basis_value = the_element.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    for (unsigned int j = 0;j < n_the_element_dof;j ++) {
	      for (unsigned int k = 0;k < n_the_element_dof;k ++) {
		local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	      }
	    }
	  }
	  local_mass_matrix.gauss_jordan();
	  local_mass_matrix.vmult(local_f1, local_rhs);
	  for (unsigned int i = 0;i < n_the_element_dof;i ++) {
	    f1(the_element_dof[i]) += local_f1(i);
	    counter[the_element_dof[i]] ++;
	  }
	  if (the_pair != end_pair) {
	    the_element_index = the_pair(1).index;
	    Element<value_type,DIM>& the_new_element = fem_space1.element(the_element_index);
	    n_the_element_dof = the_new_element.n_dof();
	    local_rhs.reinit(n_the_element_dof);
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
    }
    }
  }
  else if (&f00.femSpace() == &f1.femSpace()) {
    switch (method) {
    case MASS_ACCUMULATION: {
      const FEMSpace<value_type,DIM>& fem_space0 = f01.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> mass_accumulation(n_dof);
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
	  std::vector<value_type> f00_value = f00.value(q_point, element1);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element1);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++)
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
      const FEMSpace<value_type,DIM>& fem_space0 = f01.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      f1.Vector<value_type>::operator=(0.0);
      MassMatrix<DIM, value_type> mass_matrix(fem_space1);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> rhs(n_dof);
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
	  std::vector<value_type> f00_value = f00.value(q_point, element1);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element1);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
      const FEMSpace<value_type,DIM>& fem_space0 = f01.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      unsigned int n_dof = fem_space1.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
      ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
      ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
      int the_element_index;
      unsigned int n_the_element_dof;
      { // only used to make the variables local
	const HElement<DIM>& the_h_element = the_pair(1);
	the_element_index = the_h_element.index;
	const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	Assert(the_element.index() == the_element_index, ExcInternalError());
	n_the_element_dof = the_element.n_dof();
      }
      Vector<double> local_rhs(n_the_element_dof);
      for (;the_pair != end_pair;) {
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
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element1);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element1);
	  std::vector<value_type> f01_value = f01.value(q_point, element0);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	++ the_pair;
	if (the_pair == end_pair || the_pair(1).index != the_element_index) {
	  const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	  const std::vector<int>& the_element_dof = the_element.dof();
	  FullMatrix<double> local_mass_matrix(n_the_element_dof, n_the_element_dof);
	  Vector<value_type> local_f1(n_the_element_dof);
	  double volume = element1.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = the_element.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = the_element.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = the_element.local_to_global(quad_info.quadraturePoint());
	  std::vector<std::vector<value_type> > basis_value = the_element.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    for (unsigned int j = 0;j < n_the_element_dof;j ++) {
	      for (unsigned int k = 0;k < n_the_element_dof;k ++) {
		local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	      }
	    }
	  }
	  local_mass_matrix.gauss_jordan();
	  local_mass_matrix.vmult(local_f1, local_rhs);
	  for (unsigned int i = 0;i < n_the_element_dof;i ++) {
	    f1(the_element_dof[i]) += local_f1(i);
	    counter[the_element_dof[i]] ++;
	  }
	  if (the_pair != end_pair) {
	    the_element_index = the_pair(1).index;
	    const Element<value_type,DIM>& the_new_element = fem_space1.element(the_element_index);
	    n_the_element_dof = the_new_element.n_dof();
	    local_rhs.reinit(n_the_element_dof);
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
    }
    }
  }
  else if (&f01.femSpace() == &f1.femSpace()) {
    switch (method) {
    case MASS_ACCUMULATION: {
      const FEMSpace<value_type,DIM>& fem_space0 = f00.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> mass_accumulation(n_dof);
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
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element1);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element1);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      f1(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	      mass_accumulation(element_dof1[j]) += Jxw*basis_value1[j][l];
	    }
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++)
	f1(i) /= mass_accumulation(i);
      break;
    }
    case LEAST_SQUARE: {
      const FEMSpace<value_type,DIM>& fem_space0 = f00.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      f1.Vector<value_type>::operator=(0.0);
      MassMatrix<DIM, value_type> mass_matrix(fem_space1);
      mass_matrix.algebricAccuracy() = algebric_accuracy;
      mass_matrix.build();
      unsigned int n_dof = fem_space1.n_dof();
      Vector<double> rhs(n_dof);
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
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element1);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element1);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      rhs(element_dof1[j]) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
      }
      AMGSolver solver(mass_matrix);
      solver.solve(f1, rhs);
      break;
    }
    case LOCAL_LEAST_SQUARE: {
      const FEMSpace<value_type,DIM>& fem_space0 = f00.femSpace();
      FEMSpace<value_type,DIM>& fem_space1 = f1.femSpace();
      const RegularMesh<DIM>& regular_mesh0 = static_cast<const RegularMesh<DIM> &>(fem_space0.mesh());
      RegularMesh<DIM>& regular_mesh1 = static_cast<RegularMesh<DIM> &>(fem_space1.mesh());
      IrregularMesh<DIM>& irregular_mesh0 = regular_mesh0.irregularMesh();
      IrregularMesh<DIM>& irregular_mesh1 = regular_mesh1.irregularMesh();
      if (&(irregular_mesh0.geometryTree()) != &(irregular_mesh1.geometryTree())) {
	std::cerr << "The two FEM functions are even not on the same hierarchy geometry tree."
		  << std::endl;
	Assert(false, ExcInternalError());
      }
      unsigned int n_dof = fem_space1.n_dof();
      std::vector<int> counter(n_dof, 0);
      f1.Vector<value_type>::operator=(0.0);
      IrregularMeshPair<DIM> mesh_pair(irregular_mesh0, irregular_mesh1);
      ActiveElementPairIterator<DIM> the_pair = mesh_pair.beginActiveElementPair();
      ActiveElementPairIterator<DIM> end_pair = mesh_pair.endActiveElementPair();
      int the_element_index;
      unsigned int n_the_element_dof;
      { // only used to make the variables local
	const HElement<DIM>& the_h_element = the_pair(1);
	the_element_index = the_h_element.index;
	const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	Assert(the_element.index() == the_element_index, ExcInternalError());
	n_the_element_dof = the_element.n_dof();
      }
      Vector<double> local_rhs(n_the_element_dof);
      for (;the_pair != end_pair;) {
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
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element1);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	else {
	  double volume = element0.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	  std::vector<value_type> f00_value = f00.value(q_point, element0);
	  std::vector<value_type> f01_value = f01.value(q_point, element1);
	  std::vector<std::vector<value_type> > basis_value1 = element1.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    value_type f_value = (*f)(f00_value[l], f01_value[l]);
	    for (unsigned int j = 0;j < n_element_dof1;j ++) {
	      local_rhs(j) += Jxw*f_value*basis_value1[j][l];
	    }
	  }
	}
	++ the_pair;
	if (the_pair == end_pair || the_pair(1).index != the_element_index) {
	  const Element<value_type,DIM>& the_element = fem_space1.element(the_element_index);
	  const std::vector<int>& the_element_dof = the_element.dof();
	  FullMatrix<double> local_mass_matrix(n_the_element_dof, n_the_element_dof);
	  Vector<value_type> local_f1(n_the_element_dof);
	  double volume = element1.templateElement().volume();
	  const QuadratureInfo<DIM>& quad_info = the_element.findQuadratureInfo(algebric_accuracy);
	  std::vector<double> jacobian = the_element.local_to_global_jacobian(quad_info.quadraturePoint());
	  unsigned int n_quadrature_point = quad_info.n_quadraturePoint();
	  std::vector<Point<DIM> > q_point = the_element.local_to_global(quad_info.quadraturePoint());
	  std::vector<std::vector<value_type> > basis_value = the_element.basis_function_value(q_point);
	  for (unsigned int l = 0;l < n_quadrature_point;l ++) {
	    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
	    for (unsigned int j = 0;j < n_the_element_dof;j ++) {
	      for (unsigned int k = 0;k < n_the_element_dof;k ++) {
		local_mass_matrix(j, k) += Jxw*basis_value[j][l]*basis_value[k][l];
	      }
	    }
	  }
	  local_mass_matrix.gauss_jordan();
	  local_mass_matrix.vmult(local_f1, local_rhs);
	  for (unsigned int i = 0;i < n_the_element_dof;i ++) {
	    f1(the_element_dof[i]) += local_f1(i);
	    counter[the_element_dof[i]] ++;
	  }
	  if (the_pair != end_pair) {
	    the_element_index = the_pair(1).index;
	    const Element<value_type,DIM>& the_new_element = fem_space1.element(the_element_index);
	    n_the_element_dof = the_new_element.n_dof();
	    local_rhs.reinit(n_the_element_dof);
	  }
	}
      }
      for (unsigned int i = 0;i < n_dof;i ++) f1(i) /= counter[i];
      break;
    }
    default: {
      Assert(false, ExcNotImplemented());
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
