/**
 * @file   Functional.templates.h
 * @author Robert Lie
 * @date   Thu Oct 28 08:53:51 2004
 * 
 * @brief  
 * 
 * 
 */

#ifndef _Functional_templates_h_
#define _Functional_templates_h_

#include "Functional.h"

AFEPACK_OPEN_NAMESPACE

template <class value_type, int DIM>
value_type Functional::L1Norm(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	value_type norm = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			norm += Jxw*fabs(f_value[l]);
		}
	}
	return norm;
};



template <class value_type, int DIM>
value_type Functional::L2Norm(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	value_type norm = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			norm += Jxw*f_value[l]*f_value[l];
		}
	}
	return sqrt(fabs(norm));
};



template <class value_type, int DIM>
value_type Functional::L0Norm(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	value_type norm = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			f_value[l] = fabs(f_value[l]);
			if (f_value[l] > norm) norm = f_value[l];
		}
	}
	return norm;
};



template <class value_type, int DIM>
value_type Functional::LpNorm(FEMFunction<value_type, DIM>& f, double p, int algebric_accuracy)
{
	value_type norm = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			norm += Jxw*pow(f_value[l], p);
		}
	}
	return pow(norm, 1./p);
};



template <class value_type, int DIM>
value_type Functional::W11Seminorm(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	int i;
	value_type norm, n[DIM];
	for (i = 0;i < DIM;i ++) n[i] = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			for (i = 0;i < DIM;i ++)
				n[i] += Jxw*fabs(f_gradient[l][i]);
		}
	}
	for (norm = 0, i = 0;i < DIM;i ++) norm += n[i];
	return norm;
};



template <class value_type, int DIM>
value_type Functional::H1Seminorm(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	value_type norm = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			for (int i = 0;i < DIM;i ++)
				norm += Jxw*f_gradient[l][i]*f_gradient[l][i];
		}
	}
	return sqrt(fabs(norm));
};



template <class value_type, int DIM>
value_type Functional::W1pSeminorm(FEMFunction<value_type, DIM>& f, double p, int algebric_accuracy)
{
	value_type norm = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			for (int i = 0;i < DIM;i ++)
				norm += Jxw*pow(f_gradient[l][i], p);
		}
	}
	return pow(norm, 1./p);
};



template <class value_type, int DIM>
value_type Functional::W10Seminorm(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	int i;
	value_type norm, n[DIM];
	for (i = 0;i < DIM;i ++) n[i] = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			for (i = 0;i < DIM;i ++) {
				f_gradient[l][i] = fabs(f_gradient[l][i]);
				if (n[i] < f_gradient[l][i])
					n[i] = f_gradient[l][i];
			}
		}
	}
	for (norm = 0, i = 0;i < DIM;i ++) norm += n[i];
	return norm;
};



template <class value_type, int DIM>
value_type Functional::meanValue(FEMFunction<value_type, DIM>& f, int algebric_accuracy)
{
	value_type val = 0;
	double vol = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
		  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
		  val += Jxw*f_value[l];
		  vol += Jxw;
		}
	}
	return val/vol;
};



template <class value_type, int DIM>
value_type Functional::meanValue(const Function<value_type>& f, 
				 FEMSpace<value_type,DIM>& fem_space,
				 int algebric_accuracy)
{
	value_type val = 0;
	double vol = 0;
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		for (int l = 0;l < n_quadrature_point;l ++) {
		  double Jxw = quad_info.weight(l)*jacobian[l]*volume;
		  val += Jxw*f.value(q_point[l]);
		  vol += Jxw;
		}
	}
	return val/vol;
};


template <class value_type, int DIM>
value_type Functional::L1Error(FEMFunction<value_type, DIM>& f, 
			       const Function<value_type>& f1, 
			       int algebric_accuracy)
{
	value_type error = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			double df_value = f1.value(q_point[l]) - f_value[l];
			error += Jxw*fabs(df_value);
		}
	}
	return error;
};



template <class value_type, int DIM>
value_type Functional::L2Error(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
{
	value_type error = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			double df_value = f1.value(q_point[l]) - f_value[l];
			error += Jxw*df_value*df_value;
		}
	}
	return sqrt(fabs(error));
};


template <class value_type, int DIM>
value_type Functional::L0Error(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
{
	value_type error = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double df_value = f1.value(q_point[l]) - f_value[l];
			df_value = fabs(df_value);
			if (df_value > error) error = df_value;
		}
	}
	return error;
};


template <class value_type, int DIM>
value_type Functional::LpError(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, double p, int algebric_accuracy)
{
	value_type error = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<double> f_value = f.value(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			double df_value = f1.value(q_point[l]) - f_value[l];
			error += Jxw*pow(fabs(df_value), p);
		}
	}
	return pow(error, 1.0/p);
};



template <class value_type, int DIM>
value_type Functional::W11SemiError(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
{
	int i;
	value_type error, n[DIM];
	for (i = 0;i < DIM;i ++) n[i] = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			std::vector<value_type> f1_gradient = f1.gradient(q_point[l]);
			for (i = 0;i < DIM;i ++) {
				n[i] += Jxw*fabs(f_gradient[l][i] - f1_gradient[i]);
			}
		}
	}
	for (error = 0, i = 0;i < DIM;i ++) error += n[i];
	return error;
};



template <class value_type, int DIM>
value_type Functional::H1SemiError(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
{
	int i;
	value_type error, n[DIM];
	for (i = 0;i < DIM;i ++) n[i] = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			std::vector<value_type> f1_gradient = f1.gradient(q_point[l]);
			for (i = 0;i < DIM;i ++) {
				double df_value = f_gradient[l][i] - f1_gradient[i];
				n[i] += Jxw*df_value*df_value;
			}
		}
	}
	for (error = 0, i = 0;i < DIM;i ++) error += n[i];
	return sqrt(fabs(error));
};



template <class value_type, int DIM>
value_type Functional::W1pSemiError(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, double p, int algebric_accuracy)
{
	int i;
	value_type error, n[DIM];
	for (i = 0;i < DIM;i ++) n[i] = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			std::vector<value_type> f1_gradient = f1.gradient(q_point[l]);
			for (i = 0;i < DIM;i ++) {
				double df_value = fabs(f_gradient[l][i] - f1_gradient[i]);
				n[i] += Jxw*pow(df_value, p);
			}
		}
	}
	for (error = 0, i = 0;i < DIM;i ++) error += n[i];
	return pow(error, 1./p);
};



template <class value_type, int DIM>
value_type Functional::W10SemiError(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
{
	int i;
	value_type error, n[DIM];
	for (i = 0;i < DIM;i ++) n[i] = 0;
	FEMSpace<value_type,DIM>& fem_space = f.femSpace();
	typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
	typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > f_gradient = f.gradient(q_point, *the_element);
		for (int l = 0;l < n_quadrature_point;l ++) {
			std::vector<value_type> f1_gradient = f1.gradient(q_point[l]);
			for (i = 0;i < DIM;i ++) {
				double df_value = fabs(f_gradient[l][i] - f1_gradient[i]);
				if (n[i] < df_value) n[i] = df_value;
			}
		}
	}
	for (error = 0, i = 0;i < DIM;i ++) error += n[i];
	return error;
};


template <class value_type, int DIM>
value_type Functional::L1Norm(const Function<value_type>& f, 
			      FEMSpace<value_type,DIM>& fem_space,
			      int algebric_accuracy)
{
  value_type norm = 0;
  typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
  typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      double f_value = f.value(q_point[l]);
      norm += Jxw*fabs(f_value);
    }
  }
  return norm;
};



template <class value_type, int DIM>
value_type Functional::L2Norm(const Function<value_type>& f, 
			      FEMSpace<value_type,DIM>& fem_space,
			      int algebric_accuracy)
{
  value_type norm = 0;
  typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
  typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    double volume = the_element->templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
    std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      double f_value = f.value(q_point[l]);
      norm += Jxw*f_value*f_value;
    }
  }
  return sqrt(fabs(norm));
};

AFEPACK_CLOSE_NAMESPACE

#endif // _Functional_templates_h_

/**
 * end of file
 * 
 */
