/**
 * @file   Functional.h
 * @author Robert Lie
 * @date   Thu Oct 28 08:51:59 2004
 * 
 * @brief  
 * 
 * 
 */

#ifndef _Functional_h_
#define _Functional_h_

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
 * This namespace packs a list of functionals on a finite element function or several
 * finite element functions. Those functionals are generally related with a numerical
 * intergration on the whole domain.
 */
namespace Functional {
/**
 * \f$L^1\f$ norm of a finite element function.
 */
	template <class value_type, int DIM> value_type L1Norm(FEMFunction<value_type, DIM>&, int);
/**
 * \f$L^2\f$ norm of a finite element function.
 */
	template <class value_type, int DIM> value_type L2Norm(FEMFunction<value_type, DIM>&, int);
/**
 * \f$L^\infty\f$ norm of a finite element function.
 */
	template <class value_type, int DIM> value_type L0Norm(FEMFunction<value_type, DIM>&, int);
/**
 * \f$L^p\f$ norm of a finite element function.
 */
	template <class value_type, int DIM> value_type LpNorm(FEMFunction<value_type, DIM>&, double, int);
/**
 * \f$W^{1,1}\f$ semi-norm of a finite element function.
 */
	template <class value_type, int DIM> value_type W11Seminorm(FEMFunction<value_type, DIM>&, int);
/**
 * \f$H^1\f$ semi-norm of a finite element function.
 */
	template <class value_type, int DIM> value_type H1Seminorm(FEMFunction<value_type, DIM>&, int);
/**
 * \f$W^{1,p}\f$ semi-norm of a finite element function.
 */
	template <class value_type, int DIM> value_type W1pSeminorm(FEMFunction<value_type, DIM>&, double, int);
/**
 * \f$W^{1,\infty}\f$ semi-norm of a finite element function.
 */
	template <class value_type, int DIM> value_type W10Seminorm(FEMFunction<value_type, DIM>&, int);
	

	/**
	 * Calculate the mean value of the finite elment function \f$ \int_\Omega u_h 
	 * dx / | \Omega | \f$.
	 * 
	 */
	template <class value_type, int DIM> 
	  value_type meanValue(FEMFunction<value_type, DIM>&, 
			       int);
	/**
	 * Calculate the mean value of the analytic funcion by numberical
	 * integration in the given finite element space.
	 * 
	 */
	template <class value_type, int DIM> 
	  value_type meanValue(const Function<value_type>&, 
			       FEMSpace<value_type,DIM>&,
			       int);

	template <class value_type, int DIM>
	value_type L1Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template <class value_type, int DIM>
	value_type L2Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template <class value_type, int DIM>
	value_type L0Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template <class value_type, int DIM>
	value_type LpError(FEMFunction<value_type, DIM>&, const Function<value_type>&, double, int);
	
	template <class value_type, int DIM>
	value_type W11SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template <class value_type, int DIM>
	value_type H1SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template <class value_type, int DIM>
	value_type W1pSemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, double, int);
	template <class value_type, int DIM>
	value_type W10SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);

	/**
	 * \f$L^1\f$ norm of an analytic function.
	 */
	template <class value_type, int DIM> 
	  value_type L1Norm(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);
	/**
	 * \f$L^2\f$ norm of an analytic function.
	 */
	template <class value_type, int DIM> 
	  value_type L2Norm(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);

};

AFEPACK_CLOSE_NAMESPACE

#endif // _Functional_h_

/**
 * end of file
 * 
 */
