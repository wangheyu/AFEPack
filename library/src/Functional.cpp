/**
 * @file   Functional.cpp
 * @author Robert Lie
 * @date   Thu Oct 28 08:55:18 2004
 * 
 * @brief  
 * 
 * 
 */

#include "Functional.templates.h"

AFEPACK_OPEN_NAMESPACE

#define value_type double
#define DIM 2

namespace Functional {
	template value_type L1Norm(FEMFunction<value_type, DIM>&, int);
	template value_type L2Norm(FEMFunction<value_type, DIM>&, int);
	template value_type L0Norm(FEMFunction<value_type, DIM>&, int);
	template value_type LpNorm(FEMFunction<value_type, DIM>&, double, int);

	template value_type W11Seminorm(FEMFunction<value_type, DIM>&, int);
	template value_type H1Seminorm(FEMFunction<value_type, DIM>&, int);
	template value_type W1pSeminorm(FEMFunction<value_type, DIM>&, double, int);
	template value_type W10Seminorm(FEMFunction<value_type, DIM>&, int);

	template value_type L1Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type L2Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type L0Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type LpError(FEMFunction<value_type, DIM>&, const Function<value_type>&, double, int);
	
	template value_type W11SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type H1SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type W1pSemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, double, int);
	template value_type W10SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);

  template value_type meanValue(FEMFunction<value_type,DIM>&, int);
  template value_type meanValue(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);

  template value_type L1Norm(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);
  template value_type L2Norm(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);
};

#undef DIM

#define DIM 3

namespace Functional {
	template value_type L1Norm(FEMFunction<value_type, DIM>&, int);
	template value_type L2Norm(FEMFunction<value_type, DIM>&, int);
	template value_type L0Norm(FEMFunction<value_type, DIM>&, int);
	template value_type LpNorm(FEMFunction<value_type, DIM>&, double, int);

	template value_type W11Seminorm(FEMFunction<value_type, DIM>&, int);
	template value_type H1Seminorm(FEMFunction<value_type, DIM>&, int);
	template value_type W1pSeminorm(FEMFunction<value_type, DIM>&, double, int);
	template value_type W10Seminorm(FEMFunction<value_type, DIM>&, int);
	
	template value_type L1Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type L2Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type L0Error(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type LpError(FEMFunction<value_type, DIM>&, const Function<value_type>&, double, int);
	
	template value_type W11SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type H1SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);
	template value_type W1pSemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, double, int);
	template value_type W10SemiError(FEMFunction<value_type, DIM>&, const Function<value_type>&, int);

  template value_type meanValue(FEMFunction<value_type,DIM>&, int);
  template value_type meanValue(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);

  template value_type L1Norm(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);
  template value_type L2Norm(const Function<value_type>&, FEMSpace<value_type,DIM>&, int);
};

#undef DIM
#undef value_type

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
