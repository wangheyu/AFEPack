/**
 * @file   Operator.project.cpp
 * @author Robert Lie
 * @date   Thu Mar 30 11:46:42 2006
 * 
 * @brief  projections
 * 
 * 
 */

#include "Operator.project.templates.h"

AFEPACK_OPEN_NAMESPACE

#define value_type double
#define DIM 2
namespace Operator {
	template void L2Project(const FEMFunction<value_type, DIM>&, 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(value_type (*)(const double *), 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(value_type (*)(const Point<DIM>&), 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(const Function<value_type>&, 
                                FEMFunction<value_type, DIM>&, Method, int);

	template void L2Project(value_type (*)(const value_type&), 
                                const FEMFunction<value_type, DIM>&, 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(value_type (*)(const value_type&, const value_type&), 
                                const FEMFunction<value_type, DIM>&, 
                                const FEMFunction<value_type, DIM>&, 
                                FEMFunction<value_type, DIM>&, Method, int);
};

#undef DIM

#define DIM 3
namespace Operator {
	template void L2Project(const FEMFunction<value_type, DIM>&, 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(value_type (*)(const double *), 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(value_type (*)(const Point<DIM>&), 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(const Function<value_type>&, 
                                FEMFunction<value_type, DIM>&, Method, int);

	template void L2Project(value_type (*)(const value_type&), 
                                const FEMFunction<value_type, DIM>&, 
                                FEMFunction<value_type, DIM>&, Method, int);
	template void L2Project(value_type (*)(const value_type&, const value_type&), 
                                const FEMFunction<value_type, DIM>&, 
                                const FEMFunction<value_type, DIM>&, 
                                FEMFunction<value_type, DIM>&, Method, int);
};

#undef DIM
#undef value_type

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

