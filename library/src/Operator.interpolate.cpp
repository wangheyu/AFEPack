/**
 * @file   Operator.interpolate.cpp
 * @author Robert Lie
 * @date   Thu Mar 30 11:46:13 2006
 * 
 * @brief  interpolations
 * 
 * 
 */

#include "Operator.interpolate.templates.h"

AFEPACK_OPEN_NAMESPACE

#define value_type double
#define DIM 2
namespace Operator {
	template void L2Interpolate(const FEMFunction<value_type, DIM>&, FEMFunction<value_type, DIM>&);
	template void L2Interpolate(value_type (*)(const double *), FEMFunction<value_type, DIM>&);
	template void L2Interpolate(value_type (*)(const Point<DIM>&), FEMFunction<value_type, DIM>&);
	template void L2Interpolate(const Function<value_type>&, FEMFunction<value_type, DIM>&);
};

#undef DIM

#define DIM 3
namespace Operator {
	template void L2Interpolate(const FEMFunction<value_type, DIM>&, FEMFunction<value_type, DIM>&);
	template void L2Interpolate(value_type (*)(const double *), FEMFunction<value_type, DIM>&);
	template void L2Interpolate(value_type (*)(const Point<DIM>&), FEMFunction<value_type, DIM>&);
	template void L2Interpolate(const Function<value_type>&, FEMFunction<value_type, DIM>&);
};

#undef DIM
#undef value_type

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

