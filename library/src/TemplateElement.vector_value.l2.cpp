/**
 * @file   TemplateElement.vector_value.l2.cpp
 * @author Robert Lie
 * @date   Thu Mar 30 11:50:00 2006
 * 
 * @brief  template element with vector length=2
 * 
 * 
 */

#include "TemplateElement.templates.h"

AFEPACK_OPEN_NAMESPACE

#define vector_length 2
#define value_type nVector<vector_length,double>
#define DIM 1
	template class ShapeFunction<value_type,DIM>;

#define TDIM 1
	template class BasisFunction<value_type,DIM,TDIM>;
	template class BasisFunctionAdmin<value_type,DIM,TDIM>;
	template class TemplateElement<value_type,DIM,TDIM>;


#undef TDIM
#undef DIM

#define DIM 2
	template class ShapeFunction<value_type,DIM>;

#define TDIM 1
	template class BasisFunction<value_type,DIM,TDIM>;
	template class BasisFunctionAdmin<value_type,DIM,TDIM>;
	template class TemplateElement<value_type,DIM,TDIM>;


#undef TDIM
#define TDIM 2
	template class BasisFunction<value_type,DIM,TDIM>;
	template class BasisFunctionAdmin<value_type,DIM,TDIM>;
	template class TemplateElement<value_type,DIM,TDIM>;


#undef TDIM
#undef DIM

#define DIM 3
	template class ShapeFunction<value_type,DIM>;

#define TDIM 1
	template class BasisFunction<value_type,DIM,TDIM>;
	template class BasisFunctionAdmin<value_type,DIM,TDIM>;
	template class TemplateElement<value_type,DIM,TDIM>;


#undef TDIM
#define TDIM 2
	template class BasisFunction<value_type,DIM,TDIM>;
	template class BasisFunctionAdmin<value_type,DIM,TDIM>;
	template class TemplateElement<value_type,DIM,TDIM>;


#undef TDIM
#define TDIM 3
	template class BasisFunction<value_type,DIM,TDIM>;
	template class BasisFunctionAdmin<value_type,DIM,TDIM>;
	template class TemplateElement<value_type,DIM,TDIM>;


#undef TDIM
#undef DIM
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
