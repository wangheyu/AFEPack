/**
 * @file   TemplateElement.cpp
 * @author Robert Lie
 * @date   Wed Jun 21 08:18:41 2003
 * 
 * @brief  
 * 
 * 
 */

#include "TemplateElement.templates.h"

AFEPACK_OPEN_NAMESPACE

#define value_type double
#define DIM 1

template class BasisFunctionIdentity<DIM>;
template class ShapeFunction<value_type,DIM>;
template class TemplateDOF<DIM>;
template class UnitOutNormal<DIM>;

template bool operator==(const BasisFunctionIdentity<DIM>&, const BasisFunctionIdentity<DIM>&);

#define TDIM 1
template class BasisFunction<value_type,DIM,TDIM>;
template class CoordTransform<TDIM,DIM>;
template class BasisFunctionAdmin<value_type,DIM,TDIM>;
template class TemplateElement<value_type,DIM,TDIM>;
#undef TDIM

template struct GeometryAdditionalData<DIM>;

#undef DIM

#define DIM 2

template class BasisFunctionIdentity<DIM>;
template class ShapeFunction<value_type,DIM>;
template class TemplateDOF<DIM>;
template class UnitOutNormal<DIM>;

template bool operator==(const BasisFunctionIdentity<DIM>&, const BasisFunctionIdentity<DIM>&);
#define TDIM 1
template class BasisFunction<value_type,DIM,TDIM>;
template class CoordTransform<TDIM,DIM>;
template class BasisFunctionAdmin<value_type,DIM,TDIM>;
template class TemplateElement<value_type,DIM,TDIM>;
#undef TDIM

#define TDIM 2
template class BasisFunction<value_type,DIM,TDIM>;
template class CoordTransform<TDIM,DIM>;
template class BasisFunctionAdmin<value_type,DIM,TDIM>;
template class TemplateElement<value_type,DIM,TDIM>;
#undef TDIM

template struct GeometryAdditionalData<DIM>;

#undef DIM

#define DIM 3

template class BasisFunctionIdentity<DIM>;
template class ShapeFunction<value_type,DIM>;
template class TemplateDOF<DIM>;
template class UnitOutNormal<DIM>;

template bool operator==(const BasisFunctionIdentity<DIM>&, const BasisFunctionIdentity<DIM>&);

#define TDIM 1
template class BasisFunction<value_type,DIM,TDIM>;
template class CoordTransform<TDIM,DIM>;
template class BasisFunctionAdmin<value_type,DIM,TDIM>;
template class TemplateElement<value_type,DIM,TDIM>;

#undef TDIM
#define TDIM 2
template class BasisFunction<value_type,DIM,TDIM>;
template class CoordTransform<TDIM,DIM>;
template class BasisFunctionAdmin<value_type,DIM,TDIM>;
template class TemplateElement<value_type,DIM,TDIM>;

#undef TDIM
#define TDIM 3
template class BasisFunction<value_type,DIM,TDIM>;
template class CoordTransform<TDIM,DIM>;
template class BasisFunctionAdmin<value_type,DIM,TDIM>;
template class TemplateElement<value_type,DIM,TDIM>;
#undef TDIM

template struct GeometryAdditionalData<DIM>;

#undef DIM

#undef value_type

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

