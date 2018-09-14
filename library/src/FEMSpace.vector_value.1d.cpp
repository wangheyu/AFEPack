/**
 * @file   FEMSpace.vector_value.1d.cpp
 * @author Robert Lie
 * @date   Thu Mar 30 11:30:58 2006
 * 
 * @brief  vecter valued finite element space for DOW=1
 * 
 * 
 */

#include "FEMSpace.templates.h"

#define DOW 1
#define DIM 1
#define TDIM 1

#define vector_length 1
#include "FEMSpace.vector_value.templates.h"
#undef _FEMSpace_vector_value_templates_h_

AFEPACK_OPEN_NAMESPACE

#define value_type nVector<vector_length,double>
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

#define vector_length 2
#include "FEMSpace.vector_value.templates.h"
#undef _FEMSpace_vector_value_templates_h_

AFEPACK_OPEN_NAMESPACE

#define value_type nVector<vector_length,double>
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

#define vector_length 3
#include "FEMSpace.vector_value.templates.h"
#undef _FEMSpace_vector_value_templates_h_

AFEPACK_OPEN_NAMESPACE

#define value_type nVector<vector_length,double>
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

#undef TDIM
#undef DIM
#undef DOW

/**
 * end of file
 * 
 */
