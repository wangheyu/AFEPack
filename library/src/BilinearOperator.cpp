///////////////////////////////////////////////////////////////////////////
// BilinearOperator.cpp : by R.Lie
//

#include "BilinearOperator.templates.h"

AFEPACK_OPEN_NAMESPACE

#define DOW         DIM
#define TDIM        DIM
#define TDIM0       DIM
#define TDIM1       DIM

#define value_type0 double
#define value_type1 double
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;
template class L2InnerProduct<DIM,value_type0, value_type1,DOW,TDIM0,TDIM1>;
template class MassMatrix<DIM,value_type0,DOW,TDIM>;
template class StiffMatrix<DIM,value_type0,DOW,TDIM>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;
template class L2InnerProduct<DIM,value_type0, value_type1,DOW,TDIM0,TDIM1>;
template class MassMatrix<DIM,value_type0,DOW,TDIM>;
template class StiffMatrix<DIM,value_type0,DOW,TDIM>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;
template class L2InnerProduct<DIM,value_type0, value_type1,DOW,TDIM0,TDIM1>;
template class MassMatrix<DIM,value_type0,DOW,TDIM>;
template class StiffMatrix<DIM,value_type0,DOW,TDIM>;

#undef DIM
#undef value_type1
#undef value_type0

///////////////////////////////////////////////////////////////////////////

#define vector_length 2
#define value_type0 nVector<vector_length,double>
#define value_type1 double
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM
#undef value_type1
#undef value_type0

#define value_type0 double
#define value_type1 nVector<vector_length,double>
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM
#undef value_type1
#undef value_type0

#define value_type0 nVector<vector_length,double>
#define value_type1 nVector<vector_length,double>
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM
#undef value_type1
#undef value_type0
#undef vector_length

///////////////////////////////////////////////////////////////////////////

#define vector_length 3
#define value_type0 nVector<vector_length,double>
#define value_type1 double
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM
#undef value_type1
#undef value_type0

#define value_type0 double
#define value_type1 nVector<vector_length,double>
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM
#undef value_type1
#undef value_type0

#define value_type0 nVector<vector_length,double>
#define value_type1 nVector<vector_length,double>
#define DIM 1

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 2

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM

#define DIM 3

template class BilinearOperator<DIM,value_type0,value_type1,DOW,TDIM0,TDIM1>;

#undef DIM
#undef value_type1
#undef value_type0
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

//
// end of file
///////////////////////////////////////////////////////////////////////////
