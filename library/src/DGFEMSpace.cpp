///////////////////////////////////////////////////////////////////////////
// DGFEMSpace.cpp : by R.Lie
//

#include "DGFEMSpace.templates.h"

AFEPACK_OPEN_NAMESPACE

#define value_type double
#define DIM 2
#define DOW 2
#define TDIM 2
#define TDIM1 1

template class TemplateDGElement<TDIM1,DOW>;
template class DGElement<value_type,DIM,DOW,TDIM,TDIM1>;
template class DGFEMSpace<value_type,DIM,DOW,TDIM,TDIM1>;

template std::vector<double> unitOutNormal(const Point<DIM>&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DIM> >&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);

#undef TDIM1
#undef TDIM
#undef DOW

template struct DGElementAdditionalData<value_type,DIM>;

#undef DIM

#define DIM 3
#define DOW 3
#define TDIM 3
#define TDIM1 2

template class TemplateDGElement<TDIM1,DOW>;
template class DGElement<value_type,DIM,DOW,TDIM,TDIM1>;
template class DGFEMSpace<value_type,DIM,DOW,TDIM,TDIM1>;

template std::vector<double> unitOutNormal(const Point<DIM>&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DIM> >&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);

#undef TDIM1
#undef TDIM
#undef DOW

template struct DGElementAdditionalData<value_type,DIM>;

#undef DIM
#undef value_type

AFEPACK_CLOSE_NAMESPACE

//
// end of file
///////////////////////////////////////////////////////////////////////////
