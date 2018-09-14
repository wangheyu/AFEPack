//-*-Mode: C++;-*-
#ifndef _HybridFEMSpace_h_
#define _HybridFEMSpace_h_

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "baseclass.h"

//---- HybridFEMSpace -----------------------------------------------------------

template <DIM, vt0, vt1>
class HybridFEMSpace
{
private:
	FEMSpace<vt0, DIM, DIM, DIM> * sp0;
	FEMSpace<vt1, DIM, DOW, DIM-1> * sp1;
public:
    HybridFEMSpace();
    ~HybridFEMSpace();

};

#endif
