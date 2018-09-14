/**
 * @file   HGeometry.1d.cpp
 * @author Robert Lie
 * @date   Sun Apr 29 10:48:51 2007
 * 
 * @brief  
 * 
 * 
 */

#include "HGeometry.templates.h"

AFEPACK_OPEN_NAMESPACE

#define DOW 1
#include "HGeometry.templates.nd.h"

template <> HGeometry<0,DOW> * HGeometry<0,DOW>::parent = NULL;
template <> std::vector<HGeometry<0,DOW> *> HGeometry<0,DOW>::vertex(1, (HGeometry<0,DOW> *)NULL);
template <> std::vector<HGeometry<0,DOW> *> HGeometry<0,DOW>::boundary(1, (HGeometry<0,DOW> *)NULL);
template <> std::vector<HGeometry<0,DOW> *> HGeometry<0,DOW>::child(1, (HGeometry<0,DOW> *)NULL);

template class HGeometry<0,DOW>;
template class HGeometry<1,DOW>;
template class HGeometry<2,DOW>;
template class HGeometry<3,DOW>;

template <>
void (*HGeometry<1,DOW>::mid_point)(const Point<DOW>&,
				    const Point<DOW>&,
				    bmark_t, 
				    Point<DOW>&) = NULL;

template std::ostream& operator<<(std::ostream&, const HGeometry<0,DOW>&);
template std::ostream& operator<<(std::ostream&, const HGeometry<1,DOW>&);
template std::ostream& operator<<(std::ostream&, const HGeometry<2,DOW>&);
template std::ostream& operator<<(std::ostream&, const HGeometry<3,DOW>&);

#define DIM 1
#include "HGeometry.cpptemplate"
#undef DIM

#define DIM 2
#include "HGeometry.cpptemplate"
#undef DIM

#define DIM 3
#include "HGeometry.cpptemplate"
#undef DIM

#undef DOW

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
