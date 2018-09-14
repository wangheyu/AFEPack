/**
 * @file   Geometry.3d.cpp
 * @author Robert Lie
 * @date   Sun Apr 29 21:14:40 2007
 * 
 * @brief  
 * 
 * 
 */

#include "Geometry.templates.h"

AFEPACK_OPEN_NAMESPACE

#define DOW 3

	template class Point<DOW>;
	template class QuadratureInfo<DOW>;
	template class QuadratureInfoAdmin<DOW>;
	template class TemplateGeometry<DOW>;

	template Point<DOW> midpoint(const Point<DOW>&, const Point<DOW>&);
	template double distance(const Point<DOW>&, const Point<DOW>&);
	template Point<DOW> barycenter(const std::vector<Point<DOW> >&, const double *);
	template Point<DOW> operator+(const Point<DOW>&, const Point<DOW>&);
	template Point<DOW> operator-(const Point<DOW>&, const Point<DOW>&);

	template std::istream& operator>>(std::istream&, Point<DOW>&);
	template std::ostream& operator<<(std::ostream&, const Point<DOW>&);

	template filtering_istream& operator>>(filtering_istream&, QuadratureInfo<DOW>&);
	template std::ostream& operator<<(std::ostream&, const QuadratureInfo<DOW>&);

	template filtering_istream& operator>>(filtering_istream&, QuadratureInfoAdmin<DOW>&);
	template std::ostream& operator<<(std::ostream&, const QuadratureInfoAdmin<DOW>&);

	template filtering_istream& operator>>(filtering_istream&, TemplateGeometry<DOW>&);
	template std::ostream& operator<<(std::ostream&, const TemplateGeometry<DOW>&);

#define DIM 1
	template class Mesh<DIM,DOW>;
	template class SimplestMesh<DIM,DOW>;

	template std::istream& operator>>(std::istream&, Mesh<DIM,DOW>&);
	template std::ostream& operator<<(std::ostream&, const Mesh<DIM,DOW>&);
#undef DIM
	
#define DIM 2
	template class Mesh<DIM,DOW>;
	template class SimplestMesh<DIM,DOW>;

	template std::istream& operator>>(std::istream&, Mesh<DIM,DOW>&);
	template std::ostream& operator<<(std::ostream&, const Mesh<DIM,DOW>&);
#undef DIM
	
#define DIM 3
	template class Mesh<DIM,DOW>;
	template class SimplestMesh<DIM,DOW>;

	template std::istream& operator>>(std::istream&, Mesh<DIM,DOW>&);
	template std::ostream& operator<<(std::ostream&, const Mesh<DIM,DOW>&);
#undef DIM
#undef DOW

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
