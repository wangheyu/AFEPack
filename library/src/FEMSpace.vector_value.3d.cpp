/**
 * @file   FEMSpace.vector_value.3d.cpp
 * @author Robert Lie
 * @date   Thu Mar 30 11:32:08 2006
 * 
 * @brief  vector valued finite element space for DOW=3
 * 
 * 
 */

#include "FEMSpace.templates.h"

#define DOW 3
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
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

#undef TDIM
#undef DIM

#define DIM 2
#define TDIM 2

#define vector_length 1
#include "FEMSpace.vector_value.templates.h"
#undef _FEMSpace_vector_value_templates_h_

AFEPACK_OPEN_NAMESPACE

#define value_type nVector<vector_length,double>
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
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
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

#undef TDIM
#undef DIM

#define DIM 3
#define TDIM 3

#define vector_length 1
#include "FEMSpace.vector_value.templates.h"
#undef _FEMSpace_vector_value_templates_h_

AFEPACK_OPEN_NAMESPACE

#define value_type nVector<vector_length,double>
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
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
#undef value_type
#undef vector_length

AFEPACK_CLOSE_NAMESPACE

#undef TDIM
#undef DIM
#undef DOW

#define DIM 3
#define vector_length   DIM
#define value_type        nVector<vector_length,double>
#define vector_zero       value_type(0.0)

AFEPACK_OPEN_NAMESPACE

template <>
void FEMFunction<value_type, DIM>::writeTecplotData(const std::string& filename)
{
  std::ofstream os(filename.c_str());

  os << "Variables = "
     << "\x22" << "X\x22, "
     << "\x22" << "Y\x22, "
     << "\x22" << "Z\x22, ";
  for (unsigned int k = 0;k < vector_length;k ++) {
    os << "\x22" << "U" << k << "\x22";
  }
  os << "\n";
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  FEMSpace<value_type, DIM>& fem_space = femSpace();
  unsigned int n_dof = fem_space.n_dof();
  std::vector<bool> flag(n_dof, false);
  FEMSpace<value_type, DIM>::ElementIterator the_element = fem_space.beginElement();
  FEMSpace<value_type, DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      if (flag[element_dof[i]]) continue;
      const Point<DIM>& p = fem_space.dofInfo(element_dof[i]).interp_point;
      value_type v = value(p, *the_element);
      os << p << "\t";
      for (unsigned int j = 0;j < vector_length;j ++)
	 os << v[j] << "\n";
      flag[element_dof[i]] = true;
    }
  }
  os.close();
};

template <>
void FEMFunction<value_type,DIM>::writeOpenDXData(const std::string& filename,
						  int flag) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(12);
  os.setf(std::ios::fixed, std::ios::floatfield);
  const FEMSpace<value_type,DIM>& fem_space = femSpace();
  const Mesh<DIM,DIM>& mesh = fem_space.mesh();
  int n_node = mesh.n_point();
  std::vector<int> count(n_node, 0);
  std::vector<value_type> val(n_node, vector_zero);
  int i, j, k;
  FEMSpace<value_type,DIM>::ConstElementIterator the_element = fem_space.beginElement();
  FEMSpace<value_type,DIM>::ConstElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const GeometryBM& geo = the_element->geometry();
    for (i = 0;i < geo.n_vertex();i ++) {
      k = mesh.geometry(0, geo.vertex(i)).vertex(0);
      count[k] += 1;
      value_type v = value(mesh.point(k), *the_element);
      for (j = 0;j < vector_length;j ++)
	val[geo.vertex(i)][j] += v[j];
    }
  }
  for (i = 0;i < n_node;i ++) {
    Assert(count[i] > 0, ExcInternalError());
    for (j = 0;j < vector_length;j ++)
      val[i][j] /= count[i];
  }
	
  os << "object 1 class array type float rank 1 shape 3 item " 
     << n_node << " data follows\n";
  for (i = 0;i < n_node;i ++) {
    os << mesh.point(i) << "\n";
  }
	
  int n_element = mesh.n_geometry(DIM);
  for (i = 0, j = 0;i < n_element;i ++) {
    switch (mesh.geometry(DIM,i).n_vertex()) {
    case 4:
      j += 1;
      break;
    case 5:
      j += 2;
      break;
    case 7:
      j += 4;
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << "\nobject 2 class array type int rank 1 shape 4 item "
     << j << " data follows\n";
  for (i = 0;i < n_element;i ++) {
    switch (mesh.geometry(DIM,i).n_vertex()) {
    case 4:
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(3)).vertex(0) << "\t\n";
      break;
    case 5:
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(4)).vertex(0) << "\t\n";
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(3)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(4)).vertex(0) << "\t\n";
      break;
    case 7:
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(6)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(5)).vertex(0) << "\t\n";
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(4)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(6)).vertex(0) << "\t\n";
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(3)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(5)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(4)).vertex(0) << "\t\n";
      os << mesh.geometry(0, mesh.geometry(3,i).vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(4)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(5)).vertex(0) << "\t"
	 << mesh.geometry(0, mesh.geometry(3,i).vertex(6)).vertex(0) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << "attribute \"element type\" string \"tetrahedra\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
  os << "object 3 class array type float rank 1 shape "
     << vector_length << " item "
     << n_node << " data follows\n";
  for (i = 0;i < n_node;i ++) {
    for (j = 0;j < vector_length;j ++) {
      os << val[i][j] << "\t";
    }
    os << "\n";
  }
  os << "attribute \"dep\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction-3d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "component \"data\" value 3\n"
     << "end\n";
  os.close();
};

AFEPACK_CLOSE_NAMESPACE

#undef vector_zero
#undef value_type
#undef vector_length
#undef DIM

/**
 * end of file
 * 
 */
