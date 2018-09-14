/**
 * @file   FEMSpace.cpp
 * @author Robert Lie
 * @date   Wed Jun 14 12:06:24 2006
 * 
 * @brief  
 * 
 * 
 */


#include "FEMSpace.templates.h"

AFEPACK_OPEN_NAMESPACE

BoundaryConditionInfo::Type BoundaryConditionInfo::DIRICHLET		= 1;
BoundaryConditionInfo::Type BoundaryConditionInfo::NEUMANN		= 2;
BoundaryConditionInfo::Type BoundaryConditionInfo::ROBIN		= 3;


template <>
void FEMFunction<double, 2>::writeEasyMeshData(const std::string& filename)
{
  femSpace().mesh().writeEasyMesh(filename);
  std::ofstream os((filename + ".dat").c_str());
	
  os.precision(12);
  os.setf(std::ios::fixed, std::ios::floatfield);
  FEMSpace<double,2>& fem_space = femSpace();
  Mesh<2,2>& mesh = fem_space.mesh();
  int n_node = mesh.n_point();
  std::vector<int> count(n_node, 0);
  std::vector<double> val(n_node, 0);
  int i, j, k;
  FEMSpace<double,2>::ElementIterator the_element = fem_space.beginElement();
  FEMSpace<double,2>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    GeometryBM& geo = the_element->geometry();
    for (i = 0;i < geo.n_vertex();i ++) {
      k = mesh.geometry(0, geo.vertex(i)).vertex(0);
      count[k] += 1;
      val[k] += value(mesh.point(k), *the_element);
    }
  }
  for (i = 0;i < n_node;i ++) {
    Assert(count[i] > 0, ExcInternalError());
    val[i] /= count[i];
  }
  for (i = 0;i < n_node;i ++)	os << val[i] << "\n";
  os.close();
}
	
template <>
void FEMFunction<double, 2>::writeTecplotData(const std::string& filename)
{
  std::ofstream os(filename.c_str());

  os << "Variables = "
     << "\x22" << "X\x22, "
     << "\x22" << "Y\x22, "
     << "\x22" << "U\x22\n";
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  FEMSpace<double, 2>& fem_space = femSpace();
  unsigned int n_dof = fem_space.n_dof();
  std::vector<bool> flag(n_dof, false);
  FEMSpace<double, 2>::ElementIterator the_element = fem_space.beginElement();
  FEMSpace<double, 2>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      if (flag[element_dof[i]]) continue;
      const Point<2>& p = fem_space.dofInfo(element_dof[i]).interp_point;
      double v = value(p, *the_element);
      os << p << "\t"
	 << v << "\n";
      flag[element_dof[i]] = true;
    }
  }
  os.close();
}



template <>
void FEMFunction<double, 3>::writeTecplotData(const std::string& filename)
{
  std::ofstream os(filename.c_str());

  os << "Variables = "
     << "\x22" << "X\x22, "
     << "\x22" << "Y\x22, "
     << "\x22" << "Z\x22, "
     << "\x22" << "U\x22\n";
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  FEMSpace<double, 3>& fem_space = femSpace();
  unsigned int n_dof = fem_space.n_dof();
  std::vector<bool> flag(n_dof, false);
  FEMSpace<double, 3>::ElementIterator the_element = fem_space.beginElement();
  FEMSpace<double, 3>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      if (flag[element_dof[i]]) continue;
      const Point<3>& p = fem_space.dofInfo(element_dof[i]).interp_point;
      double v = value(p, *the_element);
      os << p << "\t"
	 << v << "\n";
      flag[element_dof[i]] = true;
    }
  }
  os.close();
}

template <>
void FEMFunction<double, 1>::writeOpenDXData(const std::string& filename,
					     int flag) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(12);
  os.setf(std::ios::fixed, std::ios::floatfield);
  const FEMSpace<double,1>& fem_space = femSpace();
  const Mesh<1,1>& mesh = fem_space.mesh();
  int n_node = mesh.n_point();
  os << "object 1 class array type float rank 1 shape 2 item " 
     << 2*n_node << " data follows\n";
  for (int i = 0;i < n_node;i ++) {
    os << mesh.point(i) << "\t0.0\n"
       << mesh.point(i) << "\t1.0\n";
  }

  int n_element = fem_space.n_element();
  int n_data = n_element;
  if (flag == 0) n_data = 2*n_node; // position dependent data
  std::vector<double> val(n_data, 0);

  os << "\nobject 2 class array type int rank 1 shape 4 item "
     << n_element << " data follows\n";
  FEMSpace<double,1>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  for (int i = 0;the_element != end_element;++ the_element,++ i) {
    GeometryBM& geo = the_element->geometry();
    os << 2*geo.vertex(0) << "\t" << 2*geo.vertex(0) + 1 << "\t"
       << 2*geo.vertex(1) << "\t" << 2*geo.vertex(1) + 1 << "\n";
    double v0 = value(mesh.point(geo.vertex(0)), *the_element);
    double v1 = value(mesh.point(geo.vertex(1)), *the_element);
    if (flag == 1) {
      val[i] = 0.5*(v0 + v1);
    } else {
      val[2*geo.vertex(0)] = v0; val[2*geo.vertex(0) + 1] = v0;
      val[2*geo.vertex(1)] = v1; val[2*geo.vertex(1) + 1] = v1;
    }
  }
  os << "attribute \"element type\" string \"quads\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object 3 class array type float rank 0 item "
     << n_data << " data follows\n";
  for (int i = 0;i < n_data;i ++) os << val[i] << "\n";
  if (flag == 0)
    os << "attribute \"dep\" string \"positions\"\n\n";
  else if (flag == 1)
    os << "attribute \"dep\" string \"connections\"\n\n";
	
  os << "object \"FEMFunction-1d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "component \"data\" value 3\n"
     << "end\n";
  os.close();
}


template <>
void FEMFunction<double, 2>::writeOpenDXData(const std::string& filename,
					     int flag) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(12);
  os.setf(std::ios::fixed, std::ios::floatfield);
  const FEMSpace<double,2>& fem_space = femSpace();
  const Mesh<2,2>& mesh = fem_space.mesh();
  int n_node = mesh.n_point();
  os << "object 1 class array type float rank 1 shape 2 item " 
     << n_node << " data follows\n";
  for (int i = 0;i < n_node;i ++) {
    os << mesh.point(i) << "\n";
  }

  int n_element = 0;
  FEMSpace<double,2>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  for (;the_element != end_element;++ the_element) {
    const GeometryBM& geo = the_element->geometry();
    switch (geo.n_vertex()) {
    case 3: n_element += 1; break;
    case 4: n_element += 2; break;
    }
  }

  int n_data = n_element;
  if (flag == 0) n_data = n_node; // position dependent data
  std::vector<int> count(n_data, 0);
  std::vector<double> val(n_data, 0);

  os << "\nobject 2 class array type int rank 1 shape 3 item "
     << n_element << " data follows\n";
  the_element = fem_space.beginElement();
  for (int j = 0;the_element != end_element;++ the_element) {
    const GeometryBM& geo = the_element->geometry();
    switch (geo.n_vertex()) {
    case 3:
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t\n";
      break;
    case 4:
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t\n";
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(3)).vertex(0) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
    if (flag == 0) {
      for (int i = 0;i < geo.n_vertex();i ++) {
	int k = mesh.geometry(0, geo.vertex(i)).vertex(0);
	count[k] += 1;
	val[k] += value(mesh.point(k), *the_element);
      }
    } else if (flag == 1) {
      Point<2> bc = (mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0)) + 
                     mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0)) +
                     mesh.point(mesh.geometry(0, geo.vertex(2)).vertex(0)));
      bc[0] /= 3; bc[1] /= 3;
      val[j ++] = value(bc, *the_element);

      if (geo.n_vertex() == 4) {
	bc = (mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0)) +
              mesh.point(mesh.geometry(0, geo.vertex(2)).vertex(0)) +
	      mesh.point(mesh.geometry(0, geo.vertex(3)).vertex(0)));
	bc[0] /= 3; bc[1] /= 3;
	val[j ++] = value(bc, *the_element);
      }
    }
  }
  os << "attribute \"element type\" string \"triangles\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  if (flag == 0) {
    for (int i = 0;i < n_node;i ++) {
      Assert(count[i] > 0, ExcInternalError());
      val[i] /= count[i];
    }
  }
	
  os << "object 3 class array type float rank 0 item "
     << n_data << " data follows\n";
  for (int i = 0;i < n_data;i ++) os << val[i] << "\n";
  if (flag == 0)
    os << "attribute \"dep\" string \"positions\"\n\n";
  else if (flag == 1)
    os << "attribute \"dep\" string \"connections\"\n\n";
	
  os << "object \"FEMFunction-2d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "component \"data\" value 3\n"
     << "end\n";
  os.close();
}


template <>
void FEMFunction<double,3>::writeOpenDXData(const std::string& filename,
					    int flag) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(12);
  os.setf(std::ios::fixed, std::ios::floatfield);
  const FEMSpace<double,3>& fem_space = femSpace();
  const Mesh<3,3>& mesh = fem_space.mesh();
  int n_node = mesh.n_point();
  os << "object 1 class array type float rank 1 shape 3 item " 
     << n_node << " data follows\n";
  for (int i = 0;i < n_node;i ++) {
    os << mesh.point(i) << "\n";
  }
  
  int n_element = 0;
  FEMSpace<double,3>::ConstElementIterator 
    the_element = fem_space.beginElement(),
    end_element = fem_space.endElement();
  for (;the_element != end_element;++ the_element) {
    const GeometryBM& geo = the_element->geometry();
    switch (geo.n_vertex()) {
    case 4: n_element += 1; break;
    case 5: n_element += 2; break;
    case 7: n_element += 4; break;
    }
  }

  int n_data = n_element;
  if (flag == 0) n_data = n_node;
  std::vector<int> count(n_data, 0);
  std::vector<double> val(n_data, 0);

  os << "\nobject 2 class array type int rank 1 shape 4 item "
     << n_element << " data follows\n";
  the_element = fem_space.beginElement();
  for (int j = 0;the_element != end_element;++ the_element) {
    const GeometryBM& geo = the_element->geometry();
    switch (geo.n_vertex()) {
    case 4:
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(3)).vertex(0) << "\t\n";
      break;
    case 5:
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(4)).vertex(0) << "\t\n";
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(3)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(4)).vertex(0) << "\t\n";
      break;
    case 7:
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(1)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(6)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(5)).vertex(0) << "\t\n";
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(2)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(4)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(6)).vertex(0) << "\t\n";
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(3)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(5)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(4)).vertex(0) << "\t\n";
      os << mesh.geometry(0, geo.vertex(0)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(4)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(5)).vertex(0) << "\t"
	 << mesh.geometry(0, geo.vertex(6)).vertex(0) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
    if (flag == 0) {
      for (int i = 0;i < geo.n_vertex();i ++) {
        int k = mesh.geometry(0, geo.vertex(i)).vertex(0);
        count[k] += 1;
        val[k] += value(mesh.point(k), *the_element);
      }
    } else if (flag == 1) {
      int k = mesh.geometry(0, geo.vertex(0)).vertex(0);
      double v = value(mesh.point(k), *the_element);
      switch (geo.n_vertex()) {
      case 4: val[j ++] = v; break;
      case 5: val[j ++] = v; val[j ++] = v; break;
      case 7: val[j ++] = v; val[j ++] = v; 
              val[j ++] = v; val[j ++] = v; break;
      }
    }
  }
  os << "attribute \"element type\" string \"tetrahedra\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  if (flag == 0) {
    for (int i = 0;i < n_data;i ++) {
      Assert(count[i] > 0, ExcInternalError());
      val[i] /= count[i];
    }
  }

  os << "object 3 class array type float rank 0 item "
     << n_data << " data follows\n";
  for (int i = 0;i < n_data;i ++) os << val[i] << "\n";
  if (flag == 0)
    os << "attribute \"dep\" string \"positions\"\n\n";
  else if (flag == 1)
    os << "attribute \"dep\" string \"connections\"\n\n";
	
  os << "object \"FEMFunction-3d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "component \"data\" value 3\n"
     << "end\n";
  os.close();
}


/////////////////////////////////////////////////////////////////////

#define value_type double
#define DOW 1
#define DIM 1
#define TDIM 1
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
template class LocalFEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef TDIM
#undef DIM
template struct ElementAdditionalData<value_type,DOW>;
#undef DOW

/////////////////////////////////////////////////////////////////////

#define DOW 2
#define DIM 1
#define TDIM 1
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
template class LocalFEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef TDIM
#undef DIM

#define DIM 2
#define TDIM 2
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
template class LocalFEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef TDIM
#undef DIM
template struct ElementAdditionalData<value_type,DOW>;
#undef DOW

/////////////////////////////////////////////////////////////////////

#define DOW 3
#define DIM 1
#define TDIM 1
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
template class LocalFEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef TDIM
#undef DIM

#define DIM 2
#define TDIM 2
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
template class LocalFEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef TDIM
#undef DIM

#define DIM 3
#define TDIM 3
template class Element<value_type,DIM,DOW,TDIM>;
template class FEMSpace<value_type,DIM,DOW,TDIM>;
template class FEMFunction<value_type,DIM,DOW,TDIM>;
template class LocalFEMFunction<value_type,DIM,DOW,TDIM>;

template class BoundaryCondition<value_type,DIM,DOW,TDIM>;
template class BoundaryConditionAdmin<value_type,DIM,DOW,TDIM>;
#undef TDIM
#undef DIM
template struct ElementAdditionalData<value_type,DOW>;
#undef DOW
#undef value_type

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

