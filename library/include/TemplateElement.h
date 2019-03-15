///////////////////////////////////////////////////////////////////////////
// TemplateElement.h : by R.Lie
//

#ifndef _TemplateElement_h_
#define _TemplateElement_h_

#include <dlfcn.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#include <base/exceptions.h>

#include "Miscellaneous.h"
#include "Geometry.h"

AFEPACK_OPEN_NAMESPACE

template <int DIM> class BasisFunctionIdentity;
template <class value_type, int DIM> class ShapeFunction;
template <class value_type, int DIM, int TDIM> class BasisFunction;
template <int DIM> class TemplateDOF;
template <int TDIM, int DIM> class CoordTransform;
template <class value_type, int DIM, int TDIM> class BasisFunctionAdmin;
template <int DIM> class UnitOutNormal;
template <class value_type, int DIM, int TDIM> class TemplateElement;

template <int DIM> bool operator==(const BasisFunctionIdentity<DIM>&, const BasisFunctionIdentity<DIM>&);

///////////////////////////////////////////////////////////////////////////

/**
 * Basis function identity information. The shape functions on neighbour
 * elements are combined together as basis function according if their
 * identities are the same. It includes three information: polynomial order
 * on the boundary of the element; interpolation operator; additional flag.
 * We call the meaning of this class "Basis Function Identity Protocol(BFIP)",
 * while it's not a very strict protocol. You can appoint the value of the
 * members as your willing, bearing the mind that those value should be
 * coincide for different shape functions of the same basis function.
 */
template <int DIM>
class BasisFunctionIdentity 
{
 public:
  enum { dim = DIM };
 public:
  unsigned int order; /**< Polynomial order. */
  int alpha[DIM]; /**< Interpolation operator. */
  unsigned int flag; /**< Additional flag. */
 public:
  /** Default contructor. */
  BasisFunctionIdentity() {};
  /** Copy contructor. */
  BasisFunctionIdentity(const BasisFunctionIdentity<DIM>& i) {
    order = i.order;
    for (int j = 0;j < DIM;j ++)
      alpha[j] = i.alpha[j];
    flag = i.flag;
  };
  /** Destructor. */
  ~BasisFunctionIdentity() {};
 public:
  /** Copy operator. */
  BasisFunctionIdentity<DIM>& operator=(const BasisFunctionIdentity<DIM>& i) {
    order = i.order;
    for (int j = 0;j < DIM;j ++)
      alpha[j] = i.alpha[j];
    flag = i.flag;
    return *this;
  };
  friend bool operator== <>(const BasisFunctionIdentity<DIM>&, const BasisFunctionIdentity<DIM>&); /**< Judge if two basis function identities are equal. */
  template <class STREAM,int GDIM> 
    friend STREAM& operator>>(STREAM&, BasisFunctionIdentity<GDIM>&);
  template <class STREAM,int GDIM> 
    friend STREAM& operator<<(STREAM&, const BasisFunctionIdentity<GDIM>&);
};

template <class STREAM, int DIM>
  STREAM& operator>>(STREAM& is, BasisFunctionIdentity<DIM>& i)
{
  is >> i.order;
  for (int j = 0;j < DIM;j ++)
    is >> i.alpha[j];
  is >> i.flag;
  return is;
}

template <class STREAM, int DIM>
  STREAM& operator<<(STREAM& os, const BasisFunctionIdentity<DIM>& i)
{
  os << i.order;
  for (int j = 0;j < DIM;j ++)
    os << i.alpha[j];
  os << i.flag;
  return os;
}

/**
 * Shape function on an element. This class will mainly used to
 */
template <class value_type, int DIM>
  class ShapeFunction
{
 public:
  typedef value_type value_t;
  enum { dim = DIM };
 private:
  void * handle; /**< Handle of object to operate on the shared library. */
  std::string library_name; /**< Shared library name. */
  std::string value_function_name; /**< Name of value function. */
  std::string gradient_function_name; /**< Name of gradient function. */
  void (*value_function)(const double *, const double **, void *); /**< Pointer to value function. */
  void (*gradient_function)(const double *, const double **, void *); /**< Pointer to gradient function. */
 public:
  std::string library_path;
 public:
  ShapeFunction(); /**< Default contructor. */
  ShapeFunction(const ShapeFunction<value_type,DIM>&); /**< Copy contructor. */
  ~ShapeFunction(); /**< Destructor. */
 public:
  ShapeFunction<value_type,DIM>& operator=(const ShapeFunction<value_type,DIM>&); /**< Copy operator. */
  void loadFunction(); /**< Open the shared library and load the functions. */
  void unloadFunction(); /**< Close the shared library. */
  value_type value(const Point<DIM>&, const std::vector<Point<DIM> >&) const; /**< Value of the basis function on a point. */
  std::vector<value_type> gradient(const Point<DIM>&, const std::vector<Point<DIM> >&) const; /**< Gradient of the basis function on a point. */
  std::vector<value_type> value(const std::vector<Point<DIM> >&, const std::vector<Point<DIM> >&) const; /**< Value of the basis function on points. */
  std::vector<std::vector<value_type> > gradient(const std::vector<Point<DIM> >&, const std::vector<Point<DIM> >&) const; /**< Value of the basis function on points. */
  value_type value(const Point<DIM>&, const double **) const; /**< Value of the basis function on a point. */
  std::vector<value_type> gradient(const Point<DIM>&, const double **) const; /**< Gradient of the basis function on a point. */
  std::vector<value_type> value(const std::vector<Point<DIM> >&, const double **) const; /**< Value of the basis function on points. */
  std::vector<std::vector<value_type> > gradient(const std::vector<Point<DIM> >&, const double **) const; /**< Gradient of the basis function on points. */

#ifdef QUADRATIC_ELEMENT_SUPPORT
 private:
  std::string hesse_function_name;
  void (*hesse_function)(const double *, const double **, void *);
 public:
  std::vector<std::vector<value_type> > 
    hesse(const Point<DIM>&, const std::vector<Point<DIM> >&) const;
  std::vector<std::vector<std::vector<value_type> > > 
    hesse(const std::vector<Point<DIM> >&, const std::vector<Point<DIM> >&) const;
  std::vector<std::vector<value_type> > 
    hesse(const Point<DIM>&, const double **) const;
  std::vector<std::vector<std::vector<value_type> > 
    hesse(const std::vector<Point<DIM> >&, const double **) const;
#endif //QUADRATIC_ELEMENT_SUPPORT

  DeclException1(ExcFileOpen, char *,
                 << "Can't open library "
                 << arg1);
  DeclException2(ExcLoadFunction, char *, char *,
                 << "Can't load function "
                 << arg1
                 << " from library "
                 << arg2);	
 public:
  template <class STREAM,class VT,int GDIM> 
    friend STREAM& operator>>(STREAM&, ShapeFunction<VT,GDIM>&); /**< Stream input. */
  template <class STREAM,class VT,int GDIM> 
    friend STREAM& operator<<(STREAM&, const ShapeFunction<VT,GDIM>&); /**< Stream output. */
};

template <class STREAM, class value_type, int DIM>
  STREAM& operator>>(STREAM& is, ShapeFunction<value_type,DIM>& f)
{
  is >> f.library_name;
  is >> f.value_function_name;
  is >> f.gradient_function_name;
#ifdef QUADRATIC_ELEMENT_SUPPORT
  is >> f.hesse_function_name;
#endif // QUADRATIC_ELEMENT_SUPPORT
  f.loadFunction();

  return is;
}

template <class STREAM, class value_type, int DIM>
  STREAM& operator<<(STREAM& os, const ShapeFunction<value_type,DIM>& f)
{
  os << f.library_name << "\t"
     << f.value_function_name << "\t"
     << f.gradient_function_name 
#ifdef QUADRATIC_ELEMENT_SUPPORT
     << "\t" << f.hesse_function_name
#endif // QUADRATIC_ELEMENT_SUPPORT
     << "\n";
  return os;
}

/**
 * Basis function
 */
template <class value_type, int DIM, int TDIM=DIM>
  class BasisFunction : public ShapeFunction<value_type, DIM>
  {
  public:
  typedef value_type value_t;
  enum { dim = DIM, tdim = TDIM };
  public:
  typedef BasisFunctionIdentity<TDIM> Identity;
  private:
  //	TemplateGeometry<TDIM> * geo;
  Point<TDIM> ip;
  Identity id;
  public:
  BasisFunction();
  explicit BasisFunction(const Point<TDIM>&);
  BasisFunction(const BasisFunction<value_type,DIM,TDIM>&);
  ~BasisFunction();
  public:
  BasisFunction<value_type,DIM,TDIM>& operator=(const BasisFunction<value_type,DIM,TDIM>&);
  void reinit(const Point<TDIM>&);
  const Point<TDIM>& interpPoint() const;
  Point<TDIM>& interpPoint();
  const Identity& identity() const;
  Identity& identity();
  public:
  template <class STREAM,class VT,int GDIM,int TGDIM> 
  friend STREAM& operator>>(STREAM&, BasisFunction<VT,GDIM,TGDIM>&);
  template <class STREAM,class VT,int GDIM,int TGDIM> 
  friend STREAM& operator<<(STREAM&, const BasisFunction<VT,GDIM,TGDIM>&);
  };

template <class STREAM, class value_type, int DIM, int TDIM> 
  STREAM& operator>>(STREAM& is, BasisFunction<value_type,DIM,TDIM>& b)
{
  is >> b.interpPoint();
  is >> b.identity();
  is >> dynamic_cast<ShapeFunction<value_type,DIM>&>(b);
  return is;
}

template <class STREAM, class value_type, int DIM, int TDIM> 
  STREAM& operator<<(STREAM& os, const BasisFunction<value_type,DIM,TDIM>& b)
{
  os << b.interpPoint() << "\n";
  os << b.identity() << "\n";
  os << dynamic_cast<const ShapeFunction<value_type,DIM>&>(b);
  return os;
}

/**
 * Index of degree of freedom to identify the relationship of a degree of freedom to
 * certain geometry.
 */
struct DOFIndex {
  int dimension; /**< Dimension of geometry associated with the degree of freedom. */
  int geometry_index; /**< Index of geometry associated with the degree of freedom. */
  int dof_index; /**< Local index of the degree of freedom on the geometry. */
};

template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM>
  struct DOFInfo {
    enum { dim = DIM, dow = DOW, tdim = TDIM };
    typedef value_type value_t;

    typedef typename Mesh<DIM,DOW>::bmark_t bmark_t;
    Point<DOW> interp_point;
    typename BasisFunction<value_type,DIM,TDIM>::Identity identity;
    bmark_t boundary_mark;
  };

/**
 * This struct is in fact the relationship of degrees of freedom to geometries. From
 * the information in this struct, you can retrieve the geometry information of certain
 * degree of freedom, and get the degree of freedom on certain geometry.
 */
struct DegreeOfFreedom
{
  unsigned int n_dof; /**< Number of degree of freedoms. */
  std::vector<std::vector<int> > n_geometry_dof; /**< Number of degree of freedoms on certain geometry. */
  std::vector<std::vector<std::vector<int> > > geometry_dof; /**< Array of degree of freedom on certain geometry. */
  std::vector<DOFIndex> dof_index; /**< Geometry index of degree of freedoms. */
};

/**
 * The degree of freedom distribution on certain template geometry. This class builds the
 * relationship of some local degree of freedom and the geometries of a template element
 * geometry. With this information, the package then can build the degree of freedom on
 * the whole mesh.
 */
template <int DIM>
class TemplateDOF : public DegreeOfFreedom
{
 public:
  enum { dim = DIM };
 private:
  TemplateGeometry<DIM> * geometry; /**< Pointer to the template element geometry. */
 public:
  TemplateDOF() : geometry(NULL) {}
  TemplateDOF(TemplateGeometry<DIM>&); /**< Default constructor. */
  TemplateDOF(const TemplateDOF<DIM>&); /**< Copy Contructor. */
  ~TemplateDOF(); /**< Destructor. */
 public:
  TemplateDOF<DIM>& operator=(const TemplateDOF<DIM>&); /**< Copy operator. */
  void reinit(TemplateGeometry<DIM>&); /**< Reinitialization. */
  void readData(const std::string&); /**< Read in the information from a file. */
  void writeData(const std::string&) const; /**< Write out the information to a file. */
 public:
  template <class STREAM,int GDIM> 
    friend STREAM& operator>>(STREAM&, TemplateDOF<GDIM>&); /**< Stream input. */
  template <class STREAM,int GDIM> 
    friend STREAM& operator<<(STREAM&, const TemplateDOF<GDIM>&); /**< Stream output. */
};

template <class STREAM, int DIM>
  STREAM& operator>>(STREAM& is, TemplateDOF<DIM>& t)
{
  int i, j, k, l, m;
	
  for (i = 0;i <= DIM;i ++) {
    for (j = 0;j < t.geometry->n_geometry(i);j ++) {
      t.n_geometry_dof[i][j] = 0;
    }
  }
  is >> i;
  for (t.n_dof = 0, j = 0;j < i;j ++) {
    is >> k >> l >> m;
    t.n_geometry_dof[k][l] += m;
    t.n_dof += m;
  }
  t.dof_index.resize(t.n_dof);
  for (k = 0, i = 0;i <= DIM;i ++) {
    for (j = 0;j < t.geometry->n_geometry(i);j ++) {
      t.geometry_dof[i][j].resize(t.n_geometry_dof[i][j]);
      for (l = 0;l < t.n_geometry_dof[i][j];l ++) {
	t.dof_index[k].dimension = i;
	t.dof_index[k].geometry_index = j;
	t.dof_index[k].dof_index = l;
	t.geometry_dof[i][j][l] = k ++;
      }
    }
  }
  return is;
}

template <class STREAM, int DIM>
  STREAM& operator<<(STREAM& os, const TemplateDOF<DIM>& t)
{
  int i, j, k;
	
  for (k = 0, i = 0;i <= DIM;i ++) {
    for (j = 0;j < t.geometry->n_geometry(i);j ++) {
      if (t.n_geometry_dof[i][j] > 0)
	k ++;
    }
  }
  os << k << "\n";
  for (i = 0;i <= DIM;i ++) {
    for (j = 0;j < t.geometry->n_geometry(i);j ++) {
      if (t.n_geometry_dof[i][j] > 0) {
	os << i << "\t"
	   << j << "\t"
	   << t.n_geometry_dof[i][j] << "\n";
      }
    }
  }
  return os;
}

/**
 * Coordinate transformation between the template element geometry and real element geometry.
 * This class only contains the information of several functions to implement such operations.
 * Those functions serious depend on the order of the vertices of the geometries. The facilities
 * provided including: get the corresponding point on the real element of a point on the
 * template element and vice visa, get the Jacobian determinant at a certain point on
 * the template element or on the real element.
 */
template <int TDIM, int DIM=TDIM>
  class CoordTransform
  {
  public:
  enum { dim = DIM, tdim = TDIM };
  typedef Point<DIM> point_t;
  typedef Point<TDIM> ref_point_t;
  private:
  void * handle; /**< Handle of the object to open the shared library. */
  std::string library_path;
  std::string library_name; /**< Name of the shared library. */
  std::string l2g_function_name; /**< Name of the function to map a point on template element to real element. */
  std::string g2l_function_name; /**< Name of the function to map a point on real element to template element. */
  std::string l2g_jacobian_function_name; /**< Jacobian determinant of the function at a point on template element. */
  std::string g2l_jacobian_function_name; /**< Jacobian determinant of the function at a point on real element. */
  void (*l2g_function)(const double *, const double **, const double **, double *); /**< Pointer to the function to map a point on template element to real element. */
  void (*g2l_function)(const double *, const double **, const double **, double *); /**< Pointer to the function to map a point on real element to template element. */
  double (*l2g_jacobian_function)(const double *, const double **, const double **); /**< Pointer to the function given the Jacobian determinant at a point on template element. */
  double (*g2l_jacobian_function)(const double *, const double **, const double **); /**< Pointer to the function given the Jacobian determinant at a point on real element. */
  public:
  CoordTransform(); /**< Default constructor. */
  CoordTransform(const CoordTransform<TDIM,DIM>&); /**< Copy constructor. */
  ~CoordTransform(); /**< Destructor. */
  public:
  CoordTransform<TDIM,DIM>& operator=(const CoordTransform<TDIM,DIM>&); /**< Copy operator. */
  void loadFunction(); /**< Open the shared library and get those function pointers. */
  void unloadFunction(); /**< Close the shared library. */
  point_t local_to_global(const ref_point_t&, 
                          const std::vector<ref_point_t >&, 
                          const std::vector<point_t >&) const; /**< Mapping a point on template element to real element. */
  ref_point_t global_to_local(const point_t&, 
                              const std::vector<ref_point_t >&, 
                              const std::vector<point_t >&) const; /**< Mapping a point on real element to template element. */
  double local_to_global_jacobian(const ref_point_t&, 
                                  const std::vector<ref_point_t >&, 
                                  const std::vector<point_t >&) const; /**< Jacobian determinant at a point on template element. */
  double global_to_local_jacobian(const point_t&,
                                  const std::vector<ref_point_t >&, 
                                  const std::vector<point_t >&) const; /**< Jacobian determinant at a point on real element. */
  std::vector<point_t > local_to_global(const std::vector<ref_point_t >&, 
                                        const std::vector<ref_point_t >&, 
                                        const std::vector<point_t >&) const; /**< Mapping points on template element to real element. */
  std::vector<ref_point_t > global_to_local(const std::vector<point_t >&, 
                                            const std::vector<ref_point_t >&, 
                                            const std::vector<point_t >&) const; /**< Mapping points on real element to template element. */
  std::vector<double> local_to_global_jacobian(const std::vector<ref_point_t >&, 
                                               const std::vector<ref_point_t >&, 
                                               const std::vector<point_t >&) const; /**< Jacobian determinant at points on template element. */
  std::vector<double> global_to_local_jacobian(const std::vector<point_t >&,
                                               const std::vector<ref_point_t >&, 
                                               const std::vector<point_t >&) const; /**< Jacobian determinant at points on real element. */
  void readData(const std::string&); /**< Read in the information from a file. */
  void writeData(const std::string&) const; /**< Write out the information to a file. */

  DeclException1(ExcFileOpen, char *,
                 << "Can't open library "
                 << arg1);
  DeclException2(ExcLoadFunction, char *, char *,
                 << "Can't load function "
                 << arg1
                 << " from library "
                 << arg2);	
  public:
  template <class STREAM,int TGDIM,int GDIM> 
  friend STREAM& operator>>(STREAM&, CoordTransform<TGDIM,GDIM>&); /**< Stream input. */
  template <class STREAM,int TGDIM,int GDIM> 
  friend STREAM& operator<<(STREAM&, const CoordTransform<TGDIM,GDIM>&); /**< Stream output. */
  };

template <class STREAM, int TDIM, int DIM>
  STREAM& operator>>(STREAM& is, CoordTransform<TDIM,DIM>& c)
{
  is >> c.library_name
     >> c.l2g_function_name
     >> c.g2l_function_name
     >> c.l2g_jacobian_function_name
     >> c.g2l_jacobian_function_name;
	
  c.loadFunction();
  return is;
}

template <class STREAM, int TDIM, int DIM>
  STREAM& operator<<(STREAM& os, const CoordTransform<TDIM,DIM>& c)
{
  os << c.library_name << "\n\t"
     << c.l2g_function_name << "\t"
     << c.g2l_function_name << "\t"
     << c.l2g_jacobian_function_name << "\t"
     << c.g2l_jacobian_function_name << "\n";
	
  return os;
}

/**
 * Administrator of a vector of basis functions. This class build the relationship
 * between those basis functions to those degree of freedom on a template element.
 * Then the complete information to calculate the value of a basis function are
 * provided until now.
 */
template <class value_type, int DIM, int TDIM=DIM>
  class BasisFunctionAdmin : public std::vector<BasisFunction<value_type,DIM,TDIM> >
  {
  public:
  typedef value_type value_t;
  enum { dim = DIM, tdim = TDIM };
  typedef BasisFunction<value_t,DIM,TDIM> basis_func_t;

  private:
  typedef BasisFunctionAdmin<value_type,DIM,TDIM> basis_func_admin_t;

  std::string library_path;
  TemplateDOF<TDIM> * df; /**< Related degree of freedom on template element. */
  public:
  BasisFunctionAdmin(); /**< Default constructor. */
  BasisFunctionAdmin(const int&); /**< Initilized with certain length. */
  BasisFunctionAdmin(const int&, TemplateDOF<TDIM>&); /**< Initilized with certain length and template DOF. */
  BasisFunctionAdmin(TemplateDOF<TDIM>&); /**< Initilized with template DOF. */
  BasisFunctionAdmin(const basis_func_admin_t&); /**< Copy constructor. */
  ~BasisFunctionAdmin(); /**< Destructor. */
  public:
  void reinit(TemplateDOF<TDIM>&); /**< Reinitialized with certain template DOF. */
  basis_func_admin_t& operator=(const basis_func_admin_t&); /**< Copy operator. */
  const TemplateDOF<TDIM>& dof() const; /**< Template degree of freedom. */
  TemplateDOF<TDIM>& dof(); /**< Template degree of freedom. */
  void readData(const std::string&); /**< Read in the information from a file. */
  void writeData(const std::string&) const; /**< Write out the information to a file. */
  public:
  template <class STREAM,class VT,int GDIM,int TGDIM> 
  friend STREAM& operator>>(STREAM&, BasisFunctionAdmin<VT,GDIM,TGDIM>&); /**< Stream input. */
  template <class STREAM,class VT,int GDIM,int TGDIM> 
  friend STREAM& operator<<(STREAM&, const BasisFunctionAdmin<VT,GDIM,TGDIM>&); /**< Stream output. */
  };

template <class STREAM, class value_type, int DIM, int TDIM>
  STREAM& operator>>(STREAM& is, BasisFunctionAdmin<value_type,DIM,TDIM>& b)
{
  unsigned int i, j, k, l;
	
  is >> i;
  if (i != b.df->n_dof) {
    std::cerr << "number of basis functions: " << i
	      << "\n is not equal to"
	      << "\nnumber of dofs: " << b.df->n_dof
	      << std::endl;
    abort();
  }
  b.resize(i);
  std::vector<std::vector<int> > count;
  j = b.df->n_geometry_dof.size();
  count.resize(j);
  for (k = 0;k < j;k ++)
    count[k].resize(b.df->n_geometry_dof[k].size(), 0);
  for (j = 0;j < i;j ++) {
    is >> k >> l;
    b[b.df->geometry_dof[k][l][count[k][l]]].library_path = b.library_path;
    is >> b[b.df->geometry_dof[k][l][count[k][l] ++]];
  }
  return is;
}

template <class STREAM, class value_type, int DIM, int TDIM>
  STREAM& operator<<(STREAM& os, const BasisFunctionAdmin<value_type,DIM,TDIM>& b)
{
  int i, j;
	
  i = b.size();
  os << i << "\n";
  for (j = 0;j < i;j ++) {
    os << "\t" << b.df->dof_index[j].dimension
       << b.df->dof_index[j].geometry_index << "\n";
    os << b[j] << "\n";
  }
  return os;
}


template <int DIM>
class UnitOutNormal
{
 public:
  enum { dim = DIM };
  typedef Point<DIM> point_t;
 private:
  void * handle; /**< Handle of the object to open the shared library. */
  std::string library_path;
  std::string library_name; /**< Name of the shared library. */
  std::string function_name;
  void (*function)(const double *, const double **, int, double *);
 public:
  UnitOutNormal();
  UnitOutNormal(const UnitOutNormal<DIM>&);
  ~UnitOutNormal();
 public:
  UnitOutNormal<DIM>& operator=(const UnitOutNormal<DIM>&);
  void loadFunction();
  void unloadFunction();
  std::vector<double> value(const point_t&, 
                            const std::vector<point_t >&, const int&) const;
  std::vector<std::vector<double> > value(const std::vector<point_t >&, 
                                          const std::vector<point_t >&, const int&) const;
  std::vector<double> value(const point_t&, 
                            const double **, const int&) const;
  std::vector<std::vector<double> > value(const std::vector<point_t >&, 
                                          const double **, const int&) const;
  void readData(const std::string&); /**< Read in the information from a file. */
  void writeData(const std::string&) const; /**< Write out the information to a file. */

  DeclException1(ExcFileOpen, char *,
                 << "Can't open library "
                 << arg1);
  DeclException2(ExcLoadFunction, char *, char *,
                 << "Can't load function "
                 << arg1
                 << " from library "
                 << arg2);	
 public:
  template <class STREAM,int GDIM> 
    friend STREAM& operator>>(STREAM&, UnitOutNormal<GDIM>&); /**< Stream input. */
  template <class STREAM,int GDIM> 
    friend STREAM& operator<<(STREAM&, const UnitOutNormal<GDIM>&); /**< Stream output. */
};

template <class STREAM, int DIM>
  STREAM& operator>>(STREAM& is, UnitOutNormal<DIM>& c)
{
  is >> c.library_name
     >> c.function_name;
	
  c.loadFunction();
  return is;
}

template <class STREAM, int DIM>
  STREAM& operator<<(STREAM& os, const UnitOutNormal<DIM>& c)
{
  os << c.library_name << "\n\t"
     << c.function_name << "\n";
	
  return os;
}

/**
 * Template element is the package of those information to implement all finite element
 * operations on a template element geometry. It's only a simple package of information
 * and provides those accesses to those information.
 */
template <class value_type, int DIM, int TDIM=DIM>
  class TemplateElement
  {
  public:
  typedef value_type value_t;
  enum { dim = DIM, tdim = TDIM };
  typedef Point<DIM> point_t;
  typedef Point<TDIM> ref_point_t;
  typedef TemplateGeometry<TDIM> geometry_t;
  typedef CoordTransform<TDIM,DIM> coord_trans_t;
  typedef BasisFunction<value_t,DIM,TDIM> basis_func_t;
  typedef TemplateElement<value_t,DIM,TDIM> template_t;
  typedef QuadratureInfoAdmin<TDIM> quad_info_t;
  typedef UnitOutNormal<DIM> unit_normal_t;

  private:
  typedef BasisFunctionAdmin<value_type,DIM,TDIM> basis_func_admin_t;

  geometry_t * geo; /**< Template element geometry. */
  TemplateDOF<TDIM> * df; /**< Degree of freedom distribution on the template element geometry. */
  coord_trans_t * ct; /**< Coordinate transformation. */
  basis_func_admin_t * bf; /**< Basis functions. */
  unit_normal_t * uon; /**< Unit Our Normal. */
  public:
  TemplateElement() {} /**< Default constructor. */
  TemplateElement(geometry_t&, 
                  TemplateDOF<TDIM>&, 
                  coord_trans_t&, 
                  basis_func_admin_t&, 
                  unit_normal_t& = *((unit_normal_t *)NULL)); 
  TemplateElement(const template_t&); /**< Copy constructor. */
  ~TemplateElement(); /**< Destructor. */
  public:
  template_t& operator=(const template_t&); /**< Copy operator. */
  void reinit(geometry_t&, 
              TemplateDOF<TDIM>&,
              coord_trans_t&,
              basis_func_admin_t&,
              unit_normal_t& = *((unit_normal_t *)NULL)); /**< Reinitialization. */
  /** Template element geometry. */
  const geometry_t& geometry() const {return *geo;}
  /** Template element geometry. */
  geometry_t& geometry() {return *geo;}
  /** DOF distribution on emplate element geometry. */
  const TemplateDOF<TDIM>& dof() const {return *df;}
  /** DOF distribution on emplate element geometry. */
  TemplateDOF<TDIM>& dof() {return *df;}
  /** Coordinate transformation. */
  const coord_trans_t& coordTransform() const {return *ct;}
  /** Coordinate transformation. */
  coord_trans_t& coordTransform() {return *ct;}
  /** Basis functions. */
  const basis_func_admin_t& basisFunction() const {return *bf;}
  /** Basis functions. */
  basis_func_admin_t& basisFunction() {return *bf;}
  /** The \p i-th basis function. */
  const basis_func_t& basisFunction(const int& i) const {return (*bf)[i];}
  /** The \p i-th basis function. */
  basis_func_t& basisFunction(const int& i) {return (*bf)[i];}

  /** Vertex array of the template element geometry. */
  const std::vector<ref_point_t >& vertexArray() const;
  /** Quadrature information on the template element geometry. */
  const quad_info_t& quadratureInfo() const {return geo->quadratureInfo();}
  /** Quadrature information on the template element geometry. */
  quad_info_t& quadratureInfo() {return geo->quadratureInfo();}
  const unit_normal_t& unitOutNormal() const {return *uon;}
  unit_normal_t& unitOutNormal() {return *uon;}
  /** Quadrature information with algebraic accuracy \p i on the template element geometry. */
  const QuadratureInfo<TDIM>& findQuadratureInfo(const int& i) const {return geo->findQuadratureInfo(i);}
  /** Volume of the template element geometry. */
  double volume() const {return geo->volume();}
  int n_dof() const {return df->n_dof;}
  };


///////////////////////////////////////////////////////////////////////////


/**
 * The following struct are used for optimization of the finite element code.
 * If the mesh on which the finite element space built on is not changed, some
 * infomation about the geometry information and the basis function on the
 * quadrature points will not vary, too. Thus these information can be stored
 * there to make the code more efficient. Note to use this facility only when
 * you have enough memory.
 */
template <int DIM>
struct GeometryAdditionalData 
{
public:
  int n_quadrature_point;
  double volume;
  Point<DIM> bc;
  std::vector<double> Jxw;
  std::vector<Point<DIM> > q_point;
};

AFEPACK_CLOSE_NAMESPACE

#endif //_TemplateElement_h_

//
// end of file
///////////////////////////////////////////////////////////////////////////
