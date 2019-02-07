/**
 * @file   FEMSpace.h
 * @author Ruo Li
 * @date   Thu Feb 17 10:58:17 2005
 * 
 * @brief  declaration of class Element, class FEMSpace, class FEMFunction
 * 
 * 
 */

#ifndef _FEMSpace_h_
#define _FEMSpace_h_

#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>

#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>

#include "Geometry.h"
#include "HGeometry.h"
#include "TemplateElement.h"
#ifdef MULTITHREAD
#include "Thread.h"
#endif

AFEPACK_OPEN_NAMESPACE

template <class value_type, int DIM, int DOW, int TDIM> class Element;
template <class value_type, int DIM, int DOW, int TDIM> class FEMSpace;
template <class value_type, int DIM, int DOW, int TDIM, typename Number> class FEMFunction;

///////////////////////////////////////////////////////////////////////////

/**
 * Element in finite element space. An element is contructed by mapping a template
 * element to a real element geometry, then all those information about degree of
 * freedom and basis function get its image on the real element. Then many operations
 * needed in finite element calculation can be implemented on class \p Element.
 */
template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM>
  class Element
  {
  public:
  enum { dim = DIM, dow = DOW, tdim = TDIM };
  typedef value_type value_t;
  typedef FEMSpace<value_t,DIM,DOW,TDIM> fe_space_t;
  typedef Element<value_t,DIM,DOW,TDIM> element_t;
  typedef TemplateElement<value_t,DOW,TDIM> template_t;

  private:
  fe_space_t * sp; /**< Pointer to finite element space it belongs to. */
  int geometry_index; /**< Index of real geometry in the mesh the finite element space on. */
  int template_element_index; /**< Index of template element. */
  std::vector<int> dof_index; /**< Information of degree of freedom on this element. */
  std::vector<std::vector<int> > geo_img; /**< Geometry image built according template geometry. */
  public:
  Element() {}
  explicit Element(fe_space_t&); /**< Default contructor. */
  Element(const element_t&);  /**< Copy constructor. */
  ~Element(); /**< Destructor. */
  public:
  void reinit(fe_space_t&, int, int); /**< Reinitialization. */
  void reinit(int, int); /**< Reinitialization. */
  void reinit(fe_space_t&, int, int, const std::vector<int>&); /**< Reinitialization. */
  element_t& operator=(const element_t&); /**< Copy operator. */
  const fe_space_t& femSpace() const; /**< Finite element space. */
  fe_space_t& femSpace(); /**< Finite element space. */
  Mesh<DIM,DOW>& mesh() const; /**< Mesh the finite element space on. */
  int index() const; /**< Geometry index. */
  int& index(); /**< Geometry index. */
  GeometryBM& geometry() const; /**< Geometry. */
  void geometry(const GeometryBM&); /**< Set geometry. */
  void geometry(int); /**< Set geometry according geometry index. */
  int templateElementIndex() const; /**< Get template element index. */
  template_t& templateElement() const; /**< Template element. */
  void templateElement(int); /**< Set template element according index. */
  TemplateGeometry<TDIM>& templateGeometry() const; /**< Template element geometry. */
  const std::vector<int>& dof() const; /**< Degree of freedom. */
  std::vector<int>& dof(); /**< Degree of freedom. */
  unsigned int n_dof() const; /**< Number of degree of freedom. */
  const std::vector<std::vector<int> >& geometryImage() const; /**< Geometry image. */
  std::vector<std::vector<int> >& geometryImage(); /**< Geometry image. */
  void buildGeometryImage(); /**< Build geometry image. */
  void lazyBuildGeometryImage(); /**< Build geometry image. */
  void buildVertexArray(std::vector<Point<DOW> >&) const; /**< Build vertex array. */
  const double ** buildVertexArray() const; /**< Build vertex array. */
	
  const typename BasisFunction<value_t,DIM,TDIM>::Identity& basis_function_identity(int i) const; /**< Basis function identity of certain basis function. */
  value_t basis_function_value(int i, const Point<DOW>&) const; /**< Value of certain basis function at a point. */
  std::vector<value_t> basis_function_value(const Point<DOW>&) const; /**< Value of all basis functions at a point. */
  std::vector<value_t> basis_function_value(int i, const std::vector<Point<DOW> >&) const; /**< Value of certain basis function at points. */
  std::vector<std::vector<value_t> > basis_function_value(const std::vector<Point<DOW> >&) const; /**< Value of all basis functions at points. */
  std::vector<value_t> basis_function_gradient(int i, const Point<DOW>&) const; /**< Gradient of certain basis function at a point. */
  std::vector<std::vector<value_t> > basis_function_gradient(const Point<DOW>&) const; /**< Gradient of all basis functions at a point. */
  std::vector<std::vector<value_t> > basis_function_gradient(int i, const std::vector<Point<DOW> >&) const; /**< Gradient of certain basis function at points. */
  std::vector<std::vector<std::vector<value_t> > > basis_function_gradient(const std::vector<Point<DOW> >&) const; /**< Gradient of all basis functions at points. */

#ifdef QUADRATIC_ELEMENT_SUPPORT
  std::vector<std::vector<value_t> > 
  basis_function_hesse(int i, 
                       const Point<DOW>&) const; 
  std::vector<std::vector<std::vector<value_t> > > 
  basis_function_hesse(const Point<DOW>&) const;
  std::vector<std::vector<std::vector<value_t> > >
  basis_function_hesse(int i, 
                       const std::vector<Point<DOW> >&) const;
  std::vector<std::vector<std::vector<std::vector<value_t> > > >
  basis_function_gradient(const std::vector<Point<DOW> >&) const;
#endif // QUADRATIC_ELEMENT_SUPPORT

  Point<DOW> local_to_global(const Point<TDIM>&) const; /**< Map a point from template element to this element. */
  Point<TDIM> global_to_local(const Point<DOW>&) const; /**< Map a point from this element to template element. */
  double local_to_global_jacobian(const Point<TDIM>&) const; /**< Jacobian determinant at a point on template element. */
  double global_to_local_jacobian(const Point<DOW>&) const; /**< Jacobian deterninant at a point on this element. */
  std::vector<Point<DOW> > local_to_global(const std::vector<Point<TDIM> >&) const; /**< Map points from template element to this element. */
  std::vector<Point<TDIM> > global_to_local(const std::vector<Point<DOW> >&) const; /**< Map points from this element to template element. */
  std::vector<double> local_to_global_jacobian(const std::vector<Point<TDIM> >&) const; /** Jacobian determinant at points on template element. */
  std::vector<double> global_to_local_jacobian(const std::vector<Point<DOW> >&) const; /**< Jacobian determinant at points on real element. */

  std::vector<double> unitOutNormal(const Point<DOW>&, int) const;
  std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DOW> >&, int) const;
	
  const QuadratureInfo<TDIM>& findQuadratureInfo(int) const; /**< Quadrature information with certain algebraic accuracy. */

  public:
  DeclException1(ExcMeshData, std::string, << "Mesh data uncompatible: " << arg1);
#ifdef __SERIALIZATION__ 
  template <class T> T *
  new_property(const property_id_t<T>& pid) const {
    return femSpace().new_property(*this, pid);
  }
  template <class T> T *
  get_property(const property_id_t<T>& pid) const {
    return femSpace().get_property(*this, pid);
  }
  template <class T> void
  free_property(const property_id_t<T>& pid) const {
    femSpace().free_property(*this, pid);
  }
#endif // __SERIALIZATION__
  };

/**
 * Finite element space is the main class used to implement the global manipulation of
 * degree of freedom. It contains a template element array, which is the templates for
 * all the element, and an element array. Other data are used to construct the information
 * of global degree of freedom. The class is designed open-minded for the user to
 * implement different finite element classes freely by re-implement the member
 * function \p buildElement.
 */
template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM>
  class FEMSpace
  {
  public:
  enum { dim = DIM, dow = DOW, tdim = TDIM };

  typedef value_type value_t;
  typedef Element<value_t,DIM,DOW,TDIM> element_t;
  typedef FEMSpace<value_t,DIM,DOW,TDIM> fe_space_t;
  typedef DOFInfo<value_t, DIM, DOW,TDIM> dof_info_t;
  typedef typename element_t::template_t template_t;
  typedef Mesh<DIM,DOW> mesh_t;
  typedef typename mesh_t::bmark_t bmark_t;
  typedef typename std::vector<element_t>::iterator ElementIterator;
  typedef typename std::vector<element_t>::const_iterator ConstElementIterator;

  private:
  mesh_t * msh; /**< Mesh the finite element space on. */
  std::vector<template_t> * tmp_ele; /**< Template element array. */
  std::vector<element_t> ele; /**< Element array. */
  DegreeOfFreedom df; /**< Degree of freedom. */
  std::vector<dof_info_t> df_in; /**< Information of degree of freedoms. */

  u_int effi_flag;		/**< the flag is used for efficiency */

  public:
  FEMSpace() : effi_flag(0xFFFF) {};  /**< Default constructor. */
  explicit FEMSpace(mesh_t&, std::vector<template_t>&); 
  FEMSpace(const fe_space_t&); /**< Copy constructor. */
  virtual ~FEMSpace(); /**< Destructor. */
  public:
  fe_space_t& operator=(const fe_space_t&); /**< Copy operator. */
  void reinit(mesh_t& = *((mesh_t *)NULL),
              std::vector<template_t>& 
              = *((std::vector<template_t> *)NULL)); /**< Reinitialization. */

  const mesh_t& mesh() const {return *msh;}
  mesh_t& mesh() {return *msh;}

  const std::vector<template_t>& templateElement() const {return *tmp_ele;}
  std::vector<template_t>& templateElement() {return *tmp_ele;}
  const template_t& templateElement(int i) const {return (*tmp_ele)[i];}
  template_t& templateElement(int i) {return (*tmp_ele)[i];}

  int n_element() const {return ele.size();};
  const std::vector<element_t>& element() const {return ele;}
  std::vector<element_t>& element() {return ele;}
  const element_t& element(int i) const {return ele[i];}
  element_t& element(int i) {return ele[i];}

  const DegreeOfFreedom& dof() const {return df;}
  DegreeOfFreedom& dof() {return df;}

  const std::vector<dof_info_t>& dofInfo() const {return df_in;}
  std::vector<dof_info_t>& dofInfo() {return df_in;}
  const dof_info_t& dofInfo(int i) const {return df_in[i];}
  dof_info_t& dofInfo(int i) {return df_in[i];}

  unsigned int n_dof() const {return dof().n_dof;}
  const std::vector<int>& geometryDof(int i, int j) const {return dof().geometry_dof[i][j];}
  const std::vector<int>& geometryDof(int i, const Geometry& g) const {return dof().geometry_dof[i][g.index()];}

  const std::vector<DOFIndex>& dofIndex() const {return dof().dof_index;}
  const DOFIndex& dofIndex(int i) const {return dof().dof_index[i];}

  const bmark_t& dofBoundaryMark(int i) const {return dofInfo(i).boundary_mark;}
  bmark_t& dofBoundaryMark(int i) {return dofInfo(i).boundary_mark;}

  virtual void buildElement(bool lazy = false); /**< Build the element in the finite element space by appointing the
                                                   template element index. Need to be re-implemented in derived class. */
  virtual void buildDofBoundaryMark(); /**< Build the boundary marker of degree of
                                          freedom according the boundary marker of the geometry the degree of freedom on. */
  void buildDof(); /**< Build degree of freedom for the finite element space. */
  void updateDofInterpPoint(); /**< Update interpolate point of DOFs when the structure of the mesh is reserved while the location of the mesh points are changed. Designed for moving mesh method. */
  /** Begin element iterator. */
  ElementIterator beginElement() {return element().begin();};
  ConstElementIterator beginElement() const {return element().begin();};
  /** End element iterator. */
  ElementIterator endElement() {return element().end();};
  ConstElementIterator endElement() const {return element().end();};
	
  DeclException1(ExcDOFData, std::string, << "DOF data uncompatible: " << arg1);

  const unsigned int& efficiencyFlag() const {return effi_flag;};
  unsigned int& efficiencyFlag() {return effi_flag;}

  public:
#ifdef __SERIALIZATION__ 
  template <class T> T *
  new_property(const element_t& e,
               const property_id_t<T>& pid) const {
    typedef RegularMesh<DIM,DOW> reg_mesh_t;
    const reg_mesh_t * mesh = dynamic_cast<const reg_mesh_t *>(msh);
    return mesh->template new_property<T,DIM>(e.geometry().index(), pid);
  }
  template <class T> T *
  get_property(const element_t& e,
               const property_id_t<T>& pid) const {
    typedef RegularMesh<DIM,DOW> reg_mesh_t;
    const reg_mesh_t * mesh = dynamic_cast<const reg_mesh_t *>(msh);
    return mesh->template get_property<T,DIM>(e.geometry().index(), pid);
  }
  template <class T> void
  free_property(const element_t& e,
               const property_id_t<T>& pid) const {
    typedef RegularMesh<DIM,DOW> reg_mesh_t;
    const reg_mesh_t * mesh = dynamic_cast<const reg_mesh_t *>(msh);
    mesh->template free_property<T,DIM>(e.geometry().index(), pid);
  }
#endif // __SERIALIZATION__

  friend class Element<value_t,DIM,DOW,TDIM>;

#ifdef MULTITHREAD
  void threadBuildElement(bool lazy, int n_thread = 1, int rank = 0);
  void threadBuildDof0(std::vector<std::vector<bool> >&, 
                       pthread_mutex_t&, int n_thread = 1, int rank = 0);
  void threadBuildDof1(std::vector<std::vector<bool> >&, 
                       pthread_mutex_t&, int n_thread = 1, int rank = 0);
#endif
  };

/**
 * Finite element function is a vector associated with a finite element space. The
 * entry of the vector is related with the corresponding degree of freedom of the
 * finite element space.
 */
template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, typename Number=double>
  class FEMFunction : public Vector<Number>
  {
  public:
  enum { dim = DIM, dow = DOW, tdim = TDIM };

  typedef value_type value_t;
  typedef Number number_t;
  typedef Element<value_t,DIM,DOW,TDIM> element_t;
  typedef FEMSpace<value_t,DIM,DOW,TDIM> fe_space_t;
  typedef FEMFunction<value_t,DIM,DOW,TDIM,number_t> fe_func_t;

  private:
  fe_space_t * sp; /**< Finite element space. */
  public:
  FEMFunction() {}
  FEMFunction(fe_space_t&); /**< Default constructor. */
  virtual ~FEMFunction(); /**< Destructor. */
  public:
  fe_space_t& femSpace(); /**< Finite element space, obsolete */
  const fe_space_t& femSpace() const { return *sp; };
  void reinit(fe_space_t&, bool = false); /**< Reinitialization. */
	
  value_type value(const Point<DOW>&, const element_t&) const; /**< Value of the function at a point in certain element. */
  std::vector<value_type> gradient(const Point<DOW>&, const element_t&) const; /**< Value of the function at points in certain element. */
  std::vector<value_type> value(const std::vector<Point<DOW> >&, const element_t&) const; /**< Gradient of the function at a point in certain element. */
  std::vector<std::vector<value_type> > gradient(const std::vector<Point<DOW> >&, const element_t&) const; /**< Gradient of the function at points in certain element. */

  /**
   * for optimization of the code when the basis function value and gradient are assumed known,
   * the following functions will make it more efficient to calculate the value and gradient of the
   * finite element function.
   */
  value_type value(const std::vector<value_type>&, const element_t&) const; /**< Value of the function at a point in certain element assuming the basis function value is known. */
  std::vector<value_type> gradient(const std::vector<std::vector<value_type> >&, const element_t&) const; /**< Value of the function at points in certain element assuming the basis function gradient is known. */
  std::vector<value_type> value(const std::vector<std::vector<value_type> >&, const element_t&) const; /**< Gradient of the function at a point in certain element assuming the basis function value is known. */
  std::vector<std::vector<value_type> > gradient(const std::vector<std::vector<std::vector<value_type> > >&, const element_t&) const; /**< Gradient of the function at points in certain element assuming the basis function gradient is known. */

  public:
  void loadData(const std::string& filename); /**< Load data written by writeData. */
  void writeData(const std::string& filename); /**< Write the data of the function in internal data format.
                                                  This function will output the value of the vector only used to be rereaded in. */
  void writeEasyMeshData(const std::string& filename); /**< Write the data of the function in easymesh data format.
                                                          This function will output mesh into easymesh data format at first and then output the function value
                                                          on the nodes in the file named filename.dat. Only valid for 2 dimension and the mesh is turned into a 
                                                          RegularMesh object, for such otherwise object this function can't handle.*/
  /** 
   * Write the data of the function in Tecplot data format. This function 
   * will output the value of the function on those interpolating point of 
   * basis functions. 
   * 
   */
  void writeTecplotData(const std::string& filename);
  /**
   * Write the data of the function in IBM Open Data Explorer data format. This 
   * function will output the value of the function on nodes of the meshand will 
   * only be avaiable for certain case becase I'm not so familiar to this data format. 
   *
   * @param filename the file name to write to
   * @param flag indicate to use position dependent data or connection dependent data,
   *             flag=0: position dependent data;
   *             flag=1: connection dependent data;
   *             default 0
   */
  void writeOpenDXData(const std::string& filename,
                       int flag = 0) const;
  };



template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, typename Number=double>
  class LocalFEMFunction : public Vector<Number>
  {
  private:
  Element<value_type,DIM,DOW,TDIM> * ele;
  public:
  LocalFEMFunction() {}
  LocalFEMFunction(Element<value_type,DIM,DOW,TDIM> &);
  virtual ~LocalFEMFunction() {};
  public:
  Element<value_type,DIM,DOW,TDIM>& element();
  void reinit(Element<value_type,DIM,DOW,TDIM> &);
  value_type value(const Point<DOW>&) const;
  std::vector<value_type> gradient(const Point<DOW>&) const;
  std::vector<value_type> value(const std::vector<Point<DOW> >&) const;
  std::vector<std::vector<value_type> > gradient(const std::vector<Point<DOW> >&) const;	
  };



/**
 * Boundary condifiton information. This struct now only contains three static variables
 * about the type of different boundary condition. Unfortunately until now, only the Dirichlet
 * boundary condition is supported. Because the Neumann boundary condition is in fact in
 * a finite element space involved the element shell, I believe it can only be implemented
 * after the library can cope with hybrid element.
 */
struct BoundaryConditionInfo
{
public:
  typedef int	Type;
  static Type	DIRICHLET;
  static Type	NEUMANN;
  static Type	ROBIN;
};



template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, typename Number=double>
  class BoundaryCondition : public BoundaryConditionInfo
  {
  public:
  typedef typename FEMSpace<value_type,DIM,DOW,TDIM>::bmark_t bmark_t;
  typedef typename BoundaryConditionInfo::Type Type;
  private:
  Type					boundary_type;
  bmark_t			boundary_mark;
  public:
  BoundaryCondition() {};
  BoundaryCondition(const Type& t, const bmark_t& m) :
  boundary_type(t), boundary_mark(m) {};
  BoundaryCondition(const BoundaryCondition<value_type,DIM>& b) :
  boundary_type(b.boundaryType()), boundary_mark(b.boundaryMark()) {};
  virtual ~BoundaryCondition() {};
  public:
  const Type& boundaryType() const {return boundary_type;};
  Type& boundaryType() {return boundary_type;};
  const bmark_t& boundaryMark() const {return boundary_mark;};
  bmark_t& boundaryMark() {return boundary_mark;};
  void reinit(const Type& t, const bmark_t& m) {
    boundary_type = t;
    boundary_mark = m;
  }
  BoundaryCondition<value_type,DIM,DOW,TDIM>& operator=(const BoundaryCondition<value_type,DIM,DOW,TDIM>& b) {
    boundary_type = b.boundaryType();
    boundary_mark = b.boundaryMark();
    return *this;
  };
  public:
  virtual Number value(const Point<DOW>&) const {return 0;};
  };



/**
 * Boundary condition adopting a function as its value method. This class is used
 * to make a general C function be usable for the library.
 */
template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, typename Number=double>
  class BoundaryFunction : public BoundaryCondition<value_type,DIM,DOW,TDIM,Number>
  {
  public:
  typedef typename BoundaryCondition<value_type,DIM,DOW,TDIM,Number>::bmark_t bmark_t;
  typedef typename BoundaryCondition<value_type,DIM,DOW,TDIM,Number>::Type Type;
  private:
  bool is_newed;
  const Function<Number> * f;
  public:
  BoundaryFunction() : is_newed(false), f(NULL) {};
  BoundaryFunction(const Type& t,
                   const bmark_t& m,
                   const Function<Number>& fun) :
  BoundaryCondition<value_type,DIM>(t,m),
  is_newed(false),
  f(&fun)	{};
  BoundaryFunction(const Type& t,
                   const bmark_t& m, 
                   value_type (*fun)(const double *),
                   std::vector<value_type> (*grad)(const double *) = NULL) :
  BoundaryCondition<value_type,DIM,DOW,TDIM,Number>(t,m) 
  {
    is_newed = true;
    f = new FunctionFunction<Number>(fun, grad);
  };
  BoundaryFunction(const BoundaryFunction<value_type,DIM,DOW,TDIM,Number>& b) :
  BoundaryCondition<value_type,DIM,DOW,TDIM,Number>(b), 
  is_newed(false), 
  f(b.f) {};
  virtual ~BoundaryFunction() {
    if (is_newed && f != NULL) delete f;
  };
  public:
  void reinit(const Type& t,
              const bmark_t& m, 
              value_type (*fun)(const double *),
              std::vector<value_type> (*grad)(const double *) = NULL) {
    BoundaryCondition<value_type,DIM,DOW,TDIM,Number>::reinit(t, m);
    if (is_newed && f != NULL) delete f;
    is_newed = true;
    f = new FunctionFunction<Number>(fun, grad);
  }
  virtual Number value(const Point<DOW>& p) const {return f->value(p);};
  virtual std::vector<Number> gradient(const Point<DOW>& p) const {return f->gradient(p);};
  };



/**
 * Boundary condition. This class provides facilities to apply a Dirichlet boundary condition
 * on the sparse matrix and right hand side vector discretized from the finite element space.
 * It can include a list of boundary condition for different material boundary, and try retrieve
 * the correct boundary condition for them respectively.
 */
template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, typename Number=double>
  class BoundaryConditionAdmin : std::vector<const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *>
  {
  public:
  typedef typename FEMSpace<value_type,DIM,DOW,TDIM>::bmark_t bmark_t;
  private:
  std::vector<int>					index_map;
  const FEMSpace<value_type,DIM,DOW,TDIM> *		fem_space;
  public:
  BoundaryConditionAdmin(const FEMSpace<value_type,DIM,DOW,TDIM>& sp = 
                         *((const FEMSpace<value_type,DIM,DOW,TDIM> *)(NULL))) {fem_space = &sp;};
  ~BoundaryConditionAdmin() {};
  public:
  void reinit(const FEMSpace<value_type,DIM,DOW,TDIM>& sp) {fem_space = &sp; index_map.clear();};
  void setFemSpace(const FEMSpace<value_type,DIM,DOW,TDIM>& sp) {fem_space = &sp;};
  const FEMSpace<value_type,DIM,DOW,TDIM>& femSpace() const {return *fem_space;};
  void apply(SparseMatrix<double>& A, Vector<double>& u, Vector<double>& f, 
             bool preserve_symmetry = true); /**< Apply the boundary conditions on the given
                                                sparse matrix, right hand side vector and the undetermined variables. */
  void clearEntry(Vector<double>& f);
  void add(const BoundaryCondition<value_type,DIM,DOW,TDIM,Number>& b);

  bool isValid(const BoundaryCondition<value_type,DIM,DOW,TDIM,Number>& bc) const {
    return (&bc != NULL);
  };
  const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> * find(const bmark_t& bm) const {
    if ((unsigned int)bm >= index_map.size() || index_map[bm] == -1)
      return ((const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *)NULL);
    else
      return ((*this)[index_map[bm]]);
  };
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
template <typename value_type, int DIM>
  struct ElementAdditionalData : public GeometryAdditionalData<DIM>
{
 public:
  std::vector<std::vector<value_type> > basis_value;
  std::vector<std::vector<std::vector<value_type> > > basis_gradient;
};

namespace Migration {

  /** 
   * 获取树结构中几何体上的一个输入流，用来输出数据。
   * 
   * @param geo 树结构中的几何体
   * @param data_id 数据 ID
   * @param is 获得的输入流
   */
  template <class HGEO, class STREAM> void
    get_import_stream(HGEO& geo,
                      const data_id_t& data_id,
                      STREAM& is) {
    is.set_buffer(details::get_buffer(geo, data_id, false));
  }

  /** 
   * 获取网格中的几何体上的一个输出流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param data_id 数据 ID
   * @param dimension 几何体维数
   * @param geo_idx 几何体序号
   * @param os 获得的输出流
   */
  template <class MESH, class STREAM> void
    get_export_stream(MESH& mesh,
                      const data_id_t& data_id,
                      u_int dimension,
                      u_int geo_idx,
                      STREAM& os) {
    os.set_buffer(details::get_buffer(mesh, data_id, dimension, geo_idx, true));
  }
  
  /** 
   * 获取网格中几何体上的一个输入流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param data_id 数据 ID
   * @param dimension 几何体维数
   * @param geo_idx 几何体序号
   * @param is 获得的输入流
   */
  template <class MESH, class STREAM> void
    get_import_stream(MESH& mesh,
                      const data_id_t& data_id,
                      u_int dimension,
                      u_int geo_idx,
                      STREAM& is) {
    is.set_buffer(details::get_buffer(mesh, data_id, dimension, geo_idx, false));
  }

  /** 
   * 获取有限元空间上的一个自由度的输出流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param sp 有限元空间
   * @param data_id 数据 ID
   * @param dof 自由度编号
   * @param os 获得的输出流
   */
  template <class MESH, class SP, class STREAM> void
    get_dof_export_stream(MESH& mesh,
                          SP& sp,
                          const data_id_t& data_id,
                          int dof,
                          STREAM& os) {
    const DOFIndex& di = sp.dofIndex(dof);
    get_export_stream(mesh, data_id, di.dimension, di.geometry_index, os);
  }

  /** 
   * 获取有限元空间上的一个自由度的输入流，用来输出数据。
   * 
   * @param mesh 正则网格
   * @param sp 有限元空间
   * @param data_id 数据 ID
   * @param dof 自由度编号
   * @param is 获得的输入流
   */
  template <class MESH, class SP, class STREAM> void
    get_dof_import_stream(MESH& mesh,
                          SP& sp,
                          const data_id_t& data_id,
                          int dof,
                          STREAM& is) {
    const DOFIndex& di = sp.dofIndex(dof);
    get_import_stream(mesh, data_id, di.dimension, di.geometry_index, is);
  }

} // namespace Migration

AFEPACK_CLOSE_NAMESPACE

#endif //_FEMSpace_h_

/**
 * end of file
 * 
 */

