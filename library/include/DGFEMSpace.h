///////////////////////////////////////////////////////////////////////////////
// DGFEMSpace.h : by R.Lie, July 7, 2003
//

#ifndef _DGFEMSpace_h_
#define _DGFEMSpace_h_

#include <lac/full_matrix.h>

#include "Miscellaneous.h"
#include "Geometry.h"
#include "TemplateElement.h"
#include "FEMSpace.h"

AFEPACK_OPEN_NAMESPACE

template <int TDIM, int DIM> class TemplateDGElement;
template <class value_type, int DIM, int DOW, int TDIM, int TDIM1> class DGElement;
template <class value_type, int DIM, int DOW, int TDIM, int TDIM1> class DGFEMSpace;

template <class value_type, int DIM, int DOW, int TDIM, int TDIM1> std::vector<double> unitOutNormal(const Point<DIM>&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);
template <class value_type, int DIM, int DOW, int TDIM, int TDIM1> std::vector<std::vector<double> > unitOutNormal(const std::vector<Point<DIM> >&, const Element<value_type,DIM,DOW,TDIM>&, const DGElement<value_type,DIM,DOW,TDIM,TDIM1>&);

///////////////////////////////////////////////////////////////////////////////

template <int TDIM, int DIM=TDIM+1>
  class TemplateDGElement
  {
  private:
  TemplateGeometry<TDIM> * geo;
  CoordTransform<TDIM,DIM> * ct;
  public:
  TemplateDGElement() {}
  TemplateDGElement(TemplateGeometry<TDIM>& _geo,
                    CoordTransform<TDIM,DIM>& _ct) :
  geo(&_geo), ct(&_ct) {};
  TemplateDGElement(const TemplateDGElement<TDIM,DIM>& other) :
  geo(other.geo), ct(other.ct) {};
  ~TemplateDGElement() {};
  public:
  void reinit(TemplateGeometry<TDIM>& _geo,
              CoordTransform<TDIM,DIM>& _ct) {geo = &_geo; ct = &_ct;};
  const TemplateGeometry<TDIM>& geometry() const {return *geo;};
  TemplateGeometry<TDIM>& geometry() {return *geo;};
  const CoordTransform<TDIM,DIM>& coordTransform() const {return *ct;};
  CoordTransform<TDIM,DIM>& coordTransform() {return *ct;};
  const std::vector<Point<TDIM> >& vertexArray() const {return geo->vertexArray();};
  const QuadratureInfoAdmin<TDIM>& quadratureInfo() const {return geo->quadratureInfo();};
  QuadratureInfoAdmin<TDIM>& quadratureInfo() {return geo->quadratureInfo();};
  const QuadratureInfo<TDIM>& findQuadratureInfo(const int& i) const {return geo->findQuadratureInfo(i);};
  double volume() const {return geo->volume();};
  };

template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, int TDIM1=DIM-1>
  class DGElement
  {
  public:
  enum { dim = DIM, dow = DOW, tdim = TDIM, tdim1 = TDIM1 };
  typedef value_type value_t;
  typedef DGFEMSpace<value_t,DIM,DOW,TDIM,TDIM1> fe_space_t;
  typedef Element<value_t,DIM,DOW,TDIM> element_t;
  typedef TemplateDGElement<TDIM1,DOW> dg_template_t;
  typedef DGElement<value_t,DIM,DOW,TDIM,TDIM1> dg_element_t;
  typedef typename Mesh<DIM,DOW>::bmark_t bmark_t;
  typedef CoordTransform<TDIM1,DOW> coord_trans_t;

  private:
  fe_space_t * sp; /**< Pointer to finite element space it belongs to. */
  int geometry_index; /**< Index of real geometry in the mesh the finite element space on. */
  int template_element_index; /**< Index of template element. */
  std::vector<std::vector<int> > geo_img; /**< Geometry image built according template geometry. */
  element_t * neigh[2]; /**< the two neighbours of this face. */
  int bnd_idx[2]; /**< the index of this face as the boundary of the neighbour. */

  public:
  DGElement() {};  /**< Default contructor. */
  DGElement(fe_space_t&);
  DGElement(const dg_element_t&);  /**< Copy constructor. */
  ~DGElement(); /**< Destructor. */
  public:
  void reinit(fe_space_t&, const int&, const int&); /**< Reinitialization. */
  void reinit(const int&, const int&); /**< Reinitialization. */
  dg_element_t& operator=(const dg_element_t&); /**< Copy operator. */
  const fe_space_t& femSpace() const; /**< Finite element space. */
  fe_space_t& femSpace(); /**< Finite element space. */
  const int& index() const; /**< Geometry index. */
  int& index(); /**< Geometry index. */
  const int& boundaryIndex(const u_int&) const;
  int& boundaryIndex(const u_int&);
  GeometryBM& geometry() const; /**< Geometry. */
  bmark_t boundaryMark() const; /**< Boundary mark of its geometry. */
  void geometry(const GeometryBM&); /**< Set geometry. */
  void geometry(const int&); /**< Set geometry according geometry index. */
  TemplateDGElement<TDIM1,DOW>& templateElement() const; /**< Template element. */
  void templateElement(const int&); /**< Set template element according index. */
  TemplateGeometry<TDIM1>& templateGeometry() const; /**< Template element geometry. */
  const std::vector<std::vector<int> >& geometryImage() const; /**< Geometry image. */
  std::vector<std::vector<int> >& geometryImage(); /**< Geometry image. */
  void buildGeometryImage(); /**< Build geometry image. */
  void buildVertexArray(std::vector<Point<DOW> >&) const; /**< Build vertex array. */
  const double ** buildVertexArray() const; /**< Build vertex array. */
	
  Point<DOW> local_to_global(const Point<TDIM1>&) const; /**< Map a point from template element to this element. */
  Point<TDIM1> global_to_local(const Point<DOW>&) const; /**< Map a point from this element to template element. */
  double local_to_global_jacobian(const Point<TDIM1>&) const; /**< Jacobian determinant at a point on template element. */
  double global_to_local_jacobian(const Point<DOW>&) const; /**< Jacobian deterninant at a point on this element. */
  std::vector<Point<DOW> > local_to_global(const std::vector<Point<TDIM1> >&) const; /**< Map points from template element to this element. */
  std::vector<Point<TDIM1> > global_to_local(const std::vector<Point<DOW> >&) const; /**< Map points from this element to template element. */
  std::vector<double> local_to_global_jacobian(const std::vector<Point<TDIM1> >&) const; /** Jacobian determinant at points on template element. */
  std::vector<double> global_to_local_jacobian(const std::vector<Point<DOW> >&) const; /**< Jacobian determinant at points on real element. */

  const element_t * p_neighbourElement(const int& i) const;
  element_t * p_neighbourElement(const int& i);
  const element_t& neighbourElement(const int& i) const;
  element_t& neighbourElement(const int& i);

  const QuadratureInfo<TDIM1>& findQuadratureInfo(const int&) const; /**< Quadrature information with certain algebraic accuracy. */
  std::vector<double> unitNormal() const;

  public:
  DeclException1(ExcMeshData, std::string, << "Mesh data uncompatible: " << arg1);
  public:
  friend std::vector<double> unitOutNormal <>(const Point<DIM>&, 
                                              const element_t&, 
                                              const dg_element_t&);
  friend std::vector<std::vector<double> > unitOutNormal <>(const std::vector<Point<DIM> >&, 
                                                            const element_t&, 
                                                            const dg_element_t&);
  friend class DGFEMSpace<value_type, DIM, DOW, TDIM, TDIM1>;
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

template <class value_type, int DIM, int DOW=DIM, int TDIM=DIM, int TDIM1=DIM-1>
  class DGFEMSpace : public FEMSpace<value_type, DIM, DOW, TDIM>
  {
  public:
  enum { dim = DIM, dow = DOW, tdim = TDIM, tdim1 = TDIM1 };
  typedef value_type value_t;
  typedef DGFEMSpace<value_t,DIM,DOW,TDIM,TDIM1> fe_space_t;
  typedef DGElement<value_t,DIM,DOW,TDIM,TDIM1> dg_element_t;
  typedef FEMSpace<value_t, DIM, DOW, TDIM> base_t;
  typedef TemplateDGElement<TDIM1,DOW> dg_template_t;
  typedef typename base_t::template_t template_t;
  typedef typename base_t::element_t element_t;

  typedef typename std::vector<dg_element_t>::iterator DGElementIterator;
  typedef typename std::vector<dg_element_t>::const_iterator ConstDGElementIterator;

  private:
  std::vector<dg_template_t> * tmp_dgele;
  std::vector<dg_element_t> dgele;
  public:
  DGFEMSpace() {}
  DGFEMSpace(Mesh<DIM,DOW>&,
	     std::vector<template_t>&,
	     std::vector<dg_template_t>&);
  DGFEMSpace(const fe_space_t&);
  virtual ~DGFEMSpace();
  public:
  void reinit(Mesh<DIM,DOW>&, 
              std::vector<template_t>& = *((std::vector<template_t> *)NULL),
              std::vector<dg_template_t>& = *((std::vector<dg_template_t> *)NULL));
  fe_space_t& operator=(const fe_space_t&);
  const std::vector<dg_template_t>& templateDGElement() const {return *tmp_dgele;};
  std::vector<dg_template_t>& templateDGElement() {return *tmp_dgele;}
  const dg_template_t& templateDGElement(const int& i) const {return (*tmp_dgele)[i];};
  dg_template_t& templateDGElement(const int& i) {return (*tmp_dgele)[i];};
  int n_DGElement() const {return dgele.size();};
  const std::vector<dg_element_t>& dgElement() const {return dgele;};
  std::vector<dg_element_t>& dgElement() {return dgele;};
  const dg_element_t& dgElement(const int& i) const {return dgele[i];};
  dg_element_t& dgElement(const int& i) {return dgele[i];};
  virtual void buildDGElement();
  DGElementIterator beginDGElement() {return dgElement().begin();};
  DGElementIterator endDGElement() {return dgElement().end();};
  ConstDGElementIterator beginDGElement() const {return dgElement().begin();};
  ConstDGElementIterator endDGElement() const {return dgElement().end();};

#ifdef __SERIALIZATION__
  template <class T> T *
  new_property(const dg_element_t& ele,
               const property_id_t<T>& pid) const {
    typedef RegularMesh<DIM,DOW> reg_mesh_t;
    const reg_mesh_t& mesh = dynamic_cast<const reg_mesh_t&>(this->mesh());
    return mesh.template new_property<T,DIM-1>(ele.geometry().index(), pid);
  }
  template <class T> T *
  get_property(const dg_element_t& ele,
               const property_id_t<T>& pid) const {
    typedef RegularMesh<DIM,DOW> reg_mesh_t;
    const reg_mesh_t& mesh = dynamic_cast<const reg_mesh_t&>(this->mesh());
    return mesh.template get_property<T,DIM-1>(ele.geometry().index(), pid);
  }
  template <class T> void
  free_property(const dg_element_t& ele,
                const property_id_t<T>& pid) const {
    typedef RegularMesh<DIM,DOW> reg_mesh_t;
    const reg_mesh_t& mesh = dynamic_cast<const reg_mesh_t&>(this->mesh());
    mesh.template free_property<T,DIM-1>(ele.geometry().index(), pid);
  }

  using base_t::new_property;
  using base_t::get_property;
  using base_t::free_property;
#endif // __SERIALIZATION__
  public:
  friend class DGElement<value_t,DIM,DOW,TDIM,TDIM1>;
  };


///////////////////////////////////////////////////////////////////////////////

/**
 * The following struct are used for optimization of the finite element code.
 * If the mesh on which the finite element space built on is not changed, some
 * infomation about the geometry information and the basis function on the
 * quadrature points will not vary, too. Thus these information can be stored
 * there to make the code more efficient. Note to use this facility only when
 * you have enough memory.
 */
template <typename value_type, int DIM>
  struct DGElementAdditionalData : public GeometryAdditionalData<DIM>
{
 public:
  std::vector<std::vector<value_type> > basis_value0;
  std::vector<std::vector<value_type> > basis_value1;
  std::vector<std::vector<double> > unit_normal;
};

AFEPACK_CLOSE_NAMESPACE

#endif // _DGFEMSpace_h_

//
// end of file
///////////////////////////////////////////////////////////////////////////////
