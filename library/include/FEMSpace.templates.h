/**
 * @file   FEMSpace.templates.h
 * @author Robert Lie
 * @date   Wed Jun 21 07:47:46 2003
 * 
 * @brief  
 * 
 * 
 */

#ifndef _FEMSpace_templates_h_
#define _FEMSpace_templates_h_

#include "FEMSpace.h"

AFEPACK_OPEN_NAMESPACE

#define TEMPLATE template <class value_type, int DIM, int DOW, int TDIM>
#define THIS Element<value_type, DIM, DOW, TDIM>

TEMPLATE
THIS::Element(fe_space_t& f) :
sp(&f)
{}

TEMPLATE
THIS::Element(const element_t& e) :
sp(e.sp),
  geometry_index(e.geometry_index),
  template_element_index(e.template_element_index),
  dof_index(e.dof_index)
{}

TEMPLATE
THIS::~Element()
{}

TEMPLATE
void THIS::reinit(fe_space_t& f, 
                  int g, 
                  int t)
{
  sp = &f;
  geometry_index = g;
  template_element_index = t;
}

TEMPLATE
void THIS::reinit(int g, int t)
{
  Assert (sp != NULL, ExcInternalError());
  geometry_index = g;
  template_element_index = t;
}

TEMPLATE
void THIS::reinit(fe_space_t& f, 
                  int g, 
                  int t, 
                  const std::vector<int>& d)
{
  sp = &f;
  geometry_index = g;
  template_element_index = t;
  dof_index = d;
}

TEMPLATE
typename THIS::element_t& THIS::operator=(const element_t& e)
{
  if (&e != NULL) {
    sp = e.sp;
    geometry_index = e.geometry_index;
    template_element_index = e.template_element_index;
    dof_index = e.dof_index;
  }
  return *this;
}

  TEMPLATE
  const typename THIS::fe_space_t& THIS::femSpace() const
  {
    return *sp;
  }

TEMPLATE
typename THIS::fe_space_t& THIS::femSpace()
{
  return *sp;
}

TEMPLATE
Mesh<DIM,DOW>& THIS::mesh() const
{
  return sp->mesh();
}

TEMPLATE
int THIS::index() const
{
  return geometry_index;
}

TEMPLATE
int& THIS::index()
{
  return geometry_index;
}

TEMPLATE
GeometryBM& THIS::geometry() const
{
  return sp->mesh().geometry(DIM, geometry_index);
}

TEMPLATE
void THIS::geometry(const GeometryBM& g)
{
  geometry_index = g.index();
}

TEMPLATE
void THIS::geometry(int i)
{
  geometry_index = i;
}

TEMPLATE
int THIS::templateElementIndex() const
{
  return template_element_index;
}

TEMPLATE
typename THIS::template_t& 
THIS::templateElement() const
{
  return sp->templateElement(template_element_index);
}

TEMPLATE
void THIS::templateElement(int i)
{
  template_element_index = i;
}

TEMPLATE
TemplateGeometry<TDIM>& THIS::templateGeometry() const
{
  return sp->templateElement(template_element_index).geometry();
}

TEMPLATE
const std::vector<int>& THIS::dof() const
{
  return dof_index;
}

TEMPLATE
std::vector<int>& THIS::dof()
{
  return dof_index;
}

TEMPLATE
unsigned int THIS::n_dof() const
{
  return dof_index.size();
}

TEMPLATE
const std::vector<std::vector<int> >& THIS::geometryImage() const
{
  return geo_img;
}

TEMPLATE
std::vector<std::vector<int> >& THIS::geometryImage()
{
  return geo_img;
}

TEMPLATE
void THIS::buildGeometryImage()
{
  Mesh<DIM,DOW>& m = mesh();
  GeometryBM& geo = this->geometry();
  template_t& t_el = templateElement();
  TemplateGeometry<TDIM>& t_geo = t_el.geometry();
	
  int i, j, k, l, i1, j1, k1, l1;
	
  geo_img.resize(TDIM+1);
  for (i = 0;i <= TDIM;i ++)
    geo_img[i].resize(t_geo.n_geometry(i), -1);
  geo_img[TDIM][0] = geo.index();
  geo_img[0] = geo.vertex();
  i1 = geo_img[0].size();
  for (i = TDIM-1;i > 0;i --) {
    // collect the geometry index in dimension \p{i}
    j = t_geo.n_geometry(i+1);
    std::vector<int> geo_index(t_geo.n_geometry(i), -1);
    for (j1 = 0, l = 0;j1 < j;j1 ++) {
      const GeometryBM& geo1 = m.geometry(i+1, geo_img[i+1][j1]);
      k = geo1.n_boundary();
      for (k1 = 0;k1 < k;k1 ++) {
	for (l1 = 0;l1 < l;l1 ++)
	  if (geo1.boundary(k1) == geo_index[l1]) break;
	if (l1 == l) {
	  assert(l < t_geo.n_geometry(i));
	  geo_index[l ++] = geo1.boundary(k1);
	}
      }
    }
    assert(l == t_geo.n_geometry(i));
    // reorder those indices according to their position in the template element
    j = t_geo.n_geometry(i);
    for (j1 = 0;j1 < j;j1 ++) {
      const GeometryBM& g = m.geometry(i, geo_index[j1]);
      k = g.n_vertex();
      std::vector<int> v(k);
      for (k1 = 0;k1 < k;k1 ++) {
	for (l1 = 0;l1 < i1;l1 ++) {
	  if (g.vertex(k1) == geo_img[0][l1]) {
	    v[k1] = l1;
	    break;
	  }
	}
	assert(l1 < i1); // ExcMeshData("the vertex is not belong to the element."));
      }
      sort(v.begin(), v.end());
      l = t_geo.n_geometry(i);
      l1 = g.index();
      for (k = 0;k < l;k ++) {
	if (v == t_geo.geometry(i, k).vertex()) {
	  //geo_img[i][j1] = l1;
	  /**
	   * 这里原来有一个严重的问题，搜索匹配有问题。
	   */
	  geo_img[i][k] = l1;
	  break;
	}
      }
      assert (k < l); //ExcMeshData("no such boundary in the template geometry."));
    }
  }

  /**
   * For save memory, if the dimension \f$i\f$ is known not used before hand,
   * we will release the memory.
   */
  const u_int& ef = sp->efficiencyFlag();
  for (i = 0;i <= TDIM;i ++) {
    if (!(ef & (1L<<i))) {
      geo_img[i].clear();
    }
  }
}

TEMPLATE
void THIS::lazyBuildGeometryImage()
{
  Geometry& geo = geometry();
  geo_img.resize(TDIM+1,std::vector<int>(1, 0));
  geo_img[TDIM].resize(1);
  geo_img[TDIM][0] = geo.index();
  geo_img[0] = geo.vertex();
}

TEMPLATE
void THIS::buildVertexArray(std::vector<Point<DOW> >& arr) const
{
  Mesh<DIM,DOW>& m = mesh();
  Geometry& geo = geometry();
  int n_vertex = geo.n_vertex();
  arr.resize(n_vertex);
  for (int i = 0;i < n_vertex;i ++)
    arr[i] = m.point(m.geometry(0,geo.vertex(i)).vertex(0));
}

TEMPLATE
const double ** THIS::buildVertexArray() const
{
  Mesh<DIM,DOW>& m = mesh();
  Geometry& geo = geometry();
  int n_vertex = geo.n_vertex();
  typedef double * double_pointer;
  const double ** arr = (const double **)new double_pointer[n_vertex];
  for (int i = 0;i < n_vertex;i ++)
    arr[i] = m.point(m.geometry(0,geo.vertex(i)).vertex(0));
  return arr;
}

TEMPLATE
const typename BasisFunction<value_type,DIM,TDIM>::Identity& THIS::basis_function_identity(int i) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.identity();
}

TEMPLATE
value_type THIS::basis_function_value(int i, const Point<DOW>& p) const
{
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.value(p, vertex_array);
}

TEMPLATE
std::vector<value_type> THIS::basis_function_value(const Point<DOW>& p) const
{
  int i, j;
  const double ** vertex_array = buildVertexArray();
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  std::vector<BasisFunction<value_type,DOW,DIM> >& basis_function = t_el.basisFunction();
  i = basis_function.size();
  std::vector<value_type> value(i);
  for (j = 0;j < i;j ++)
    value[j] = basis_function[j].value(p, vertex_array);
  delete[] vertex_array;
  return value;
}

TEMPLATE
std::vector<value_type> 
THIS::basis_function_value(int i, const std::vector<Point<DOW> >& p) const
{
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.value(p, vertex_array);
}

TEMPLATE
std::vector<std::vector<value_type> > 
THIS::basis_function_value(const std::vector<Point<DOW> >& p) const
{
  int i, i1;
  const double ** vertex_array = buildVertexArray();
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  std::vector<BasisFunction<value_type,DOW,DIM> >& basis_function = t_el.basisFunction();
  i = basis_function.size();
  std::vector<std::vector<value_type> > value(i);
  for (i1 = 0;i1 < i;i1 ++)
    value[i1] = basis_function[i1].value(p, vertex_array);
  delete[] vertex_array;
  return value;
}

TEMPLATE
std::vector<value_type> THIS::basis_function_gradient(int i, const Point<DOW>& p) const
{
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.gradient(p, vertex_array);
}

TEMPLATE
std::vector<std::vector<value_type> >
THIS::basis_function_gradient(const Point<DOW>& p) const
{
  int i, j;
  const double ** vertex_array = buildVertexArray();
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  std::vector<BasisFunction<value_type,DOW,DIM> >& basis_function = t_el.basisFunction();
  i = basis_function.size();
  std::vector<std::vector<value_type> > value(i);
  for (j = 0;j < i;j ++)
    value[j] = basis_function[j].gradient(p, vertex_array);
  delete[] vertex_array;
  return value;
}

TEMPLATE
std::vector<std::vector<value_type> >
THIS::basis_function_gradient(int i, const std::vector<Point<DOW> >& p) const
{
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.gradient(p, vertex_array);
}

TEMPLATE
std::vector<std::vector<std::vector<value_type> > >
THIS::basis_function_gradient(const std::vector<Point<DOW> >& p) const
{
  int i, i1;
  const double ** vertex_array = buildVertexArray();
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  std::vector<BasisFunction<value_type,DOW,DIM> >& basis_function = t_el.basisFunction();
  i = basis_function.size();
  std::vector<std::vector<std::vector<value_type> > > value(i);
  for (i1 = 0;i1 < i;i1 ++)
    value[i1] = basis_function[i1].gradient(p, vertex_array);
  delete[] vertex_array;
  return value;
}

#ifdef QUADRATIC_ELEMENT_SUPPORT

TEMPLATE
std::vector<std::vector<value_type> >
THIS::basis_function_hesse(int i, 
                           const Point<DOW>& p) const
{
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.hesse(p, vertex_array);
}

TEMPLATE
std::vector<std::vector<std::vector<value_type> > >
THIS::basis_function_hesse(const Point<DOW>& p) const
{
  int i, j;
  const double ** vertex_array = buildVertexArray();
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  std::vector<BasisFunction<value_type,DOW,DIM> >& basis_function = t_el.basisFunction();
  i = basis_function.size();
  std::vector<std::vector<value_type> > value(i);
  for (j = 0;j < i;j ++)
    value[j] = basis_function[j].hesse(p, vertex_array);
  delete[] vertex_array;
  return value;
}

TEMPLATE
std::vector<std::vector<std::vector<value_type> > >
THIS::basis_function_hesse(int i, 
                           const std::vector<Point<DOW> >& p) const
{
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  BasisFunction<value_type,DOW,DIM>& basis_function = t_el.basisFunction(i);
  return basis_function.hesse(p, vertex_array);
}

TEMPLATE
std::vector<std::vector<std::vector<std::vector<value_type> > > >
THIS::basis_function_hesse(const std::vector<Point<DOW> >& p) const
{
  int i, i1;
  const double ** vertex_array = buildVertexArray();
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  std::vector<BasisFunction<value_type,DOW,DIM> >& basis_function = t_el.basisFunction();
  i = basis_function.size();
  std::vector<std::vector<std::vector<value_type> > > value(i);
  for (i1 = 0;i1 < i;i1 ++)
    value[i1] = basis_function[i1].hesse(p, vertex_array);
  delete[] vertex_array;
  return value;
}

#endif // QUADRATIC_ELEMENT_SUPPORT

TEMPLATE
Point<DOW> THIS::local_to_global(const Point<TDIM>& lp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.local_to_global(lp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
Point<TDIM> THIS::global_to_local(const Point<DOW>& gp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.global_to_local(gp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
double THIS::local_to_global_jacobian(const Point<TDIM>& lp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.local_to_global_jacobian(lp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
double THIS::global_to_local_jacobian(const Point<DOW>& gp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.global_to_local_jacobian(gp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
std::vector<Point<DOW> >THIS::local_to_global(const std::vector<Point<TDIM> >& lp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.local_to_global(lp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
std::vector<Point<TDIM> > THIS::global_to_local(const std::vector<Point<DOW> >& gp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.global_to_local(gp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
std::vector<double> THIS::local_to_global_jacobian(const std::vector<Point<TDIM> >& lp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.local_to_global_jacobian(lp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
std::vector<double> THIS::global_to_local_jacobian(const std::vector<Point<DOW> >& gp) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  const CoordTransform<TDIM,DOW>& coord_transform = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return coord_transform.global_to_local_jacobian(gp, t_el.vertexArray(), vertex_array);
}

TEMPLATE
std::vector<double> THIS::unitOutNormal(const Point<DOW>& p, int n) const
{
  const double ** vertex_array = buildVertexArray();
  std::vector<double> ret(templateElement().unitOutNormal().value(p, vertex_array, n));
  delete[] vertex_array;
  return ret;
}

TEMPLATE
std::vector<std::vector<double> > THIS::unitOutNormal(const std::vector<Point<DOW> >& p, int n) const
{
  const double ** vertex_array = buildVertexArray();
  std::vector<std::vector<double> > ret(templateElement().unitOutNormal().value(p, vertex_array, n));
  delete[] vertex_array;
  return ret;
}

TEMPLATE
const QuadratureInfo<TDIM>& THIS::findQuadratureInfo(int i) const
{
  TemplateElement<value_type,DOW,TDIM>& t_el = templateElement();
  return t_el.findQuadratureInfo(i);
}

#undef THIS
#undef TEMPLATE

///////////////////////////////////////////////////////////////////////////

#define TEMPLATE template <class value_type, int DIM, int DOW, int TDIM>
#define THIS FEMSpace<value_type,DIM,DOW,TDIM>

TEMPLATE
THIS::FEMSpace(Mesh<DIM,DOW>& m,
               std::vector<template_t>& t) :
msh(&m),
  tmp_ele(&t),
  effi_flag(0xFFFF)
{}

TEMPLATE
THIS::FEMSpace(const fe_space_t& f) :
msh(f.msh),
  tmp_ele(f.tmp_ele),
  df(f.df),
  df_in(f.df_in),
  effi_flag(f.effi_flag)
{}

TEMPLATE
THIS::~FEMSpace()
{}

TEMPLATE
typename THIS::fe_space_t& THIS::operator=(const fe_space_t& f)
{
  if (&f != NULL) {
    msh = f.msh;
    tmp_ele = f.tmp_ele;
    df = f.df;
    df_in = f.df_in;
    effi_flag = f.effi_flag;
  }
  return *this;
}

  TEMPLATE
  void THIS::reinit(Mesh<DIM,DOW>& m,
                    std::vector<template_t>& t)
  {
    msh = &m;
    tmp_ele = &t;
  }

#ifdef MULTITHREAD

TEMPLATE
void THIS::buildElement(bool lazy)
{
  int n_thread = getThread();
  ThreadManager thread;
  for (int rank = 1;rank < n_thread;rank ++) {
    thread.start(encapsulate(&THIS::
			     threadBuildElement).collectArgs(this, lazy, n_thread, rank));
  }
  threadBuildElement(lazy, n_thread, 0);
  thread.join(encapsulate(&THIS::threadBuildElement));
}

TEMPLATE
void THIS::threadBuildElement(bool lazy,
                              int n_thread,
                              int rank)
{
  /*   std::cerr << "Entrying threadBuildElement with " */
  /* 	    << "lazy = " << lazy */
  /* 	    << ", n_thread = " << n_thread */
  /* 	    << ", rank = " << rank */
  /* 	    << std::endl; */
  int n = n_element()/n_thread;
  ElementIterator it_el = beginElement();
  it_el += n*rank;
  ElementIterator end_el = it_el;
  if (rank + 1 == n_thread)
    end_el = endElement();
  else
    end_el += n;
  if (lazy) {
    for (;it_el < end_el;++ it_el) {
      it_el->lazyBuildGeometryImage();
    }
  }
  else {
    for (;it_el < end_el;++ it_el) {
      it_el->buildGeometryImage();
    }
  }
}

#else

TEMPLATE
void THIS::buildElement(bool lazy)
{
  ElementIterator it_el = beginElement();
  ElementIterator end_el = endElement();
  if (lazy) {
    for (;it_el != end_el;++ it_el) {
      it_el->lazyBuildGeometryImage();
    }
  }
  else {
    for (;it_el != end_el;++ it_el) {
      it_el->buildGeometryImage();
    }
  }
}

#endif

TEMPLATE
void THIS::buildDofBoundaryMark()
{
  int i, j, k;
	
  const DegreeOfFreedom& d = dof();
  int n = n_dof();
  for (i = 0;i < n;i ++) {
    j = d.dof_index[i].dimension;
    k = d.dof_index[i].geometry_index;
    dofBoundaryMark(i) = mesh().boundaryMark(j, k);
  }
}

#ifdef MULTITHREAD

TEMPLATE
void THIS::buildDof()
{
  int i, j, k, l, m, n, i1, j1;
  std::vector<std::vector<bool> > flag;

  std::cerr << "Building degree of freedom for the FEM space ..." << std::endl;
  dof().n_geometry_dof.resize(TDIM+1);
  dof().geometry_dof.resize(TDIM+1);
  flag.resize(TDIM+1);
  for (i = 0;i <= TDIM;i ++) {
    /**
     * If the dimension \f$i\f$ is known not to be used before hand,
     * we will not allocate memeory for it
     * 
     */
    if (!(efficiencyFlag() & (1L<<i))) continue;

    j = mesh().n_geometry(i);
    dof().n_geometry_dof[i].resize(j, 0);
    flag[i].resize(j, false);
    dof().geometry_dof[i].resize(j);
  }

  dof().n_dof = 0;

  pthread_mutex_t mutex;
  pthread_mutex_init(&mutex, NULL);
  int n_thread = getThread();
  ThreadManager thread;
  for (int rank = 1;rank < n_thread;rank ++) {
    thread.start(encapsulate(&THIS::
			     threadBuildDof0).collectArgs(this, flag, mutex, n_thread, rank));
  }
  threadBuildDof0(flag, mutex, n_thread, 0);
  thread.join(encapsulate(&THIS::threadBuildDof0));
  pthread_mutex_destroy(&mutex);
  
  dof().dof_index.resize(dof().n_dof);
  dofInfo().resize(dof().n_dof);
	
  pthread_mutex_init(&mutex, NULL);
  for (int rank = 1;rank < n_thread;rank ++) {
    thread.start(encapsulate(&THIS::
			     threadBuildDof1).collectArgs(this, flag, mutex, n_thread, rank));
  }
  threadBuildDof1(flag, mutex, n_thread, 0);
  thread.join(encapsulate(&THIS::threadBuildDof1));
  pthread_mutex_destroy(&mutex);
  
  std::cerr << "\ttotal " << dof().n_dof << " degree of freedom found." << std::endl;
}

TEMPLATE
void THIS::threadBuildDof0(std::vector<std::vector<bool> >& flag,
                           pthread_mutex_t& mutex,
                           int n_thread,
                           int rank)
{
  int i, j, k, l, m, n;
  n = n_element()/n_thread;
  ElementIterator it_el = beginElement();
  it_el += n*rank;
  ElementIterator end_el = it_el;
  if (rank + 1 == n_thread)
    end_el = endElement();
  else
    end_el += n;
  for (;it_el < end_el;++ it_el) {
    TemplateElement<value_type,DOW,TDIM>& t_el = it_el->templateElement();
    const TemplateGeometry<TDIM>& t_el_geo = t_el.geometry();
    const TemplateDOF<TDIM>& t_dof = t_el.dof();
    const std::vector<std::vector<int> >& el_geo_img = it_el->geometryImage();
    it_el->dof().resize(t_dof.n_dof, -1);
    for (i = 0;i <= TDIM;i ++) {
      /**
       * If the dimension \f$i\f$ is known not to be used before hand,
       * we will not operate on this dimension at all
       * 
       */
      if (!(efficiencyFlag() & (1L<<i))) continue;

      for (j = 0;j < t_el_geo.n_geometry(i);j ++) {
	l = el_geo_img[i][j];
	m = t_dof.n_geometry_dof[i][j];
	pthread_mutex_lock(&mutex);
	if (flag[i][l]) {
	  Assert(dof().n_geometry_dof[i][l] == m, ExcDOFData("DOF number mismatch."));
	}
	else {
	  flag[i][l] = true;
	  dof().n_geometry_dof[i][l] = m;
	  dof().geometry_dof[i][l].resize(m,0);
	  for (n = 0;n < m;n ++) {
	    dof().geometry_dof[i][l][n] = dof().n_dof ++;
	  }
	}
	pthread_mutex_unlock(&mutex);
      }
    }
  }
}


TEMPLATE
void THIS::threadBuildDof1(std::vector<std::vector<bool> >& flag,
                           pthread_mutex_t& mutex,
                           int n_thread,
                           int rank)
{
  int i, j, k, l, m, n, i1, j1;
  n = n_element()/n_thread;
  ElementIterator it_el = beginElement();
  it_el += n*rank;
  ElementIterator end_el = it_el;
  if (rank + 1 == n_thread)
    end_el = endElement();
  else
    end_el += n;
  for (;it_el < end_el;++ it_el) {
    TemplateElement<value_type,DOW,TDIM>& t_el = it_el->templateElement();
    std::vector<int>& el_dof = it_el->dof();
    const TemplateGeometry<TDIM>& t_el_geo = t_el.geometry();
    const TemplateDOF<TDIM>& t_dof = t_el.dof();
    const std::vector<BasisFunction<value_type,DOW,TDIM> >& basis_function = t_el.basisFunction();
    const std::vector<std::vector<int> >& el_geo_img = it_el->geometryImage();
    double h = (mesh().point(mesh().geometry(0, it_el->geometry().vertex(0)).vertex(0))
                - mesh().point(mesh().geometry(0, it_el->geometry().vertex(1)).vertex(0))).length();
    for (i = 0;i <= TDIM;i ++) {
      /**
       * If the dimension \f$i\f$ is known not to be used before hand,
       * we will not operate on this dimension at all
       * 
       */
      if (!(efficiencyFlag() & (1L<<i))) continue;

      for (j = 0;j < t_el_geo.n_geometry(i);j ++) {
	l = el_geo_img[i][j];
	m = t_dof.n_geometry_dof[i][j];
	pthread_mutex_lock(&mutex);
	if (!flag[i][l]) {
	  pthread_mutex_unlock(&mutex);
	  for (n = 0;n < m;n ++) {
	    k = t_dof.geometry_dof[i][j][n];
	    Point<DOW> ip = it_el->local_to_global(basis_function[k].interpPoint());
	    const typename BasisFunction<value_type,DIM,TDIM>::Identity& id = it_el->basis_function_identity(k);
	    for (i1 = 0;i1 < m;i1 ++) {
	      j1 = dof().geometry_dof[i][l][i1];
	      const Point<DOW>& ip1 = dofInfo(j1).interp_point;
	      const typename BasisFunction<value_type,DIM,TDIM>::Identity& id1 = dofInfo(j1).identity;
	      if ((ip - ip1).length() < (1.e-6)*h && id == id1) {
		el_dof[k] = j1;
		break;
	      }
	    }
	    Assert(i1 < m, ExcDOFData("DOF data mismatch."));
	  }
	}
	else {
	  flag[i][l] = false;
	  for (n = 0;n < m;n ++) {
	    k = t_dof.geometry_dof[i][j][n];
	    j1 = dof().geometry_dof[i][l][n];
	    el_dof[k] = j1;
	    dof().dof_index[j1].dimension = i;
	    dof().dof_index[j1].geometry_index = l;
	    dof().dof_index[j1].dof_index = n;
	    dofInfo(j1).interp_point = it_el->local_to_global(basis_function[k].interpPoint());
	    dofInfo(j1).identity = it_el->basis_function_identity(k);
	  }
	  pthread_mutex_unlock(&mutex);
	}
      }
    }
  }
}

#else

TEMPLATE
void THIS::buildDof()
{
  int i, j, k, l, m, n, i1, j1;
  std::vector<std::vector<bool> > flag;

  std::cerr << "Building degree of freedom for the FEM space ..." << std::endl;
  dof().n_geometry_dof.resize(TDIM+1);
  dof().geometry_dof.resize(TDIM+1);
  flag.resize(TDIM+1);
  for (i = 0;i <= TDIM;i ++) {
    /**
     * If the dimension \f$i\f$ is known not to be used before hand,
     * we will not allocate memeory for it
     * 
     */
    if (!(efficiencyFlag() & (1L<<i))) continue;

    j = mesh().n_geometry(i);
    dof().n_geometry_dof[i].resize(j, 0);
    flag[i].resize(j, false);
    dof().geometry_dof[i].resize(j);
  }

  dof().n_dof = 0;
  ElementIterator it_el = beginElement();
  ElementIterator end_el = endElement();
  for (;it_el != end_el;++ it_el) {
    TemplateElement<value_type,DOW,TDIM>& t_el = it_el->templateElement();
    const TemplateGeometry<TDIM>& t_el_geo = t_el.geometry();
    const TemplateDOF<TDIM>& t_dof = t_el.dof();
    const std::vector<std::vector<int> >& el_geo_img = it_el->geometryImage();
    it_el->dof().resize(t_dof.n_dof);
    for (i = 0;i <= TDIM;i ++) {
      /**
       * If the dimension \f$i\f$ is known not to be used before hand,
       * we will not operate on this dimension at all
       * 
       */
      if (!(efficiencyFlag() & (1L<<i))) continue;

      for (j = 0;j < t_el_geo.n_geometry(i);j ++) {
	l = el_geo_img[i][j];
	m = t_dof.n_geometry_dof[i][j];
	if (flag[i][l]) {
	  Assert(dof().n_geometry_dof[i][l] == m, ExcDOFData("DOF number mismatch."));
	}
	else {
	  flag[i][l] = true;
	  dof().n_geometry_dof[i][l] = m;
	  dof().geometry_dof[i][l].resize(m);
	  for (n = 0;n < m;n ++) {
	    dof().geometry_dof[i][l][n] = dof().n_dof ++;
	  }
	}
      }
    }
  }

  dof().dof_index.resize(dof().n_dof);
  dofInfo().resize(dof().n_dof);
	
  it_el = beginElement();
  for (;it_el != end_el;++ it_el) {
    TemplateElement<value_type,DOW,TDIM>& t_el = it_el->templateElement();
    std::vector<int>& el_dof = it_el->dof();
    const TemplateGeometry<TDIM>& t_el_geo = t_el.geometry();
    const TemplateDOF<TDIM>& t_dof = t_el.dof();
    const std::vector<BasisFunction<value_type,DOW,TDIM> >& basis_function = t_el.basisFunction();
    const std::vector<std::vector<int> >& el_geo_img = it_el->geometryImage();
    double h = (mesh().point(mesh().geometry(0, it_el->geometry().vertex(0)).vertex(0))
                - mesh().point(mesh().geometry(0, it_el->geometry().vertex(1)).vertex(0))).length();
    for (i = 0;i <= TDIM;i ++) {
      /**
       * If the dimension \f$i\f$ is known not to be used before hand,
       * we will not operate on this dimension at all
       * 
       */
      if (!(efficiencyFlag() & (1L<<i))) continue;

      for (j = 0;j < t_el_geo.n_geometry(i);j ++) {
	l = el_geo_img[i][j];
	m = t_dof.n_geometry_dof[i][j];
	if (!flag[i][l]) {
	  for (n = 0;n < m;n ++) {
	    k = t_dof.geometry_dof[i][j][n];
	    Point<DOW> ip = it_el->local_to_global(basis_function[k].interpPoint());
	    const typename BasisFunction<value_type,DIM,TDIM>::Identity& id = it_el->basis_function_identity(k);
	    for (i1 = 0;i1 < m;i1 ++) {
	      j1 = dof().geometry_dof[i][l][i1];
	      const Point<DOW>& ip1 = dofInfo(j1).interp_point;
	      const typename BasisFunction<value_type,DIM,TDIM>::Identity& id1 = dofInfo(j1).identity;
	      if ((ip - ip1).length() < (1.e-6)*h && id == id1) {
		el_dof[k] = j1;
		break;
	      }
	    }
	    Assert(i1 < m, ExcDOFData("DOF data mismatch."));
	  }
	}
	else {
	  flag[i][l] = false;
	  for (n = 0;n < m;n ++) {
	    k = t_dof.geometry_dof[i][j][n];
	    j1 = dof().geometry_dof[i][l][n];
	    el_dof[k] = j1;
	    dof().dof_index[j1].dimension = i;
	    dof().dof_index[j1].geometry_index = l;
	    dof().dof_index[j1].dof_index = n;
	    dofInfo(j1).interp_point = it_el->local_to_global(basis_function[k].interpPoint());
	    dofInfo(j1).identity = it_el->basis_function_identity(k);
	  }
	}
      }
    }
  }
  std::cerr << "\ttotal " << dof().n_dof << " degree of freedom found." << std::endl;
}

#endif

TEMPLATE
void THIS::updateDofInterpPoint()
{
  ElementIterator it_el = beginElement();
  ElementIterator end_el = endElement();
  for (;it_el != end_el;++ it_el) {
    TemplateElement<value_type,DOW,TDIM>& t_el = it_el->templateElement();
    std::vector<int>& el_dof = it_el->dof();
    int n_el_dof = el_dof.size();
    const std::vector<BasisFunction<value_type,DOW,TDIM> >& basis_function = t_el.basisFunction();
    for (int j = 0;j < n_el_dof;j ++) {
      dofInfo(el_dof[j]).interp_point = it_el->local_to_global(basis_function[j].interpPoint());
    }
  }
}

#undef THIS
#undef TEMPLATE

///////////////////////////////////////////////////////////////////////////

#define TEMPLATE template <class value_type, int DIM, int DOW, int TDIM, typename Number>
#define THIS FEMFunction<value_type,DIM,DOW,TDIM,Number>

TEMPLATE
THIS::FEMFunction(fe_space_t & f) :
sp(&f)
{
  if (sp != NULL)
    Vector<Number>::reinit(sp->n_dof());
}

TEMPLATE
THIS::~FEMFunction()
{}

TEMPLATE
typename THIS::fe_space_t& THIS::femSpace()
{
  return *sp;
}

TEMPLATE
void THIS::reinit(fe_space_t& f, bool is_bare)
{
  sp = &f;
  if ((sp != NULL) && (! is_bare))
    Vector<Number>::reinit(sp->n_dof());
}

TEMPLATE
value_type THIS::value(const Point<DOW>& p, const element_t& e) const
{
  int i, j;
  value_type val;
  const std::vector<int>& dof = e.dof();
  std::vector<value_type> basis_value = e.basis_function_value(p);
  j = dof.size();
  for (i = 0, val = 0;i < j;i ++) {
    val += (*this)(dof[i])*basis_value[i];
  }
  return val;
}

TEMPLATE
std::vector<value_type> THIS::gradient(const Point<DOW>& p,
                                       const element_t& e) const
{
  int i, j, k;
  std::vector<value_type> val(DOW, 0);
  const std::vector<int>& dof = e.dof();
  std::vector<std::vector<value_type> > basis_gradient = e.basis_function_gradient(p);
  j = dof.size();
  for (i = 0;i < j;i ++) {
    for (k = 0;k < DOW;k ++) {
      val[k] += (*this)(dof[i])*basis_gradient[i][k];
    }
  }
  return val;
}

TEMPLATE
std::vector<value_type> THIS::value(const std::vector<Point<DOW> >& p,
                                    const element_t& e) const
{
  int i, j, i1, j1;
  i = p.size();
  std::vector<value_type> val(i, 0);
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  std::vector<std::vector<value_type> > basis_value = e.basis_function_value(p);
  for (i1 = 0;i1 < i;i1 ++)
    for (j1 = 0;j1 < j;j1 ++)
      val[i1] += (*this)(dof[j1])*basis_value[j1][i1];
  return val;
}

TEMPLATE
std::vector<std::vector<value_type> > 
THIS::gradient(const std::vector<Point<DOW> >& p,
               const element_t& e) const
{
  int i, j, i1, j1, m;
  i = p.size();
  std::vector<std::vector<value_type> > val(i, std::vector<value_type>(DOW, 0));
  const std::vector<int>& dof = e.dof();
  std::vector<std::vector<std::vector<value_type> > > basis_gradient = e.basis_function_gradient(p);
  j = dof.size();
  for (i1 = 0;i1 < i;i1 ++) {
    for (j1 = 0;j1 < j;j1 ++)
      for (m = 0;m < DOW;m ++)
	val[i1][m] += (*this)(dof[j1])*basis_gradient[j1][i1][m];
  }
  return val;
}


/**
 * for optimization. Assume the basis function value and gradient are known,
 * the following functions will speedup the calculation to get the value and
 * gradient of the finite element function.
 */
TEMPLATE
value_type THIS::value(const std::vector<value_type>& basis_value, const element_t& e) const
{
  int i, j;
  value_type val;
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i = 0, val = 0;i < j;i ++) {
    val += (*this)(dof[i])*basis_value[i];
  }
  return val;
}

TEMPLATE
std::vector<value_type> THIS::gradient(const std::vector<std::vector<value_type> >& basis_gradient,
                                       const element_t& e) const
{
  int i, j, k;
  std::vector<value_type> val(DOW, 0);
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i = 0;i < j;i ++) {
    for (k = 0;k < DOW;k ++) {
      val[k] += (*this)(dof[i])*basis_gradient[i][k];
    }
  }
  return val;
}

TEMPLATE
std::vector<value_type> THIS::value(const std::vector<std::vector<value_type> >& basis_value,
                                    const element_t& e) const
{
  int i, j, i1, j1;
  i = basis_value[0].size();
  std::vector<value_type> val(i, 0);
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i1 = 0;i1 < i;i1 ++)
    for (j1 = 0;j1 < j;j1 ++)
      val[i1] += (*this)(dof[j1])*basis_value[j1][i1];
  return val;
}

TEMPLATE
std::vector<std::vector<value_type> > 
THIS::gradient(const std::vector<std::vector<std::vector<value_type> > >& basis_gradient,
               const element_t& e) const
{
  int i, j, i1, j1, m;
  i = basis_gradient[0].size();
  std::vector<std::vector<value_type> > val(i, std::vector<value_type>(DIM, 0));
  const std::vector<int>& dof = e.dof();
  j = dof.size();
  for (i1 = 0;i1 < i;i1 ++) {
    for (j1 = 0;j1 < j;j1 ++)
      for (m = 0;m < DOW;m ++)
	val[i1][m] += (*this)(dof[j1])*basis_gradient[j1][i1][m];
  }
  return val;
}

// end of the optimization code

TEMPLATE
void THIS::loadData(const std::string& filename)
{
  std::ifstream is(filename.c_str());
  Vector<Number>::block_read(is);
  is.close();
}


TEMPLATE
void THIS::writeData(const std::string& filename)
{
  std::ofstream os(filename.c_str());
  Vector<Number>::block_write(os);
  os.close();
}

#undef THIS
#undef TEMPLATE

///////////////////////////////////////////////////////////////////////////

template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::LocalFEMFunction(Element<value_type,DIM,DOW,TDIM> & e)
{
  ele = &e;
  if (ele != NULL)
    Vector<value_type>::reinit(ele->n_dof());
}


template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  Element<value_type,DIM,DOW,TDIM>& LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::element()
{
  return *ele;
}


template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  void LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::reinit(Element<value_type,DIM,DOW,TDIM> & e)
{
  ele = &e;
  if (ele != NULL)
    Vector<value_type>::reinit(ele->n_dof());
}


template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  value_type LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::value(const Point<DOW>& p) const
{
  int i, j;
  value_type val;
  j = Vector<Number>::size();
  std::vector<value_type> basis_value = ele->basis_function_value(p);
  for (i = 0, val = 0;i < j;i ++) {
    val += (*this)(i)*basis_value[i];
  }
  return val;
}



template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  std::vector<value_type> LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::value(const std::vector<Point<DOW> >& p) const
{
  int i, j, i1, j1;
  i = p.size();
  std::vector<value_type> val(i, 0);
  j = Vector<Number>::size();
  std::vector<std::vector<value_type> > basis_value = ele->basis_function_value(p);
  for (i1 = 0;i1 < i;i1 ++)
    for (j1 = 0;j1 < j;j1 ++)
      val[i1] += (*this)(j1)*basis_value[j1][i1];
  return val;
}


template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  std::vector<value_type> LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::gradient(const Point<DOW>& p) const
{
  int i, j, k;
  std::vector<value_type> val(DOW, 0);
  j = Vector<Number>::size();
  std::vector<std::vector<value_type> > basis_gradient = ele->basis_function_gradient(p);
  for (i = 0;i < j;i ++) {
    for (k = 0;k < DOW;k ++) {
      val[k] += (*this)(i)*basis_gradient[i][k];
    }
  }
  return val;
}

template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  std::vector<std::vector<value_type> > 
  LocalFEMFunction<value_type,DIM,DOW,TDIM,Number>::gradient(const std::vector<Point<DOW> >& p) const
{
  int i, j, i1, j1, m;
  i = p.size();
  std::vector<std::vector<value_type> > val(i);
  j = Vector<Number>::size();
  std::vector<std::vector<std::vector<value_type> > > basis_gradient = ele->basis_function_gradient(p);
  for (i1 = 0;i1 < i;i1 ++) {
    val[i1].resize(DOW, 0);
    for (j1 = 0;j1 < j;j1 ++)
      for (m = 0;m < DOW;m ++)
	val[i1][m] += (*this)(j1)*basis_gradient[j1][i1][m];
  }
  return val;
}


///////////////////////////////////////////////////////////////////////////


template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  void BoundaryConditionAdmin<value_type,DIM,DOW,TDIM,Number>::add(const BoundaryCondition<value_type,DIM,DOW,TDIM,Number>& b)
{
  if (b.boundaryType() != BoundaryCondition<value_type,DIM,DOW,TDIM,Number>::DIRICHLET) {
    std::cerr << "Now we can only apply Dirichlet boundary condition." << std::endl;
    Assert(false, ExcInternalError());
  }
  if (b.boundaryMark() < 0) {
    std::cerr << "We now require a boundary mark to be a positive number."
	      << std::endl;
    Assert(false, ExcInternalError());
  }
  typename std::vector<const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *>::iterator 
    the_bc = std::vector<const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *>::begin();
  typename std::vector<const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *>::iterator 
    end_bc = std::vector<const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *>::end();
  for (;the_bc != end_bc;the_bc ++) {
    if ((*the_bc)->boundaryMark() == b.boundaryMark()) {
      std::cerr << "There is a boundary condition for the same boundary mark("
		<< b.boundaryMark()
		<< ") already."
		<< std::endl;
      Assert(false, ExcInternalError());
    }
  }
  this->push_back(&b);
  int i = index_map.size();
  if (i <= b.boundaryMark())
    for (;i <= b.boundaryMark();i ++)
      index_map.push_back(-1);
  index_map[b.boundaryMark()] = std::vector<const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> *>::size() - 1;
}

template <class value_type, int DIM, int DOW, int TDIM, typename Number>
void BoundaryConditionAdmin<value_type,DIM,DOW,TDIM,Number>::apply(SparseMatrix<double>& A, 
								   Vector<double>& u, 
								   Vector<double>& f, 
								   bool preserve_symmetry)
{
    unsigned int i, j, k, l, n_dof;
    n_dof = fem_space->n_dof();
    Assert (A.n() == A.m(), ExcDimensionMismatch(A.n(), A.m()));
    Assert (A.n() == f.size(), ExcDimensionMismatch(A.n(), f.size()));
    Assert (A.n() == u.size(), ExcDimensionMismatch(A.n(), u.size()));
    Assert (A.n() == n_dof, ExcDimensionMismatch(A.n(), n_dof));
    typename FEMSpace<double,DIM,DOW,TDIM>::bmark_t bm;
//  const SparsityPattern& spA = A.get_sparsity_pattern();
//  const std::size_t * row_start = spA.get_rowstart_indices();
//  const unsigned int * column_index = spA.get_column_numbers();
    for (i = 0;i < n_dof;i ++)
    {
	bm = fem_space->dofBoundaryMark(i);
	const BoundaryCondition<value_type,DIM,DOW,TDIM,Number> * bc = find(bm);
	//if (!isValid(bc)) continue;
	if (bc == NULL) continue;

	SparseMatrix<double>::iterator row_iterator = A.begin(i);
	SparseMatrix<double>::iterator row_end = A.end(i);
	double diag = row_iterator->value();
	double bnd_value = bc->value(fem_space->dofInfo(i).interp_point);
	f(i) = diag * bnd_value;
	for (++row_iterator; row_iterator != row_end; ++row_iterator)
	{
	    row_iterator->value() = 0.0;
	}

	if (preserve_symmetry)
	{
	    row_iterator = A.begin(i);
	    for (++row_iterator; row_iterator != row_end; ++row_iterator)
	    {
		k = row_iterator->column();
		SparseMatrix<double>::iterator col_iterator = A.begin(k);
		SparseMatrix<double>::iterator col_end = A.end(k);
		for (++col_iterator; col_iterator != col_end; ++col_iterator)
		    if (col_iterator->column() == i)
			break;
		if (col_iterator == col_end)
		{
		    std::cerr << "Boundary condition applying error!" << std::endl;
		    exit(-1);
		}
		f(k) -= col_iterator->value() * bnd_value;
		col_iterator->value() = 0.0;
	    }
	}
    
	// f(i) = A.diag_element(i)*bc.value(fem_space->dofInfo(i).interp_point);
	// for (j = row_start[i]+1;j < row_start[i+1];j ++) A.global_entry(j) = 0.0;
	// if (preserve_symmetry) {
	//   for (j = row_start[i]+1;j < row_start[i+1];j ++) {
	// 	k = column_index[j];
	// 	const unsigned int * p = std::find(&column_index[row_start[k]+1],
	// 					   &column_index[row_start[k+1]], i);
	// 	if (p != &column_index[row_start[k+1]]) {
	// 	  l = p - &column_index[row_start[0]];
	// 	  f(k) -= A.global_entry(l) * f(i) / A.diag_element(i);
	// 	  A.global_entry(l) = 0.;
	// 	}
	//   }
	// }
    }
}



template <class value_type, int DIM, int DOW, int TDIM, typename Number>
  void BoundaryConditionAdmin<value_type,DIM,DOW,TDIM,Number>::clearEntry(Vector<double>& f)
{
  unsigned int i, n_dof;
  n_dof = fem_space->n_dof();
  typename FEMSpace<double,DIM,DOW,TDIM>::bmark_t bm;
  for (i = 0;i < n_dof;i ++) {
    bm = fem_space->dofBoundaryMark(i);
    if (bm == 0) continue;
    f(i) = 0.0;
  }
}

AFEPACK_CLOSE_NAMESPACE

#endif //_FEMSpace_templates_h_

/**
 * end of file
 * 
 */

