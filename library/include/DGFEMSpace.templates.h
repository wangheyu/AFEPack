///////////////////////////////////////////////////////////////////////////
// DGFEMSpace.templates.h : by R.Lie
//

#ifndef _DGFEMSpace_templates_h_
#define _DGFEMSpace_templates_h_

#include "DGFEMSpace.h"

AFEPACK_OPEN_NAMESPACE

#define TEMPLATE template <class value_type, int DIM, int DOW, int TDIM, int TDIM1>
#define THIS DGElement<value_type,DIM,DOW,TDIM,TDIM1>

TEMPLATE
THIS::DGElement(fe_space_t& f) : sp(&f) 
{
  neigh[0] = NULL;
  neigh[1] = NULL;
};

TEMPLATE
THIS::DGElement(const dg_element_t& e) :
sp(e.sp),
  geometry_index(e.geometry_index),
  template_element_index(e.template_element_index),
  geo_img(e.geo_img)
{
  neigh[0] = e.neigh[0];
  neigh[1] = e.neigh[1];
};

TEMPLATE
THIS::~DGElement()
{};

TEMPLATE
void THIS::reinit(fe_space_t& f, 
                  const int& g, 
                  const int& t)
{
  sp = &f;
  geometry_index = g;
  template_element_index = t;
};

TEMPLATE
void THIS::reinit(const int& g, const int& t)
{
  Assert (sp != NULL, ExcInternalError());
  geometry_index = g;
  template_element_index = t;
};

TEMPLATE
typename THIS::dg_element_t& THIS::operator=(const dg_element_t& e)
{
  if (&e != NULL) {
    sp = e.sp;
    geometry_index = e.geometry_index;
    template_element_index = e.template_element_index;
    geo_img = e.geo_img;
    neigh[0] = e.neigh[0];
    neigh[1] = e.neigh[1];
  }
  return *this;
};

TEMPLATE
const typename THIS::fe_space_t& THIS::femSpace() const
{
  return *sp;
};

TEMPLATE
typename THIS::fe_space_t& THIS::femSpace()
{
  return *sp;
};

TEMPLATE
const int& THIS::index() const
{
  return geometry_index;
};

TEMPLATE
int& THIS::index()
{
  return geometry_index;
};

TEMPLATE
const int& THIS::boundaryIndex(const u_int& i) const
{
  return bnd_idx[i];
};

TEMPLATE
int& THIS::boundaryIndex(const u_int& i)
{
  return bnd_idx[i];
};

TEMPLATE
GeometryBM& THIS::geometry() const
{
  return sp->mesh().geometry(TDIM1, geometry_index);
};

TEMPLATE
void THIS::geometry(const GeometryBM& g)
{
  geometry_index = g.index();
};

TEMPLATE
void THIS::geometry(const int& i)
{
  geometry_index = i;
};

TEMPLATE
typename THIS::bmark_t THIS::boundaryMark() const
{
  return sp->mesh().geometry(TDIM1, geometry_index).boundaryMark();
};

TEMPLATE
TemplateDGElement<TDIM1,DOW>& THIS::templateElement() const
{
  return sp->templateDGElement(template_element_index);
};

TEMPLATE
void THIS::templateElement(const int& i)
{
  template_element_index = i;
};

TEMPLATE
TemplateGeometry<TDIM1>& THIS::templateGeometry() const
{
  return sp->templateDGElement(template_element_index).geometry();
};

TEMPLATE
const std::vector<std::vector<int> >& THIS::geometryImage() const
{
  return geo_img;
};

TEMPLATE
std::vector<std::vector<int> >& THIS::geometryImage()
{
  return geo_img;
};

TEMPLATE
void THIS::buildGeometryImage()
{
  Mesh<DIM,DOW>& m = sp->mesh();
  Geometry& geo = this->geometry();
  TemplateDGElement<TDIM1,DOW>& t_el = templateElement();
  TemplateGeometry<TDIM1>& t_geo = t_el.geometry();
	
  int i, j, k, l, i1, j1, k1, l1;
	
  geo_img.resize(TDIM1+1);
  for (i = 0;i <= TDIM1;i ++)
    geo_img[i].resize(t_geo.n_geometry(i), -1);
  geo_img[TDIM1][0] = geo.index();
  geo_img[0] = geo.vertex();
  i1 = geo_img[0].size();
  for (i = TDIM1-1;i > 0;i --) {
    // collect the geometry index in dimension \p{i}
    j = t_geo.n_geometry(i+1);
    std::vector<int> geo_index(t_geo.n_geometry(i), -1);
    for (j1 = 0, l = 0;j1 < j;j1 ++) {
      Geometry& geo1 = m.geometry(i+1, geo_img[i+1][j1]);
      k = geo1.n_boundary();
      for (k1 = 0;k1 < k;k1 ++) {
        for (l1 = 0;l1 < l;l1 ++)
          if (geo1.boundary(k1) == geo_index[l1]) break;
        if (l1 == l) {
          Assert(l < t_geo.n_geometry(i), ExcInternalError());
          geo_index[l ++] = geo1.boundary(k1);
        }
      }
    }
    Assert(l == t_geo.n_geometry(i), ExcInternalError());
    // reorder those indices according their position in the template element
    j = t_geo.n_geometry(i);
    for (j1 = 0;j1 < j;j1 ++) {
      const Geometry& g = m.geometry(i, geo_index[j1]);
      k = g.n_vertex();
      std::vector<int> v(k);
      for (k1 = 0;k1 < k;k1 ++) {
        for (l1 = 0;l1 < i1;l1 ++) {
          if (g.vertex(k1) == geo_img[0][l1]) {
            v[k1] = l1;
            break;
          }
        }
        Assert(l1 < i1, ExcMeshData("the vertex is not belong to the element."));
      }
      sort(v.begin(), v.end());
      l = t_geo.n_geometry(i);
      l1 = g.index();
      for (k = 0;k < l;k ++) {
        if (v == t_geo.geometry(i, k).vertex()) {
          geo_img[i][j1] = l1;
          break;
        }
      }
      Assert(k < l,ExcMeshData("no such boundary in the template geometry."));
    }
  }
};

TEMPLATE
void THIS::buildVertexArray(std::vector<Point<DOW> >& arr) const
{
  Mesh<DIM,DOW>& m = sp->mesh();
  Geometry& geo = geometry();
  int n_vertex = geo.n_vertex();
  arr.resize(n_vertex);
  for (int i = 0;i < n_vertex;i ++)
    arr[i] = m.point(m.geometry(0,geo.vertex(i)).vertex(0));
};

TEMPLATE
const double ** THIS::buildVertexArray() const
{
  Mesh<DIM,DOW>& m = sp->mesh();
  Geometry& geo = geometry();
  int n_vertex = geo.n_vertex();
  typedef double * double_pointer;
  const double ** arr = (const double **)new double_pointer[n_vertex];
  for (int i = 0;i < n_vertex;i ++)
    arr[i] = m.point(m.geometry(0,geo.vertex(i)).vertex(0));
  return arr;
};

TEMPLATE
Point<DOW> THIS::local_to_global(const Point<TDIM1>& lp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.local_to_global(lp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
Point<TDIM1> THIS::global_to_local(const Point<DOW>& gp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.global_to_local(gp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
double THIS::local_to_global_jacobian(const Point<TDIM1>& lp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.local_to_global_jacobian(lp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
double THIS::global_to_local_jacobian(const Point<DOW>& gp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.global_to_local_jacobian(gp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
std::vector<Point<DOW> >THIS::local_to_global(const std::vector<Point<TDIM1> >& lp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.local_to_global(lp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
std::vector<Point<TDIM1> > THIS::global_to_local(const std::vector<Point<DOW> >& gp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.global_to_local(gp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
std::vector<double> THIS::local_to_global_jacobian(const std::vector<Point<TDIM1> >& lp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.local_to_global_jacobian(lp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
std::vector<double> THIS::global_to_local_jacobian(const std::vector<Point<DOW> >& gp) const
{
  dg_template_t& t_el = templateElement();
  const coord_trans_t& ct = t_el.coordTransform();
  std::vector<Point<DOW> > vertex_array;
  buildVertexArray(vertex_array);
  return ct.global_to_local_jacobian(gp, t_el.vertexArray(), vertex_array);
};

TEMPLATE
const typename THIS::element_t& THIS::neighbourElement(const int& i) const
{
  Assert(i >= 0 && i < 2, ExcInternalError());
  return *neigh[i];
};

TEMPLATE
typename THIS::element_t& THIS::neighbourElement(const int& i)
{
  Assert(i >= 0 && i < 2, ExcInternalError());
  return *neigh[i];
};

TEMPLATE
const typename THIS::element_t * THIS::p_neighbourElement(const int& i) const
{
  return neigh[i];
};

TEMPLATE
typename THIS::element_t * THIS::p_neighbourElement(const int& i)
{
  return neigh[i];
};

TEMPLATE
const QuadratureInfo<TDIM1>& THIS::findQuadratureInfo(const int& i) const
{
  dg_template_t& t_el = templateElement();
  return t_el.findQuadratureInfo(i);
};

TEMPLATE 
std::vector<double> 
unitOutNormal(const Point<DIM>& p, 
              const Element<value_type,DIM,DOW,TDIM>& ele, 
              const DGElement<value_type,DIM,DOW,TDIM,TDIM1>& dgele)
{
  if (dgele.neigh[0] == &ele)
    return ele.unitOutNormal(p, dgele.bnd_idx[0]);
  else if (dgele.neigh[1] == &ele)
    return ele.unitOutNormal(p, dgele.bnd_idx[1]);
};

TEMPLATE 
std::vector<std::vector<double> > 
unitOutNormal(const std::vector<Point<DIM> >& p, 
              const Element<value_type,DIM,DOW,TDIM>& ele, 
              const DGElement<value_type,DIM,DOW,TDIM,TDIM1>& dgele)
{
  if (dgele.neigh[0] == &ele)
    return ele.unitOutNormal(p, dgele.bnd_idx[0]);
  else if (dgele.neigh[1] == &ele)
    return ele.unitOutNormal(p, dgele.bnd_idx[1]);
};

#undef THIS
#undef TEMPLATE

///////////////////////////////////////////////////////////////////////////

#define TEMPLATE template <class value_type, int DIM, int DOW, int TDIM, int TDIM1>
#define THIS DGFEMSpace<value_type,DIM,DOW,TDIM,TDIM1>

TEMPLATE
THIS::DGFEMSpace(Mesh<DIM,DOW>& m,
                 std::vector<template_t>& t,
                 std::vector<dg_template_t>& tdg) :
base_t(m, t), tmp_dgele(&tdg)
{};

TEMPLATE
THIS::DGFEMSpace(const fe_space_t& f) :
base_t(f), tmp_dgele(f.tmp_dgele)
{};

TEMPLATE
THIS::~DGFEMSpace()
{};

TEMPLATE
THIS& THIS::operator=(const THIS& f)
{
  base_t::operator=(f);
  if (&f != NULL) 
    tmp_dgele = f.tmp_dgele;
  return *this;
};

TEMPLATE
void THIS::reinit(Mesh<DIM,DOW>& m,
                  std::vector<template_t>& t,
                  std::vector<dg_template_t>& tdg)
{
  base_t::reinit(m,t);
  tmp_dgele = &tdg;
};

TEMPLATE
void THIS::buildDGElement()
{
  int i, j;
  std::vector<int> index(base_t::mesh().n_geometry(DIM-1), -1);
  DGElementIterator it_el = beginDGElement();
  DGElementIterator end_el = endDGElement();
  for (i = 0;it_el != end_el;++ it_el) {
    it_el->buildGeometryImage();
    it_el->neigh[0] = NULL;
    it_el->neigh[1] = NULL;
    index[it_el->index()] = i ++;
  }
  typename base_t::ElementIterator 
    the_ele = base_t::beginElement(),
    end_ele = base_t::endElement();
  for (;the_ele != end_ele;++ the_ele) {
    GeometryBM& g = the_ele->geometry();
    for (i = 0;i < g.n_boundary();i ++) {
      j = index[g.boundary(i)];
      if (j != -1) {
        if (dgele[j].neigh[0] == NULL) {
          dgele[j].neigh[0] = &(*the_ele);
          dgele[j].bnd_idx[0] = i;
        }
        else if (dgele[j].neigh[1] == NULL) {
          dgele[j].neigh[1] = &(*the_ele);
          dgele[j].bnd_idx[1] = i;
        }
        else {
          Assert(false, ExcInternalError());
        }
      }
    }
  }
};

#undef THIS
#undef TEMPLATE

AFEPACK_CLOSE_NAMESPACE

#endif //_DGFEMSpace_templates_h_

//
// end of file
///////////////////////////////////////////////////////////////////////////
