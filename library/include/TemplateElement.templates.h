/**
 * @file   TemplateElement.templates.h
 * @author Robert Lie
 * @date   Wed Jun 21 08:19:19 2003
 * 
 * @brief  
 * 
 * 
 */

#ifndef _TemplateElement_templates_h_
#define _TemplateElement_templates_h_

#include "TemplateElement.h"

AFEPACK_OPEN_NAMESPACE

template <int DIM>
bool operator==(const BasisFunctionIdentity<DIM>& i0, const BasisFunctionIdentity<DIM>& i1)
{
  if (i0.order != i1.order) return false;
  for (int i = 0;i < DIM;i ++)
    if (i0.alpha[i] != i1.alpha[i])
      return false;
  if (i0.flag != i1.flag ) return false;
  return true;
}

///////////////////////////////////////////////////////////////////////////

template <class value_type, int DIM>
  ShapeFunction<value_type,DIM>::ShapeFunction() :
  handle(NULL)
{}

template <class value_type, int DIM>
  ShapeFunction<value_type,DIM>::ShapeFunction(const ShapeFunction<value_type,DIM>& s) :
  handle(NULL),
  library_name(s.library_name),
  value_function_name(s.value_function_name),
  gradient_function_name(s.gradient_function_name)
{
  loadFunction();
}

template <class value_type, int DIM>
  ShapeFunction<value_type,DIM>::~ShapeFunction()
{
  unloadFunction();
}

template <class value_type, int DIM>
  ShapeFunction<value_type,DIM>& ShapeFunction<value_type,DIM>::operator=(const ShapeFunction<value_type,DIM>& s)
{
  if (&s != NULL) {
    library_name = s.library_name;
    value_function_name = s.value_function_name;
    gradient_function_name = s.gradient_function_name;
#ifdef QUADRATIC_ELEMENT_SUPPORT
    hesse_function_name = s.hesse_function_name;
#endif // QUADRATIC_ELEMENT_SUPPORT
  }
  return *this;
}

  template <class value_type, int DIM>
  void ShapeFunction<value_type,DIM>::loadFunction()
  {
    unloadFunction();
	
    std::string temp;
    if (library_path.length() == 0)
      temp = library_name;
    else
      temp = library_path + "/" + library_name;
    handle = AFEPackDLOpen(temp);
    if (handle == NULL) return;

    void * symbol = dlsym(handle, value_function_name.c_str());
    Assert(symbol, ExcLoadFunction(value_function_name.c_str(), library_name.c_str()));
    value_function = (void (*)(const double *, const double **, void *))symbol;

    symbol = dlsym(handle, gradient_function_name.c_str());
    Assert(symbol, ExcLoadFunction(gradient_function_name.c_str(), library_name.c_str()));
    gradient_function = (void (*)(const double *, const double **, void *))symbol;

#ifdef QUADRATIC_ELEMENT_SUPPORT
    symbol = dlsym(handle, hesse_function_name.c_str());
    Assert(symbol, ExcLoadFunction(hesse_function_name.c_str(), library_name.c_str()));
    hesse_function = (void (*)(const double *, const double **, void *))symbol;
#endif // QUADRATIC_ELEMENT_SUPPORT
  }

template <class value_type, int DIM>
    void ShapeFunction<value_type,DIM>::unloadFunction()
{
  if (handle != NULL) {
    dlclose(handle);
    handle = NULL;
  }
}

template <class value_type, int DIM>
  inline value_type ShapeFunction<value_type,DIM>::value(const Point<DIM>& p, 
                                                         const std::vector<Point<DIM> >& v) const
{
  value_type val;
  int k = v.size();
  const double * v1[k];
  for (int i = 0;i < k;i ++) v1[i] = v[i];
  (*value_function)(p, v1, (void *)(&val));
  return val;
}

template <class value_type, int DIM>
  inline std::vector<value_type> ShapeFunction<value_type,DIM>::gradient(const Point<DIM>& p, 
									 const std::vector<Point<DIM> >& v) const
{
  int k = v.size();
  const double * v1[k];
  for (int i = 0;i < k;i ++) v1[i] = v[i];
  std::vector<value_type> val(DIM);
  (*gradient_function)(p, v1, (void *)(&val[0]));
  return val;
}

template <class value_type, int DIM>
  inline std::vector<value_type> 
  ShapeFunction<value_type,DIM>::value(const std::vector<Point<DIM> >& p, 
                                       const std::vector<Point<DIM> >& v) const
{
  int k = v.size();
  int l = p.size();
  const double * v1[k];
  for (int i = 0;i < k;i ++) v1[i] = v[i];
  std::vector<value_type> val(l);
  for (int i = 0;i < l;i ++) {
    (*value_function)(p[i], v1, (void *)(&val[i]));
  }
  return val;
}

template <class value_type, int DIM>
  inline std::vector<std::vector<value_type> >
  ShapeFunction<value_type,DIM>::gradient(const std::vector<Point<DIM> >& p, 
                                          const std::vector<Point<DIM> >& v) const
{
  int k = v.size();
  int l = p.size();
  const double * v1[k];
  for (int i = 0;i < k;i ++) v1[i] = v[i];
  std::vector<std::vector<value_type> > val(l,std::vector<value_type>(DIM));
  for (int i = 0;i < l;i ++) {
    (*gradient_function)(p[i], v1, (void *)(&val[i][0]));
  }
  return val;
}

template <class value_type, int DIM>
  inline value_type ShapeFunction<value_type,DIM>::value(const Point<DIM>& p, const double ** v1) const
{
  value_type val;
  (*value_function)(p, v1, (void *)(&val));
  return val;
}

template <class value_type, int DIM>
  inline std::vector<value_type> 
  ShapeFunction<value_type,DIM>::gradient(const Point<DIM>& p, const double ** v1) const
{
  std::vector<value_type> val(DIM);
  (*gradient_function)(p, v1, (void *)(&val[0]));
  return val;
}

template <class value_type, int DIM>
  inline std::vector<value_type> 
  ShapeFunction<value_type,DIM>::value(const std::vector<Point<DIM> >& p, const double ** v1) const
{
  int l = p.size();
  std::vector<value_type> val(l);
  for (int i = 0;i < l;i ++) {
    (*value_function)(p[i], v1, (void *)(&val[i]));
  }
  return val;
}

template <class value_type, int DIM>
  inline std::vector<std::vector<value_type> >
  ShapeFunction<value_type,DIM>::gradient(const std::vector<Point<DIM> >& p, const double ** v1) const
{
  int l = p.size();
  std::vector<std::vector<value_type> > val(l, std::vector<value_type>(DIM));
  for (int i = 0;i < l;i ++) {
    (*gradient_function)(p[i], v1, (void *)(&val[i][0]));
  }
  return val;
}

#ifdef QUADRATIC_ELEMENT_SUPPORT

template <class value_type, int DIM>
  inline std::vector<std::vector<value_type> > 
  ShapeFunction<value_type,DIM>::hesse(const Point<DIM>& p, 
                                       const std::vector<Point<DIM> >& v) const
{
  value_type val[DIM*DIM];
  const double ** v1;
	
  int i, k;
  k = v.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) v1[i] = v[i];
  (*hesse_function)(p, v1, (void *)(&val));
  delete[] v1;
	
  std::vector<std::vector<value_type> > 
    result(DIM, 
	   std::vector<value_type>(DIM));
  for (i = 0;i < DIM;i ++) 
    for (k = 0;k < DIM;k ++)
      result[i][k] = val[i*DIM+k];
  return result;
}


template <class value_type, int DIM>
  inline std::vector<std::vector<std::vector<value_type> > >
  ShapeFunction<value_type,DIM>::hesse(const std::vector<Point<DIM> >& p, 
                                       const std::vector<Point<DIM> >& v) const
{
  value_type val[DIM*DIM];
  const double ** v1;
	
  int i, j, k, l;
  k = v.size();
  l = p.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) v1[i] = v[i];
  std::vector<std::vector<std::vector<value_type> > > 
    result(l,
	   std::vector<std::vector<value_type> >(DIM, 
						 std::vector<value_type>(DIM)));
  for (i = 0;i < l;i ++) {
    (*hesse_function)(p[i], v1, (void *)(&val));
    for (j = 0;j < DIM;j ++)
      for (m = 0;m < DIM;m ++)
	result[i][j][m] = val[j*DIM+m];
  }
  delete[] v1;	
  return result;
}

template <class value_type, int DIM>
  inline std::vector<std::vector<value_type> >
  ShapeFunction<value_type,DIM>::hesse(const Point<DIM>& p, 
                                       const double ** v1) const
{
  value_type val[DIM*DIM];
  (*hesse_function)(p, v1, (void *)(&val));
  std::vector<std::vector<value_type> > result(DIM, std::vector<value_type>(DIM));
  for (int i = 0;i < DIM;i ++)
    for (int j = 0;j < DIM;j ++)
      result[i][j] = val[i*DIM+j];
  return result;
}

template <class value_type, int DIM>
  inline std::vector<std::vector<std::vector<value_type> > >
  ShapeFunction<value_type,DIM>::hesse(const std::vector<Point<DIM> >& p, 
                                       const double ** v1) const
{
  value_type val[DIM*DIM];
	
  int i, j, l;
  l = p.size();
  std::vector<std::vector<std::vector<value_type> > > 
    result(l,
	   std::vector<std::vector<value_type> >(DIM, 
						 std::vector<value_type>(DIM)));
  for (i = 0;i < l;i ++) {
    (*hesse_function)(p[i], v1, (void *)(&val));
    for (j = 0;j < DIM;j ++)
      for (m = 0;m < DIM;m ++)
	result[i][j][m] = val[j*DIM+m];
  }
  return result;
}

#endif // QUADRATIC_ELEMENT_SUPPORT

///////////////////////////////////////////////////////////////////////////

template <class value_type, int DIM, int TDIM>
  BasisFunction<value_type,DIM,TDIM>::BasisFunction()
{}

template <class value_type, int DIM, int TDIM>
  BasisFunction<value_type,DIM,TDIM>::BasisFunction(const Point<TDIM>& p) :
  ip(p)
{}

template <class value_type, int DIM, int TDIM>
  BasisFunction<value_type,DIM,TDIM>::BasisFunction(const BasisFunction<value_type,DIM,TDIM>& b) :
  ip(b.ip), id(b.id)
{}

template <class value_type, int DIM, int TDIM>
  BasisFunction<value_type,DIM,TDIM>::~BasisFunction()
{}

template <class value_type, int DIM, int TDIM>
  BasisFunction<value_type,DIM,TDIM>& BasisFunction<value_type,DIM,TDIM>::operator=(const BasisFunction<value_type,DIM,TDIM>& b)
{
  if (&b != NULL) {
    ip = b.ip;
    id = b.id;
  }
  return *this;
}

  template <class value_type, int DIM, int TDIM>
  void BasisFunction<value_type,DIM,TDIM>::reinit(const Point<TDIM>& p)
  {
    ip = p;
  }

template <class value_type, int DIM, int TDIM>
    const Point<TDIM>& BasisFunction<value_type,DIM,TDIM>::interpPoint() const
{
  return ip;
}

template <class value_type, int DIM, int TDIM>
  Point<TDIM>& BasisFunction<value_type,DIM,TDIM>::interpPoint()
{
  return ip;
}

template <class value_type, int DIM, int TDIM>
  const typename BasisFunction<value_type,DIM,TDIM>::Identity& BasisFunction<value_type,DIM,TDIM>::identity() const
{
  return id;
}

template <class value_type, int DIM, int TDIM>
  typename BasisFunction<value_type,DIM,TDIM>::Identity& BasisFunction<value_type,DIM,TDIM>::identity()
{
  return id;
}

///////////////////////////////////////////////////////////////////////////

template <int DIM>
TemplateDOF<DIM>::TemplateDOF(TemplateGeometry<DIM>& g) : 
geometry(&g)
{
  if (geometry == NULL) return;
	
  int i;
  n_geometry_dof.resize(DIM+1);
  geometry_dof.resize(DIM+1);
  for (i = 0;i <= DIM;i ++) {
    n_geometry_dof[i].resize(geometry->n_geometry(i));
    geometry_dof[i].resize(geometry->n_geometry(i));
  }
}

template <int DIM>
TemplateDOF<DIM>::TemplateDOF(const TemplateDOF<DIM>& t) :
geometry(t.geometry)
{
  if (geometry == NULL) return;
	
  int i;
  n_geometry_dof.resize(DIM+1);
  geometry_dof.resize(DIM+1);
  for (i = 0;i <= DIM;i ++) {
    n_geometry_dof[i].resize(geometry->n_geometry(i));
    geometry_dof[i].resize(geometry->n_geometry(i));
  }
}

template <int DIM>
TemplateDOF<DIM>::~TemplateDOF()
{}

template <int DIM>
TemplateDOF<DIM>& TemplateDOF<DIM>::operator=(const TemplateDOF<DIM>& t)
{
  if (&t != NULL) {
    n_dof = t.n_dof;
    n_geometry_dof = t.n_geometry_dof;
    geometry_dof = t.geometry_dof;
    dof_index = t.dof_index;
    geometry = t.geometry;
  }
  return *this;
}

  template <int DIM>
  void TemplateDOF<DIM>::reinit(TemplateGeometry<DIM>& g)
  {
    geometry = &g;
    if (geometry == NULL) return;
	
    int i;
    n_geometry_dof.resize(DIM+1);
    geometry_dof.resize(DIM+1);
    for (i = 0;i <= DIM;i ++) {
      n_geometry_dof[i].resize(geometry->n_geometry(i));
      geometry_dof[i].resize(geometry->n_geometry(i));
    }
    dof_index.clear();
  }

template <int DIM>
void TemplateDOF<DIM>::readData(const std::string& s)
{
  std::string library_path = FindAFEPackLibraryFilePath(s);
  std::string filename(library_path + "/" + s);
  ExpandString(filename);
  library_path = filename.substr(0, filename.rfind('/'));

  filtering_istream is;
  OpenAFEPackLibraryFile(filename, is);
  is >> *this;
}

template <int DIM>
void TemplateDOF<DIM>::writeData(const std::string& s) const
{
  std::ofstream os(s.c_str());
  os << *this;
  os.close();
}

///////////////////////////////////////////////////////////////////////////

template <int TDIM, int DIM>
  CoordTransform<TDIM,DIM>::CoordTransform() :
  handle(NULL)
{}

template <int TDIM, int DIM>
  CoordTransform<TDIM,DIM>::CoordTransform(const CoordTransform<TDIM,DIM>& c) :
  handle(NULL),
  library_name(c.library_name),
  l2g_function_name(c.l2g_function_name),
  g2l_function_name(c.g2l_function_name),
  l2g_jacobian_function_name(c.l2g_jacobian_function_name),
  g2l_jacobian_function_name(c.g2l_jacobian_function_name)
{
  loadFunction();
}

template <int TDIM, int DIM>
  CoordTransform<TDIM,DIM>::~CoordTransform()
{
  unloadFunction();
}

template <int TDIM, int DIM>
  CoordTransform<TDIM,DIM>& CoordTransform<TDIM,DIM>::operator=(const CoordTransform<TDIM,DIM>& c)
{
  if (&c != NULL) {
    library_name = c.library_name;
    l2g_function_name = c.l2g_function_name;
    g2l_function_name = c.g2l_function_name;
    l2g_jacobian_function_name = c.l2g_jacobian_function_name;
    g2l_jacobian_function_name = c.g2l_jacobian_function_name;
  }
  loadFunction();
  return *this;
}

  template <int TDIM, int DIM>
  void CoordTransform<TDIM,DIM>::loadFunction()
  {
    unloadFunction();
	
    std::string temp;
    if (library_path.length() == 0)
      temp = library_name;
    else
      temp = library_path + "/" + library_name;
    handle = AFEPackDLOpen(temp);
    if (handle == NULL) return;

    void * symbol = dlsym(handle, l2g_function_name.c_str());
    Assert(symbol, ExcLoadFunction(l2g_function_name.c_str(), library_name.c_str()));
    l2g_function = (void (*)(const double *, const double **, const double **, double *))symbol;

    symbol = dlsym(handle, g2l_function_name.c_str());
    Assert(symbol, ExcLoadFunction(g2l_function_name.c_str(), library_name.c_str()));
    g2l_function = (void (*)(const double *, const double **, const double **, double *))symbol;
	
    symbol = dlsym(handle, l2g_jacobian_function_name.c_str());
    Assert(symbol, ExcLoadFunction(l2g_jacobian_function_name.c_str(), library_name.c_str()));
    l2g_jacobian_function = (double (*)(const double *, const double **, const double **))symbol;

    symbol = dlsym(handle, g2l_jacobian_function_name.c_str());
    Assert(symbol, ExcLoadFunction(g2l_jacobian_function_name.c_str(), library_name.c_str()));
    g2l_jacobian_function = (double (*)(const double *, const double **, const double **))symbol;
  }

template <int TDIM, int DIM>
    void CoordTransform<TDIM,DIM>::unloadFunction()
{
  if (handle != NULL) {
    dlclose(handle);
    handle = NULL;
  }
}

template <int TDIM, int DIM>
  Point<DIM> CoordTransform<TDIM,DIM>::local_to_global(
                                                       const Point<TDIM>& lp, 
                                                       const std::vector<Point<TDIM> >& lv,
                                                       const std::vector<Point<DIM> >& gv) const
{
  double val[DIM];
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  (*l2g_function)(lp, v1, v2, val);
  delete[] v1;
  delete[] v2;
  return Point<DIM>(val);
}

template <int TDIM, int DIM>
  Point<TDIM> CoordTransform<TDIM,DIM>::global_to_local(
							const Point<DIM>& gp,
							const std::vector<Point<TDIM> >& lv,
							const std::vector<Point<DIM> >& gv) const
{
  double val[TDIM];
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  (*g2l_function)(gp, v1, v2, val);
  delete[] v1;
  delete[] v2;
  return Point<TDIM>(val);
}

template <int TDIM, int DIM>
  double CoordTransform<TDIM,DIM>::local_to_global_jacobian(
							    const Point<TDIM>& lp, 
							    const std::vector<Point<TDIM> >& lv,
							    const std::vector<Point<DIM> >& gv) const
{
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  double val = (*l2g_jacobian_function)(lp, v1, v2);
  delete[] v1;
  delete[] v2;
  return val;
}

template <int TDIM, int DIM>
  double CoordTransform<TDIM,DIM>::global_to_local_jacobian(
							    const Point<DIM>& gp, 
							    const std::vector<Point<TDIM> >& lv,
							    const std::vector<Point<DIM> >& gv) const
{
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  double val = (*g2l_jacobian_function)(gp, v1, v2);
  delete[] v1;
  return val;
}

template <int TDIM, int DIM>
  std::vector<Point<DIM> > CoordTransform<TDIM,DIM>::local_to_global(
								     const std::vector<Point<TDIM> >& lp, 
								     const std::vector<Point<TDIM> >& lv,
								     const std::vector<Point<DIM> >& gv) const
{
  double val[DIM];
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  k = lp.size();
  std::vector<Point<DIM> > ret(k);
  for (i = 0;i < k;i ++) {
    (*l2g_function)(lp[i], v1, v2, val);
    ret[i] = Point<DIM>(val);
  }
  delete[] v1;
  delete[] v2;
  return ret;
}

template <int TDIM, int DIM>
  std::vector<Point<TDIM> > CoordTransform<TDIM,DIM>::global_to_local(
								      const std::vector<Point<DIM> >& gp,
								      const std::vector<Point<TDIM> >& lv,
								      const std::vector<Point<DIM> >& gv) const
{
  double val[TDIM];
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  k = gp.size();
  std::vector<Point<TDIM> > ret(k);
  for (i = 0;i < k;i ++) {
    (*g2l_function)(gp[i], v1, v2, val);
    ret[i] = Point<TDIM>(val);
  }
  delete[] v1;
  delete[] v2;
  return ret;
}

template <int TDIM, int DIM>
  std::vector<double> CoordTransform<TDIM,DIM>::local_to_global_jacobian(
									 const std::vector<Point<TDIM> >& lp, 
									 const std::vector<Point<TDIM> >& lv,
									 const std::vector<Point<DIM> >& gv) const
{
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  k = lp.size();
  std::vector<double> ret(k);
  for (i = 0;i < k;i ++)
    ret[i] = (*l2g_jacobian_function)(lp[i], v1, v2);
  delete[] v1;
  delete[] v2;
  return ret;
}

template <int TDIM, int DIM>
  std::vector<double> CoordTransform<TDIM,DIM>::global_to_local_jacobian(
									 const std::vector<Point<DIM> >& gp, 
									 const std::vector<Point<TDIM> >& lv,
									 const std::vector<Point<DIM> >& gv) const
{
  const double ** v1, ** v2;
	
  int i, k;
  k = lv.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  v2 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) {
    v1[i] = lv[i];
    v2[i] = gv[i];
  }
  k = gp.size();
  std::vector<double> ret(k);
  for (i = 0;i < k;i ++)
    ret[i] = (*g2l_jacobian_function)(gp[i], v1, v2);
  delete[] v1;
  delete[] v2;
  return ret;
}

template <int TDIM, int DIM>
  void CoordTransform<TDIM,DIM>::readData(const std::string& s)
{
  library_path = FindAFEPackLibraryFilePath(s);
  std::string filename(library_path + "/" + s);
  ExpandString(filename);
  library_path = filename.substr(0, filename.rfind('/'));

  filtering_istream is;
  OpenAFEPackLibraryFile(filename, is);
  is >> *this;
}

template <int TDIM, int DIM>
  void CoordTransform<TDIM,DIM>::writeData(const std::string& s) const
{
  std::ofstream os(s.c_str());
  os << *this;
}

///////////////////////////////////////////////////////////////////////////

template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>::BasisFunctionAdmin()
{}


template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>::BasisFunctionAdmin(const int& n) :
  std::vector<BasisFunction<value_type,DIM,TDIM> >(n)
{}

template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>::BasisFunctionAdmin(const int& n,
                                                              TemplateDOF<TDIM>& t) :
  std::vector<BasisFunction<value_type,DIM,TDIM> >(n),
  df(&t)
{}

template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>::BasisFunctionAdmin(TemplateDOF<TDIM>& t) :
  df(&t)
{}

template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>::BasisFunctionAdmin(const BasisFunctionAdmin<value_type,DIM,TDIM>& b) :
  df(b.df)
{}

template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>::~BasisFunctionAdmin()
{}

template <class value_type, int DIM, int TDIM>
  void BasisFunctionAdmin<value_type,DIM,TDIM>::reinit(TemplateDOF<TDIM>& t)
{
  df = &t;
}

template <class value_type, int DIM, int TDIM>
  BasisFunctionAdmin<value_type,DIM,TDIM>& BasisFunctionAdmin<value_type,DIM,TDIM>::operator=(const BasisFunctionAdmin<value_type,DIM,TDIM>& b)
{
  df = b.df;
  return *this;
}

  template <class value_type, int DIM, int TDIM>
  const TemplateDOF<TDIM>& BasisFunctionAdmin<value_type,DIM,TDIM>::dof() const
  {
    return *df;
  }

template <class value_type, int DIM, int TDIM>
    TemplateDOF<TDIM>& BasisFunctionAdmin<value_type,DIM,TDIM>::dof()
{
  return *df;
}

template <class value_type, int DIM, int TDIM>
  void BasisFunctionAdmin<value_type,DIM,TDIM>::readData(const std::string& s)
{
  library_path = FindAFEPackLibraryFilePath(s);
  std::string filename(library_path + "/" + s);
  ExpandString(filename);
  library_path = filename.substr(0, filename.rfind('/'));

  filtering_istream is;
  OpenAFEPackLibraryFile(filename, is);
  is >> *this;
}

template <class value_type, int DIM, int TDIM>
  void BasisFunctionAdmin<value_type,DIM,TDIM>::writeData(const std::string& s) const
{
  std::ofstream os(s.c_str());
  os << *this;
  os.close();
}

///////////////////////////////////////////////////////////////////////////

template <int DIM>
UnitOutNormal<DIM>::UnitOutNormal() :
handle(NULL)
{}

template <int DIM>
UnitOutNormal<DIM>::UnitOutNormal(const UnitOutNormal<DIM>& c) :
handle(NULL),
  library_name(c.library_name),
  function_name(c.function_name)
{
  loadFunction();
}

template <int DIM>
UnitOutNormal<DIM>::~UnitOutNormal()
{
  unloadFunction();
}

template <int DIM>
UnitOutNormal<DIM>& UnitOutNormal<DIM>::operator=(const UnitOutNormal<DIM>& c)
{
  if (&c != NULL) {
    library_name = c.library_name;
    function_name = c.function_name;
  }
  return *this;
}

  template <int DIM>
  void UnitOutNormal<DIM>::loadFunction()
  {
    unloadFunction();
	
    std::string temp;
    if (library_path.length() == 0)
      temp = library_name;
    else
      temp = library_path + "/" + library_name;
    handle = AFEPackDLOpen(temp);
    if (handle == NULL) return;

    void * symbol = dlsym(handle, function_name.c_str());
    Assert(symbol, ExcLoadFunction(function_name.c_str(), library_name.c_str()));
    function = (void (*)(const double *, const double **, int, double *))symbol;
  }

template <int DIM>
void UnitOutNormal<DIM>::unloadFunction()
{
  if (handle != NULL) {
    dlclose(handle);
    handle = NULL;
  }
}


template <int DIM>
std::vector<double> UnitOutNormal<DIM>::value(const Point<DIM>& p,
					      const std::vector<Point<DIM> >& v, const int& n) const
{
  double val[DIM];
  const double ** v1;
	
  int i, k;
  k = v.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) v1[i] = v[i];
  (*function)(p, v1, n, val);
  delete[] v1;
  return std::vector<double>(&val[0], &val[DIM]);
}

template <int DIM>
std::vector<std::vector<double> > UnitOutNormal<DIM>::value(const std::vector<Point<DIM> >& p,
							    const std::vector<Point<DIM> >& v, const int& n) const
{
  double val[DIM];
  const double ** v1;
	
  int i, k;
  k = v.size();
  typedef double * double_pointer;
  v1 = (const double **)new double_pointer[k];
  for (i = 0;i < k;i ++) v1[i] = v[i];
  k = p.size();
  std::vector<std::vector<double> > ret(k, std::vector<double>(DIM));
  for (i = 0;i < k;i ++) {
    (*function)(p[i], v1, n, val);
    std::copy(&val[0], &val[0]+DIM, ret[i].begin());
  }
  delete[] v1;
  return ret;
}

template <int DIM>
std::vector<double> UnitOutNormal<DIM>::value(const Point<DIM>& p,
					      const double ** v1, const int& n) const
{
  double val[DIM];
  (*function)(p, v1, n, val);
  std::vector<double> ret(DIM);
  return std::vector<double>(&val[0], &val[DIM]);
}

template <int DIM>
std::vector<std::vector<double> > UnitOutNormal<DIM>::value(const std::vector<Point<DIM> >& p,
							    const double ** v1, const int& n) const
{
  double val[DIM];
	
  int i, k;
  k = p.size();
  std::vector<std::vector<double> > ret(k, std::vector<double>(DIM, 0.0));
  for (i = 0;i < k;i ++) {
    (*function)(p[i], v1, n, val);
    std::copy(&val[0], &val[0]+DIM, ret[i].begin());
  }
  return ret;
}

template <int DIM>
void UnitOutNormal<DIM>::readData(const std::string& s)
{
  library_path = FindAFEPackLibraryFilePath(s);
  std::string filename(library_path + "/" + s);
  ExpandString(filename);
  library_path = filename.substr(0, filename.rfind('/'));

  filtering_istream is;
  OpenAFEPackLibraryFile(filename, is);
  is >> *this;
}

template <int DIM>
void UnitOutNormal<DIM>::writeData(const std::string& s) const
{
  std::ofstream os(s.c_str());
  os << *this;
}

///////////////////////////////////////////////////////////////////////////

template <class value_type, int DIM, int TDIM>
  TemplateElement<value_type,DIM,TDIM>::TemplateElement(
                                                        TemplateGeometry<TDIM>& g,
                                                        TemplateDOF<TDIM>& d,
                                                        CoordTransform<TDIM,DIM>& c,
                                                        BasisFunctionAdmin<value_type,DIM,TDIM>& b,
                                                        UnitOutNormal<DIM>& u) :
  geo(&g),
  df(&d),
  ct(&c),
  bf(&b),
  uon(&u)
{}

template <class value_type, int DIM, int TDIM>
  TemplateElement<value_type,DIM,TDIM>::TemplateElement(const TemplateElement<value_type,DIM,TDIM>& t) :
  geo(t.geo),
  df(t.df),
  ct(t.ct),
  bf(t.bf)
{}

template <class value_type, int DIM, int TDIM>
  TemplateElement<value_type,DIM,TDIM>::~TemplateElement()
{}

template <class value_type, int DIM, int TDIM>
  TemplateElement<value_type,DIM,TDIM>& 
  TemplateElement<value_type,DIM,TDIM>::operator=(const TemplateElement<value_type,DIM,TDIM>& t)
{
  if (&t != NULL) {
    geo = t.geo;
    df = t.df;
    ct = t.ct;
    bf = t.bf;
  }
  return *this;
}

  template <class value_type, int DIM, int TDIM>
  void TemplateElement<value_type,DIM,TDIM>::reinit(
                                                    TemplateGeometry<TDIM>& g,
                                                    TemplateDOF<TDIM>& d,
                                                    CoordTransform<TDIM,DIM>& c,
                                                    BasisFunctionAdmin<value_type,DIM,TDIM>& b,
                                                    UnitOutNormal<DIM>& u)
  {
    geo = &g;
    df = &d;
    ct = &c;
    bf = &b;
    uon = &u;
  }

template <class value_type, int DIM, int TDIM>
    const std::vector<Point<TDIM> >& TemplateElement<value_type,DIM,TDIM>::vertexArray() const
{
  return geo->vertexArray();
}

AFEPACK_CLOSE_NAMESPACE

#endif //_TemplateElement_templates_h_

/**
 * end of file
 * 
 */

