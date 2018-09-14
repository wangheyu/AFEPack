/**
 * @file   Geometry.templates.h
 * @author Robert Lie
 * @date   Fri Jun 23 15:42:15 2003
 * 
 * @brief  
 * 
 * 
 */

#ifndef _Geometry_templates_h_
#define _Geometry_templates_h_

#include <set>
#include "Geometry.h"

AFEPACK_OPEN_NAMESPACE

#define TEMPLATE template <int DIM>
#define THIS Point<DIM>

TEMPLATE
THIS::Point()
{
  for (int i = 0;i < DIM;i ++)
    x[i] = 0;
}

TEMPLATE
THIS::Point(const double * data)
{
  for (int i = 0;i < DIM;i ++)
    x[i] = data[i];
}

TEMPLATE
THIS::Point(const THIS& p)
{
  int i;
  for (i = 0;i < DIM;i ++)
    x[i] = p.x[i];
}

TEMPLATE
THIS::Point(double d, ...)
{
  x[0] = d;
  va_list ap;
  va_start(ap, d);
  for (int i = 1;i < DIM;i ++) {
    x[i] = va_arg(ap, double);
  }
  va_end(ap);
}

TEMPLATE
THIS::~Point()
{}

TEMPLATE
THIS::operator const double *() const
{
  return x;
}

TEMPLATE
THIS::operator double *()
{
  return x;
}

TEMPLATE
THIS& THIS::operator=(const THIS& p)
{
  int i;
  for (i = 0;i < DIM;i ++)
    x[i] = p.x[i];
  return *this;
}

  TEMPLATE
  const double& THIS::operator[](int i) const
  {
    return x[i];
  }

TEMPLATE
double& THIS::operator[](int i)
{
  return x[i];
}

TEMPLATE
double THIS::length() const
{
  double v = 0.;
  for (int i = 0;i < DIM;i ++)
    v += x[i]*x[i];
  return sqrt(v);
}

TEMPLATE
THIS& THIS::operator+=(const THIS& p)
{
  for (int i = 0;i < DIM;i ++)
    x[i] += p.x[i];
  return *this;
}

TEMPLATE
THIS& THIS::operator-=(const THIS& p)
{
  for (int i = 0;i < DIM;i ++)
    x[i] -= p.x[i];
  return *this;
}

TEMPLATE
THIS& THIS::operator*=(const double& s)
{
  for (int i = 0;i < DIM;i ++)
    x[i] *= s;
  return *this;
}

TEMPLATE
THIS& THIS::operator/=(const double& s)
{
  for (int i = 0;i < DIM;i ++)
    x[i] /= s;
  return *this;
}

TEMPLATE
THIS midpoint(const THIS& p1, const THIS& p2)
{
  double p[DIM];
  for (int i = 0;i < DIM;i ++)
    p[i] = (p1.x[i] + p2.x[i])/2.;
  return p;
}

TEMPLATE
double distance(const THIS& p1, const THIS& p2)
{
  double d = 0.0;
  for (int i = 0;i < DIM;i ++)
    d += (p1.x[i] - p2.x[i])*(p1.x[i] - p2.x[i]);
  return sqrt(d);
}


TEMPLATE
THIS barycenter(const std::vector<THIS >& p, const double * w)
{
  int i, j, k;
  double bc[DIM];
  k = p.size();
  if (w == NULL) {
    for (i = 0;i < DIM;i ++) {
      bc[i] = 0;
      for (j = 0;j < k;j ++)
	bc[i] += p[j][i];
      bc[i] /= k;
    }
  }
  else {
    double sw = 0;
    for (j = 0;j < k;j ++) sw += w[j];
    for (i = 0;i < DIM;i ++) {
      bc[i] = 0;
      for (j = 0;j < k;j ++)
	bc[i] += w[j]*p[j][i];
      bc[i] /= sw;
    }
  }
  return bc;
}

TEMPLATE
THIS operator+(const THIS& p1, const THIS& p2)
{
  double p[DIM];
  for (int i = 0;i < DIM;i ++)
    p[i] = p1.x[i] + p2.x[i];
  return p;
}

TEMPLATE
THIS operator-(const THIS& p1, const THIS& p2)
{
  double p[DIM];
  for (int i = 0;i < DIM;i ++)
    p[i] = p1.x[i] - p2.x[i];
  return p;
}

TEMPLATE
std::istream& operator>>(std::istream& is, THIS& p)
{
  for (int i = 0;i < DIM;i ++)
    is >> p.x[i];
  return is;
}

TEMPLATE
std::ostream& operator<<(std::ostream& os, const THIS& p)
{
  for (int i = 0;i < DIM;i ++)
    os << p.x[i] << "\t";
  return os;
}

#undef THIS
#undef TEMPLATE

///////////////////////////////////////////////////////////////////////////

#define TEMPLATE template <int DIM, int DOW>
#define THIS Mesh<DIM,DOW>

TEMPLATE
THIS::Mesh() :
geo(DIM+1)
{}

TEMPLATE
THIS::Mesh(const THIS& m) :
pnt(m.pnt),
  geo(m.geo)
{}

TEMPLATE
THIS::~Mesh()
{}

TEMPLATE
THIS& THIS::operator=(const THIS& m)
{
  if (&m != NULL) {
    pnt = m.pnt;
    geo = m.geo;
  }
  return *this;
}

  TEMPLATE
  unsigned int THIS::n_point() const
  {
    return pnt.size();
  }

TEMPLATE
unsigned int THIS::n_geometry(int n) const
{
  return geo[n].size();
}

TEMPLATE
const std::vector<Point<DOW> >& THIS::point() const
{
  return pnt;
}

TEMPLATE
std::vector<Point<DOW> >& THIS::point()
{
  return pnt;
}

TEMPLATE
const Point<DOW>& THIS::point(int i) const
{
  return pnt[i];
}

TEMPLATE
Point<DOW>& THIS::point(int i)
{
  return pnt[i];
}

TEMPLATE
const std::vector<std::vector<GeometryBM> >& THIS::geometry() const
{
  return geo;
}

TEMPLATE
std::vector<std::vector<GeometryBM> >& THIS::geometry()
{
  return geo;
}

TEMPLATE
const std::vector<GeometryBM>& THIS::geometry(int n) const
{
  return geo[n];
}

TEMPLATE
std::vector<GeometryBM>& THIS::geometry(int n)
{
  return geo[n];
}

TEMPLATE
const GeometryBM& THIS::geometry(int i, int j) const
{
  return geo[i][j];
}

TEMPLATE
GeometryBM& THIS::geometry(int i, int j)
{
  return geo[i][j];
}

TEMPLATE
typename THIS::bmark_t THIS::boundaryMark(int i, int j) const
{
  return geo[i][j].bm;
}

TEMPLATE
typename THIS::bmark_t& THIS::boundaryMark(int i, int j)
{
  return geo[i][j].bm;
}

TEMPLATE
void THIS::renumerateElement()
{
  int i, j, k, l, m, n;

  std::cerr << "Renumerating element of the mesh ..." << std::endl;
  int n_ele = n_geometry(DIM);
  std::list<int> element_index;
  std::vector<std::list<int>::iterator> element_index_iterator(n_ele);
  for (i = 0;i < n_ele;i ++) {
    element_index_iterator[i] = 
      element_index.insert(element_index.end(), 
                           i);
  }

  std::vector<std::list<std::pair<int,std::list<int>::iterator> > >
    element_to_node(n_point());
  for (i = 0;i < n_ele;i ++) {
    GeometryBM& ele = geometry(DIM,i);
    std::list<int>::iterator& the_it = element_index_iterator[i];
    for (j = 0;j < ele.n_vertex();j ++) {
      element_to_node[ele.vertex(j)].push_back
        (std::pair<int,std::list<int>::iterator>(i, the_it));
    }
  }

  std::vector<int> value(n_geometry(DIM), 0);
  std::vector<int> new_index(n_geometry(DIM));
  std::list<std::list<int>::iterator> involved_element_index;
  for (i = 0, n = -1;i < n_ele;i ++) {
    if (involved_element_index.empty()) {
      j = element_index.front();
      element_index.pop_front();
      value[j] ++;
    }
    else {
      std::list<std::list<int>::iterator>::iterator 
        it = involved_element_index.begin(),
        end = involved_element_index.end(),
        the_it = it;
      int n_the_used_vtx = value[**the_it];
      for (;it != end;++ it) {
        int& idx = **it;
        int& n_used_vtx = value[idx];
        if (geometry(DIM,idx).n_vertex() == n_used_vtx) {
          the_it = it; break;
        } else if (n_the_used_vtx < n_used_vtx) {
          the_it = it;
          n_the_used_vtx = n_used_vtx;
        }
      }
      j = **the_it;
      element_index.erase(*the_it);
      involved_element_index.erase(the_it);
    }

    GeometryBM& ele = geometry(DIM,j);
    for (k = 0;k < ele.n_vertex();k ++) {
      m = ele.vertex(k);
      std::list<std::pair<int,std::list<int>::iterator> >::iterator
        it = element_to_node[m].begin(),
        end = element_to_node[m].end();
      for (;it != end;++ it) {
	if (value[it->first] == 0) {
	  involved_element_index.push_back(it->second);
	}
	value[it->first] ++;
      }
    }
    new_index[i] = j;
    if (100*i/n_ele > n) {
      n = 100*i/n_ele;
      std::cerr << "\r" << n << "% OK!";
    }
  }

  std::vector<GeometryBM> tmp_ele(geometry(DIM));
  for (i = 0;i < n_ele;i ++) {
    GeometryBM& ele = geometry(DIM,i);
    ele = tmp_ele[new_index[i]];
    ele.index() = i;
  }
  std::cerr << " OK!" << std::endl;
}

TEMPLATE
void THIS::readData(const std::string& s)
{
  std::cerr << "Reading mesh data file " << s << " ..." << std::endl;
  std::ifstream is(s.c_str());
  is >> *this;
  is.close();
}

TEMPLATE
void THIS::writeData(const std::string& s) const
{
  std::cerr << "Writing mesh data file " << s << " ..." << std::endl;
  std::ofstream os(s.c_str());
  os << *this;
  os.close();
}

TEMPLATE
void THIS::readData1d(const std::string& s)
{
  Assert(DIM == 1, ExcMeshData("this method can only be used for 1d."));
  Assert(DOW == 1, ExcMeshData("this method can only be used for 1d."));
  int i, j;
  std::ifstream is(s.c_str());
  is >> i;
  std::vector<double> buffer(i);
  for (j = 0;j < i;j ++) is >> buffer[j];
  is.close();
  sort(buffer.begin(), buffer.end());
  point().resize(i);
  for (j = 0;j < i;j ++) point(j)[0] = buffer[j];
  geometry().resize(2);
  geometry(0).resize(i);
  boundaryMark(0, 0) = 1;
  boundaryMark(0, i-1) = 1;
  for (j = 0;j < i;j ++) {
    geometry(0, j).index() = j;
    geometry(0, j).vertex().resize(1,j);
    geometry(0, j).boundary().resize(1,j);
  }
  geometry(1).resize(i-1);
  for (j = 0;j < i-1;j ++) {
    geometry(1, j).index() = j;
    geometry(1, j).vertex().resize(2);
    geometry(1, j).vertex(0) = j;
    geometry(1, j).vertex(1) = j+1;
    geometry(1, j).boundary().resize(2);
    geometry(1, j).boundary(0) = j;
    geometry(1, j).boundary(1) = j+1;
  }
}


TEMPLATE
std::istream& operator>>(std::istream& is, THIS& m)
{
  int i, j, k, l;

  std::cerr << "\tReading points ... " << std::flush;	
  is >> i;
  m.pnt.resize(i);
  for (j = 0;j < i;j ++)
    is >> m.pnt[j];
  std::cerr << i << " OK!" << std::endl;
	
  for (i = 0;i <= DIM;i ++) {
    std:: cerr << "\tReading " << i << "-dim geometries ... " << std::flush;
    GeometryBM temp;
    std::vector<GeometryBM>& g = m.geo[i];
    is >> j;
    g.resize(j);
    for (k = 0;k < j;k ++) {
      is >> temp;
      l = temp.index();
      Assert(l < j, (typename THIS::ExcMeshData("geometry index error.")));
      g[l] = temp;
    }
    std::cerr << j << " OK!" << std::endl;
  }
  return is;
}

TEMPLATE
std::ostream& operator<<(std::ostream& os, const THIS& m)
{
  int i, j, k;
	
  os.precision(16); /// 设置输出流格式旗标
  os.setf(std::ios::scientific, std::ios::floatfield);
  std::cerr << "\tWriting points ... " << std::flush;
  i = m.pnt.size();
  os << i << "\n";
  for (j = 0;j < i;j ++)
    os << m.pnt[j] << "\n";
  std::cerr << i << " OK!" << std::endl;
	
  for (i = 0;i <= DIM;i ++) {
    std::cerr << "\tWriting " << i << "-dim geometries ... " << std::flush;
    const std::vector<GeometryBM>& g = m.geo[i];
    j = g.size();
    os << "\n" << j << "\n";
    for (k = 0;k < j;k ++) {
      os << g[k];
    }
    std::cerr << j << " OK!" << std::endl;
  }
  return os;
}

///////////////////////////////////////////////////////////////////////////

template <int DIM>
QuadratureInfo<DIM>::QuadratureInfo()
{}

template <int DIM>
QuadratureInfo<DIM>::QuadratureInfo(const QuadratureInfo<DIM>& q) :
alg_acc(q.alg_acc),
  pnt(q.pnt),
  wei(q.wei)
{}

template <int DIM>
QuadratureInfo<DIM>::~QuadratureInfo()
{}

template <int DIM>
QuadratureInfo<DIM>& QuadratureInfo<DIM>::operator=(const QuadratureInfo<DIM>& q)
{
  if (&q != NULL) {
    alg_acc = q.alg_acc;
    pnt = q.pnt;
    wei = q.wei;
  }
  return *this;
}

  template <int DIM>
  inline int QuadratureInfo<DIM>::n_quadraturePoint() const
  {
    return pnt.size();
  }

template <int DIM>
inline int QuadratureInfo<DIM>::algebricAccuracy() const
{
  return alg_acc;
}

template <int DIM>
inline int& QuadratureInfo<DIM>::algebricAccuracy()
{
  return alg_acc;
}

template <int DIM>
inline const std::vector<Point<DIM> >& QuadratureInfo<DIM>::quadraturePoint() const
{
  return pnt;
}

template <int DIM>
inline std::vector<Point<DIM> >& QuadratureInfo<DIM>::quadraturePoint()
{
  return pnt;
}

template <int DIM>
inline const Point<DIM>& QuadratureInfo<DIM>::quadraturePoint(int i) const
{
  return pnt[i];
}


template <int DIM>
inline Point<DIM>& QuadratureInfo<DIM>::quadraturePoint(int i)
{
  return pnt[i];
}

template <int DIM>
inline const std::vector<double>& QuadratureInfo<DIM>::weight() const
{
  return wei;
}

template <int DIM>
inline std::vector<double>& QuadratureInfo<DIM>::weight()
{
  return wei;
}

template <int DIM>
inline const double& QuadratureInfo<DIM>::weight(int i) const
{
  return wei[i];
}

template <int DIM>
inline double& QuadratureInfo<DIM>::weight(int i)
{
  return wei[i];
}

template <int DIM>
filtering_istream& operator>>(filtering_istream& is, QuadratureInfo<DIM>& q)
{
  int i, j;
  is >> q.alg_acc;
  is >> i;
  q.pnt.resize(i);
  q.wei.resize(i);
  for (j = 0;j < i;j ++) {
    is >> q.pnt[j];
    is >> q.wei[j];
  }
  return is;
}

template <int DIM>
std::ostream& operator<<(std::ostream& os, const QuadratureInfo<DIM>& q)
{
  int i, j;
  os << q.alg_acc << "\n";
  i = q.n_quadraturePoint();
  os << i << "\n";
  for (j = 0;j < i;j ++) {
    os << q.pnt[j];
    os << q.wei[j] << "\n";
  }
  return os;
}

///////////////////////////////////////////////////////////////////////////

template <int DIM>
QuadratureInfoAdmin<DIM>::QuadratureInfoAdmin()
{}

template <int DIM>
QuadratureInfoAdmin<DIM>::QuadratureInfoAdmin(const QuadratureInfoAdmin<DIM>& q) :
acc_tbl(q.acc_tbl)
{}

template <int DIM>
QuadratureInfoAdmin<DIM>::~QuadratureInfoAdmin()
{}

template <int DIM>
QuadratureInfoAdmin<DIM>& QuadratureInfoAdmin<DIM>::operator=(const QuadratureInfoAdmin<DIM>& q)
{
  acc_tbl = q.acc_tbl;
  return *this;
}

  template <int DIM>
  const QuadratureInfo<DIM>& QuadratureInfoAdmin<DIM>::find(int i) const
  {
    int j = i;
    if (acc_tbl[j] == -1) {
      for (;acc_tbl[j] == -1 && j < acc_tbl.size();j ++);
      if ((unsigned int)j > acc_tbl.size()) {
        std::cerr << "no such quadrature info, algebric accuracy: "
                  << i << std::endl;
        abort();
      }
      /*     std::cerr << "no such quadrature info, algebric accuracy: " */
      /* 	      << (unsigned int)i */
      /* 	      << "\na more accurate one with algebric accuracy " */
      /* 	      << j << " found!" << std::endl; */
    }
    return std::vector<QuadratureInfo<DIM> >::operator[](acc_tbl[j]);
  }

template <int DIM>
QuadratureInfo<DIM>& QuadratureInfoAdmin<DIM>::find(int i)
{
  int j = i, n = acc_tbl.size();
  if (j < n || acc_tbl[j] == -1) {
    for (;acc_tbl[j] == -1 && j < n;j ++);
  }
  if (j == n) {
    std::cerr << "no such quadrature info, algebric accuracy: "
	      << (unsigned int)i << std::endl;
    abort();
  }
  return std::vector<QuadratureInfo<DIM> >::operator[](acc_tbl[j]);
}

template <int DIM>
filtering_istream& operator>>(filtering_istream& is, QuadratureInfoAdmin<DIM>& q)
{
  int i, j, k;
	
  is >> i;
  q.resize(i);
  for (j = 0, k = -1;j < i;j ++) {
    is >> q[j];
    k = std::max(k, q[j].algebricAccuracy());
  }
  q.acc_tbl.resize(k+1, -1);
  for (j = 0;j < i;j ++) {
    q.acc_tbl[q[j].algebricAccuracy()] = j;
  }
  for (j = k;j >= 0;j --) {
    if (q.acc_tbl[j] != -1)
      i = q.acc_tbl[j];
    else
      q.acc_tbl[j] = i;
  }
  return is;
}

template <int DIM>
std::ostream& operator<<(std::ostream& os, const QuadratureInfoAdmin<DIM>& q)
{
  int i, j;
	
  i = q.size();
  os << i << "\n";
  for (j = 0;j < i;j ++)
    os << q[j] << "\n";
  return os;
}

///////////////////////////////////////////////////////////////////////////

template <int DIM>
TemplateGeometry<DIM>::TemplateGeometry() :
handle(NULL)
{}

template <int DIM>
TemplateGeometry<DIM>::TemplateGeometry(const TemplateGeometry<DIM>& t) :
handle(NULL),
  library_name(t.library_name),
  volume_function_name(t.volume_function_name),
  quad_info(t.quad_info)
{}

template <int DIM>
TemplateGeometry<DIM>::~TemplateGeometry()
{
  unloadFunction();
}

template <int DIM>
TemplateGeometry<DIM>& TemplateGeometry<DIM>::operator=(const TemplateGeometry<DIM>& t)
{
  if (&t != NULL) {
    library_name = t.library_name;
    volume_function_name = t.volume_function_name;
    quad_info = t.quad_info;
    loadFunction();
  }
  return *this;
}

  template <int DIM>
  void TemplateGeometry<DIM>::loadFunction()
  {
    unloadFunction();
    std::string temp;
    if (library_path.length() == 0)
      temp = library_name;
    else
      temp = library_path + "/" + library_name;
    handle = AFEPackDLOpen(temp);

    void * symbol = dlsym(handle, volume_function_name.c_str());
    Assert(symbol, ExcLoadFunction(volume_function_name.c_str(), library_name.c_str()));
    volume_function = (double (*)(const double **))symbol;
  }

template <int DIM>
void TemplateGeometry<DIM>::unloadFunction()
{
  if (handle != NULL) {
    dlclose(handle);
    handle = NULL;
  };
}

template <int DIM>
const QuadratureInfoAdmin<DIM>& TemplateGeometry<DIM>::quadratureInfo() const
{
  return quad_info;
}

template <int DIM>
QuadratureInfoAdmin<DIM>& TemplateGeometry<DIM>::quadratureInfo()
{
  return quad_info;
}

template <int DIM>
inline const QuadratureInfo<DIM>& TemplateGeometry<DIM>::findQuadratureInfo(int i) const
{
  return quad_info.find(i);
}

template <int DIM>
inline double TemplateGeometry<DIM>::volume() const
{
  int i, k;
  const double ** v;
  k = Mesh<DIM,DIM>::n_point();
  v = (const double **)new double[k];
  for (i = 0;i < k;i ++)
    v[i] = Mesh<DIM,DIM>::point(i);
  double val = (*volume_function)(v);
  delete[] v;
  return val;
}

template <int DIM>
void TemplateGeometry<DIM>::readData(const std::string& s)
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
void TemplateGeometry<DIM>::writeData(const std::string& s) const
{
  std::ofstream os(s.c_str());
  os << *this;
  os.close();
}

template <int DIM>
inline const std::vector<Point<DIM> >& TemplateGeometry<DIM>::vertexArray() const
{
  return Mesh<DIM,DIM>::point();
}

template <int DIM>
filtering_istream& operator>>(filtering_istream& is, TemplateGeometry<DIM>& t)
{
  int i, j, k, l, m;

  is >> t.library_name
     >> t.volume_function_name;
  t.loadFunction();
	
  is >> i;
  t.point().resize(i);
  for (j = 0;j < i;j ++) {
    is >> t.point(j);
  }
	
  for (i = 0;i <= DIM;i ++) {
    Geometry temp;
    std::vector<GeometryBM>& g = t.geometry(i);
    is >> j;
    g.resize(j);
    for (k = 0;k < j;k ++) {
      is >> temp;
      // check if the vertices of the geometry are sorted
      const std::vector<int>& v = temp.vertex();
      m = v.size() - 1;
      for (l = 0;l < m;l ++)
	Assert(v[l] < v[l+1], (typename TemplateGeometry<DIM>::ExcTemplateGeometryData("geometry vertices not sorted.")));
      l = temp.index();
      Assert(l < j, (typename TemplateGeometry<DIM>::ExcTemplateGeometryData("geometry index error.")));
      dynamic_cast<Geometry&>(g[l]) = temp;
    }
  }
	
  is >> t.quad_info;
  return is;
}

template <int DIM>
std::ostream& operator<<(std::ostream& os, const TemplateGeometry<DIM>& t)
{
  int i, j, k;
	
  os << t.library_name << "\n\t"
     << t.volume_function_name << "\n";

  os.setf(std::ios::scientific);
  i = t.point().size();
  os << i << "\n";
  for (j = 0;j < i;j ++) {
    os << t.point(j) << "\n";
  }
  os << "\n";
	
  for (i = 0;i <= DIM;i ++) {
    const std::vector<GeometryBM>& g = t.geometry(i);
    j = g.size();
    os << j << "\n";
    for (k = 0;k < j;k ++) {
      os << dynamic_cast<const Geometry&>(g[k]) << "\n";
    }
    os << "\n";
  }
  os << "\n";

  os << t.quad_info << "\n";
  return os;
}

#undef THIS
#undef TEMPLATE

///////////////////////////////////////////////////////////////////////////

template <int DIM, int DOW> void 
  SimplestMesh<DIM,DOW>::generateMesh(Mesh<DIM,DOW>& m)
{
  int i, j, k, l, n, p;
  std::cerr << "Generate mesh structure from the simplest mesh ..." << std::endl;
  
  /**
   * 为每个单元做单元片，分为两个步骤：为每个顶点做单元片，然后将每个
   * 单元的顶点上的单元片合并起来。
   */
  n = n_element();
  std::vector<std::vector<u_int> > pnt_patch(n_point());
  for (i = 0;i < n;++ i) {
    for (j = 0;j < element(i).vertex.size();++ j) {
      pnt_patch[elementVertex(i, j)].push_back(i);
    }
  }
  std::vector<std::set<u_int, std::less<u_int> > > ele_patch(n);
  for (i = 0;i < n;++ i) {
    for (j = 0;j < element(i).vertex.size();++ j) {
      ele_patch[i].insert(pnt_patch[elementVertex(i, j)].begin(),
                          pnt_patch[elementVertex(i, j)].end());
    }
  }
  pnt_patch.clear();

  std::vector<std::vector<std::vector<int> > >
    ele_geo(n, std::vector<std::vector<int> >(DIM + 1));

  GeometryBM g;
  for (i = 0;i <= DIM;i ++) m.geometry(i).clear();
  for (i = 0, p = -1;i < n;i ++) {
    std::vector<std::vector<int> >& geo_img = ele_geo[i];

    const TemplateGeometry<DIM>& t_geo = (*tg)[element(i).template_element];
    geo_img[0].resize(t_geo.n_point(), -1);
    g.vertex().resize(1);
    g.boundary().resize(1);
    for (j = 0;j < t_geo.n_point();j ++) {
      g.vertex(0) = elementVertex(i, j);
      g.boundary(0) = elementVertex(i, j);

      bool is_found = false;
      int geo_idx;
      std::set<u_int, std::less<u_int> >::iterator
        the_ele = ele_patch[i].begin(), end_ele = ele_patch[i].end();
      for (;the_ele != end_ele;++ the_ele) {
        u_int ele_idx = *the_ele;
        if (ele_idx >= i) continue;
        for (l = 0;l < ele_geo[ele_idx][0].size();++ l) {
          geo_idx = ele_geo[ele_idx][0][l];
          if (geo_idx >=0 && m.geometry(0, geo_idx).vertex(0) == g.vertex(0)) {
            is_found = true;
            break;
          }
        }
        if (is_found) break;
      }
      if (!is_found) {
        geo_idx = m.n_geometry(0);
        g.index() = geo_idx;
        m.geometry(0).push_back(g);
      }
      geo_img[0][j] = geo_idx;
    }
    for (j = 1;j <= DIM;j ++) {
      geo_img[j].resize(t_geo.n_geometry(j));
      for (k = 0;k < t_geo.n_geometry(j);k ++) {
	g.vertex().resize(t_geo.geometry(j,k).n_vertex());
	g.boundary().resize(t_geo.geometry(j,k).n_boundary());
	for (l = 0;l < g.n_vertex();l ++)
	  g.vertex(l) = geo_img[0][t_geo.geometry(j,k).vertex(l)];
	for (l = 0;l < g.n_boundary();l ++)
	  g.boundary(l) = geo_img[j-1][t_geo.geometry(j,k).boundary(l)];

        bool is_found = false;
        int geo_idx;
        std::set<u_int, std::less<u_int> >::iterator
          the_ele = ele_patch[i].begin(), end_ele = ele_patch[i].end();
        for (;the_ele != end_ele;++ the_ele) {
          u_int ele_idx = *the_ele;
          if (ele_idx >= i) continue;
          for (l = 0;l < ele_geo[ele_idx][j].size();++ l) {
            geo_idx = ele_geo[ele_idx][j][l];
            if (geo_idx >= 0 && isSame(m.geometry(j, geo_idx), g)) {
              is_found = true;
              break;
            }
          }
          if (is_found) break;
        }

        if (!is_found) {
	  geo_idx = m.n_geometry(j);
	  g.index() = geo_idx;
	  m.geometry(j).push_back(g);
	}
        geo_img[j][k] = geo_idx;
      }
    }
    if (100*i/n > p) {
      p = 100*i/n;
      std::cerr << "\r" << p << "% OK!" << std::flush;
    }
  }
  std::cerr << "\r";

  /// 构造逆映射指标	
  int n_vtx = m.n_geometry(0);
  std::vector<int> index(pnt.size());
  for (k = 0;k < m.n_geometry(0);k ++) {
    index[m.geometry(0, k).vertex(0)] = k;
  }
  /// 对所有各维几何体的顶点进行校正
  for (j = 1;j <= DIM;j ++) {
    for (k = 0;k < m.n_geometry(j);k ++) {
      GeometryBM& geo = m.geometry(j, k);
      for (l = 0;l < geo.n_vertex();l ++) {
	n = geo.vertex(l);
	geo.vertex(l) = index[m.geometry(0,n).vertex(0)];
      }
    }
  }
  /// 对线的边界进行校正
  for (k = 0;k < m.n_geometry(1);k ++) {
    GeometryBM& geo = m.geometry(1, k);
    for (l = 0;l < geo.n_boundary();l ++) {
      n = geo.boundary(l);
      geo.boundary(l) = index[m.geometry(0,n).vertex(0)];
    }
  }
  /// 对顶点进行校正
  m.point().resize(n_vtx);
  for (k = 0;k < m.n_geometry(0);k ++) {
    m.point(k) = pnt[m.geometry(0, k).vertex(0)];
    m.geometry(0, k).vertex(0) = k;
    m.geometry(0, k).boundary(0) = k;
  }
}

template <int DIM, int DOW>
void Mesh<DIM,DOW>::renumerateElementHSFC(void (*f)(const double *, double*))
{
  std::cerr << "Renumerating element of the mesh using Hibert space curve filling ..." 
            << std::flush;
  int n_ele = n_geometry(DIM);
  std::vector<std::vector<double> > x(DOW, 
                                      std::vector<double>(n_ele, 0.0));
  for (int i = 0;i < n_ele;++ i) {
    GeometryBM& ele = geometry(DIM, i);

    Point<DOW> pnt;
    int n_vtx = ele.n_vertex();
    for (int j = 0;j < n_vtx;++ j) {
      int vtx_idx = geometry(0, ele.vertex(j)).vertex(0);
      pnt += point(vtx_idx);
    }
    pnt /= n_vtx;
    if (f != NULL) {
      Point<DOW> pnt1(pnt);
      f(&(pnt1[0]), &(pnt[0]));
    }
    for (int k = 0;k < DOW;++ k) {
        x[k][i] = pnt[k];
    }
  }
  std::vector<int> new_index(n_ele);
  switch (DOW) {
  case 2:
    hsfc_renumerate(n_ele, &x[0][0], &x[1][0], &new_index[0]);
    break;
  case 3:
    hsfc_renumerate(n_ele, &x[0][0], &x[1][0], &x[2][0], &new_index[0]);
    break;
  }

  std::vector<GeometryBM> tmp_ele(geometry(DIM));
  for (int i = 0;i < n_ele;i ++) {
    GeometryBM& ele = geometry(DIM,i);
    ele = tmp_ele[new_index[i]];
    ele.index() = i;
  }
  std::cerr << " OK!" << std::endl;
}

AFEPACK_CLOSE_NAMESPACE

#endif // Geometry_templates_h_

//
// end of file
///////////////////////////////////////////////////////////////////////////
