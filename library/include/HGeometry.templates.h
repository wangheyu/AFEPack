/**
 * @file   HGeometry.templates.h
 * @author Robert Lie
 * @date   Wed Jun 21 07:58:17 2003
 * 
 * @brief  
 * 
 * 
 */


#ifndef _HGeometry_templates_h_
#define _HGeometry_templates_h_

#include "HGeometry.h"

AFEPACK_OPEN_NAMESPACE

#define TEMPLATE template <int DIM, int DOW>
#define THIS HGeometry<DIM,DOW>

TEMPLATE
THIS::HGeometry() : HGeometryBase()
{
  index = 0;
  for (int i = 0;i < THIS::n_vertex;i ++)
    vertex[i] = NULL;
  for (int i = 0;i < THIS::n_boundary;i ++)
    boundary[i] = NULL;
  parent = NULL;
  for (int i = 0;i < THIS::n_child;i ++)
    child[i] = NULL;
  bmark = 0;
}

TEMPLATE
bool THIS::isRefined() const
{
  return (child[0] != NULL);
}

TEMPLATE
void THIS::refine()
{}

TEMPLATE
void THIS::checkIntegrity() const
{
  /**
   * check if every vertex of its boundary is a vertex of this geometry
   * 
   */
  int k;
  for (int i = 0;i < THIS::n_boundary;i ++) {
    const bound_t * b = boundary[i];
    for (int j = 0;j < bound_t::n_vertex;j ++) {
      for (k = 0;k < THIS::n_vertex;k ++) {
	if (b->vertex[j] == vertex[k]) {
	  break;
        }
      }
      Assert(k < THIS::n_vertex, ExcInternalError());
    }
  }
  /**
   * if this geometry is refined, its children will check integrity
   * recursively.
   * 
   */
  if (isRefined()) {
    for (int i = 0;i < THIS::n_child;i ++) {
      child[i]->checkIntegrity();
    }
  }
}

TEMPLATE
bool THIS::isIncludePoint(const Point<DOW>& p) const
{
  return true;
}

#undef THIS /// HGeometry<DIM,DOW>

/////////////////////////////////////////////////////////////////////

#define THIS HGeometryTree<DIM,DOW>

TEMPLATE
void THIS::checkIntegrity()
{
  RootIterator the_ele = beginRootElement();
  RootIterator end_ele = endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    the_ele->checkIntegrity();
  }
}

TEMPLATE
void THIS::clear()
{
  Tools tools;
  RootIterator 
    the_ele = beginRootElement(),
    end_ele = endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    tools.clearIndex(*the_ele);
  }

  the_ele = beginRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    tools.incrIndex(*the_ele);
  }

  the_ele = beginRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    tools.decrIndex(*the_ele);
  }

  root_element.clear();
}

TEMPLATE
void THIS::readEasyMesh(const std::string& filename)
{
  std::cerr << "THIS::readEasyMesh is only avaiable for 2-dimensional case."
	    << std::endl;
  abort();
}
	
TEMPLATE
void THIS::readMesh(const std::string& filename)
{
  std::cerr << "Reading in mesh data file " 
            << filename 
            << " as geometry tree root ..." 
            << std::endl;
  std::ifstream is(filename.c_str());

  /// 读入顶点的坐标
  u_int n_point;
  is >> n_point;
  std::cerr << "\t# points: " << n_point << std::endl;
  std::vector<Point<DOW> > point(n_point);
  for (u_int i = 0;i < n_point;i ++) is >> point[i];
  is >> n_point;
  std::vector<HGeometry<0,DOW> *> geo_0d(n_point, (HGeometry<0,DOW> *)NULL);
  for (u_int i = 0;i < n_point;i ++) {
    u_int j, k;
    is >> j; geo_0d[j] = new HGeometry<0,DOW>();
    is >> k >> k; *dynamic_cast<Point<DOW> *>(geo_0d[j]) = point[k];
    is >> k >> k >> geo_0d[j]->bmark;
  }
  point.clear();
	
#define GDIM 1
  u_int n_geo_1d;
  std::vector<HGeometry<GDIM,DOW>*> geo_1d;
  if (DIM >= 1) {/// 读入一维几何体的信息
    is >> n_geo_1d;
    std::cerr << "\t# 1D-geometry: " << n_geo_1d << std::endl;
    geo_1d.resize(n_geo_1d, NULL);
    for (u_int i = 0;i < n_geo_1d;i ++) {
      u_int j, k, l;
      is >> j >> k; geo_1d[j] = new HGeometry<GDIM,DOW>();
      for (k = 0;k < GDIM + 1;k ++) {
        is >> l; geo_1d[j]->vertex[k] = geo_0d[l];
      }
      is >> k;
      for (k = 0;k < GDIM + 1;k ++) {
        is >> l;
      }
      is >> geo_1d[j]->bmark;
    }
  }
#undef GDIM

#define GDIM 2
  u_int n_geo_2d;
  std::vector<HGeometry<GDIM,DOW>*> geo_2d;
  if (DIM >= 2) {/// 读入二维几何体的信息
    is >> n_geo_2d;
    std::cerr << "\t# 2D-geometry: " << n_geo_2d << std::endl;
    geo_2d.resize(n_geo_2d, NULL);
    for (u_int i = 0;i < n_geo_2d;i ++) {
      u_int j, k, l;
      is >> j >> k; geo_2d[j] = new HGeometry<GDIM,DOW>();
      for (k = 0;k < GDIM + 1;k ++) {
        is >> l; geo_2d[j]->vertex[k] = geo_0d[l];
      }
      is >> k;
      for (k = 0;k < GDIM + 1;k ++) {
        is >> l; geo_2d[j]->boundary[k] = geo_1d[l];
      }
      is >> geo_2d[j]->bmark;
    }
  }
#undef GDIM

#define GDIM 3
  u_int n_geo_3d;
  std::vector<HGeometry<GDIM,DOW>*> geo_3d;
  if (DIM >= 3) { /// 读入三维几何体的信息
    is >> n_geo_3d;
    std::cerr << "\t# 3D-geometry: " << n_geo_3d << std::endl;
    geo_3d.resize(n_geo_3d, NULL);
    for (u_int i = 0;i < n_geo_3d;i ++) {
      u_int j, k, l;
      is >> j >> k; geo_3d[j] = new HGeometry<GDIM,DOW>();
      for (k = 0;k < GDIM + 1;k ++) {
        is >> l; geo_3d[j]->vertex[k] = geo_0d[l];
      }
      is >> k;
      for (k = 0;k < GDIM + 1;k ++) {
        is >> l; geo_3d[j]->boundary[k] = geo_2d[l];
      }
      is >> geo_3d[j]->bmark;
    }
#undef GDIM
  }
  is.close();

  if (DIM == 1) {
    for (u_int i = 0;i < n_geo_1d;++ i) {
      this->rootElement().push_back((HGeometry<DIM,DOW> *)geo_1d[i]);
    }
  } else if (DIM == 2) {
    for (u_int i = 0;i < n_geo_2d;++ i) {
      this->rootElement().push_back((HGeometry<DIM,DOW> *)geo_2d[i]);
    }
  } else if (DIM == 3) {
    for (u_int i = 0;i < n_geo_3d;++ i) {
      this->rootElement().push_back((HGeometry<DIM,DOW> *)geo_3d[i]);
    }
  }
}

#undef THIS /// HGeometryTree<DIM,DOW>

///////////////////////////////////////////////////////////////////////////

#define THIS HElement<DIM,DOW>

TEMPLATE
THIS::HElement() 
: value(-1), parent(NULL), child(THIS::n_child, NULL) {}

TEMPLATE
THIS::HElement(const element_t& ele) :
index(ele.index), value(ele.value), 
  indicator(ele.indicator),
  h_element(ele.h_element) {}

TEMPLATE
THIS::~HElement()
{}

TEMPLATE
typename THIS::element_t& THIS::operator=(const element_t& ele) {
  index = ele.index;
  value = ele.value;
  indicator = ele.indicator;
  h_element = ele.h_element;
  return *this;
}

  TEMPLATE
  void THIS::refine()
  {}

TEMPLATE
bool THIS::isRefined() const
{
  return (child.size() > 0 && child[0] != NULL);
}

TEMPLATE
void THIS::checkIntegrity() const
{
  if (!isRefined()) return;
  for (int i = 0;i < THIS::n_child;i ++) {
    child[i]->checkIntegrity();
  }
}

TEMPLATE
bool THIS::isIncludePoint(const Point<DOW>& p) const
{
  return h_element->isIncludePoint(p);
}

#undef THIS /// HElement<DIM,DOW>

/////////////////////////////////////////////////////////////////////////

#define THIS IrregularMesh<DIM,DOW>

TEMPLATE
void THIS::semiregularize()
{
  if (! geometry_tree->lock()) { // 对几何遗传树进行加锁
    std::cerr << "The hierarchy geometry tree is locked, aborting ...";
    abort();
  }

  std::cerr << "Semiregularizing the mesh ...  " << std::flush;
  int round = 0;
  const char * timer = "-/|\\";
  bool flag;
  int n_element_refined = 0;
  prepareSemiregularize();
  do {
    std::cerr << "\b" << timer[round] << std::flush;
    round = (round + 1)%4;
    flag = false;
    semiregularizeHelper(flag, n_element_refined);
  } while (flag);
  std::cerr << "\bOK!\n"
	    << "\t" << n_element_refined 
	    << " elements refined in semiregularization."
	    << std::endl;
}

TEMPLATE
void THIS::semiregularizeHelper(bool& flag, 
                                int& n_element_refined)
{
  RootIterator the_ele = beginRootElement();
  RootIterator end_ele = endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    semiregularizeHelper(flag, *the_ele, n_element_refined);
  }
}

TEMPLATE
void THIS::semiregularizeHelper(bool& flag, 
                                element_t& element, 
                                int& n_element_refined)
{
  if (element.value == 0) { // if this element is a leaf element
    h_element_t& h_element = *(element.h_element);
    flag |= Tools().semiregularizeBoundary(h_element);
    if (! Tools().isSemiregular(h_element)) { // if this element is not semiregular
      flag = true;
      element.refine();
      element.value = 1;
      for (int i = 0;i < element_t::n_child;i ++) {
	element.child[i]->value = 0;
        Tools().setGeometryUsed(*(h_element.child[i]));
      }
      n_element_refined ++;
    }
  } else { // if its not a leaf element
    assert (element.value == 1);
    for (int i = 0;i < element_t::n_child;i ++) { // then check its children
      semiregularizeHelper(flag, *element.child[i], n_element_refined);
    }
  }
}

TEMPLATE
void THIS::prepareSemiregularize()
{
  Tools tools;
  RootIterator the_root_element = beginRootElement();
  RootIterator end_root_element = endRootElement();
  for (;the_root_element != end_root_element;++ the_root_element) {
    tools.setGeometryUnusedRecursively(*(the_root_element->h_element));
  }
  RootFirstIterator 
    the_ele = beginRootFirstElement(),
    end_ele = endRootFirstElement();
  for (;the_ele != end_ele;++ the_ele) {
    tools.setGeometryUsed(*(the_ele->h_element));
  }
}

TEMPLATE
THIS::IrregularMesh(tree_t& h_geometry_tree)
{
  setGeometryTree(&h_geometry_tree);
  regular_mesh = NULL;
}

TEMPLATE void 
THIS::reinit(tree_t& h_geometry_tree, 
             bool is_bare)
{
  if (is_bare) {
    geometry_tree = &h_geometry_tree;
  } else {
    clear();
    setGeometryTree(&h_geometry_tree);
  }
}

TEMPLATE
void THIS::globalRefine(unsigned int i)
{
  MeshAdaptor<DIM, DOW> mesh_adaptor(*this);
  mesh_adaptor.globalRefine(i);
}

TEMPLATE
void THIS::randomRefine(double percent)
{
  MeshAdaptor<DIM, DOW> mesh_adaptor(*this);
  mesh_adaptor.randomRefine(percent);
}

TEMPLATE
void THIS::setGeometryTree(tree_t * h_geometry_tree)
{
  std::cerr << "Constructing the root mesh from hierarchy geometry tree ..." << std::endl;
  int i, j;
  geometry_tree = h_geometry_tree;

  // construct the root elements at first	
  std::cerr << "\tconstructing elements ..." << std::flush;
  unsigned int n_root_element = geometry_tree->n_rootElement();
  std::vector<element_t *> h_element(n_root_element);
  typename tree_t::RootIterator 
    the_h_element = geometry_tree->beginRootElement(),
    end_h_element = geometry_tree->endRootElement();
  for (i = 0;the_h_element != end_h_element;++ the_h_element) 
    the_h_element->index = i ++;
  the_h_element = geometry_tree->beginRootElement();
  for (i = 0;the_h_element != end_h_element;++ the_h_element, ++ i) {
    element_t * element = new element_t();
    element->value = 0;
    element->h_element = &*the_h_element;
    root_element.push_back(element);
    h_element[i] = element;
  }
  std::cerr << " OK!" << std::endl;
}

TEMPLATE
THIS::~IrregularMesh()
{
  clear();
}

TEMPLATE
void THIS::clear()
{
  if (geometry_tree != NULL)
    geometry_tree = NULL;

  RootIterator the_ele = beginRootElement();
  RootIterator end_ele = endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    this->deleteTree(&*the_ele);
  }
  root_element.clear();

  if (regular_mesh != NULL) {
    delete regular_mesh;
    regular_mesh = NULL;
  }
}

TEMPLATE
void THIS::checkIntegrity()
{
  RootIterator the_ele = beginRootElement();
  RootIterator end_ele = endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    the_ele->checkIntegrity();
  }
}

TEMPLATE
void THIS::deleteTree(element_t * element)
{
  if (element->isRefined()) {
    for (int i = 0;i < element_t::n_child;i ++) {
      this->deleteTree(element->child[i]);
    }
  }
  delete element;
}

TEMPLATE
void THIS::copyTree(const element_t * src,
                    element_t * dst)
{
  dst->index = src->index;
  dst->value = src->value;
  dst->indicator = src->indicator;
  dst->h_element = src->h_element;
  if (src->isRefined()) {
    dst->refine();
    for (int i = 0;i < element_t::n_child;i ++) {
      this->copyTree(src->child[i], dst->child[i]);
    }
  }
}

TEMPLATE
void THIS::copyNonnegtiveSubtree(const element_t * src,
                                 element_t * dst)
{
  assert (src->value == 0 || src->value == 1);
  dst->value = src->value;
  dst->index = src->index;
  dst->h_element = src->h_element;
  if (src->value == 1) {
    dst->refine();
    for (int i = 0;i < element_t::n_child;i ++) {
      this->copyNonnegtiveSubtree(src->child[i], dst->child[i]);
    }
  }
}

TEMPLATE
void THIS::writeFormatted(const std::string& filename)
{
  std::ofstream os(filename.c_str());
  os << (*this);
  os.close();
}



TEMPLATE
THIS::IrregularMesh() :
geometry_tree(NULL), regular_mesh(NULL)
{}


TEMPLATE
THIS::IrregularMesh(const ir_mesh_t& mesh)
{
  if (mesh.geometry_tree != NULL) {
    setGeometryTree(mesh.geometry_tree);
    copyNonnegtiveSubtree(mesh);
  }
  regular_mesh = NULL;
}


	
TEMPLATE
void THIS::copyTree(const ir_mesh_t& mesh)
{
  ConstRootIterator the_ele = mesh.beginRootElement();
  ConstRootIterator end_ele = mesh.endRootElement();
  RootIterator the_ele_1 = beginRootElement();
  for (;the_ele != end_ele;++ the_ele, ++ the_ele_1) {
    this->copyTree(&*the_ele, &*the_ele_1);
  }
}



TEMPLATE
void THIS::copyNonnegtiveSubtree(const ir_mesh_t& mesh)
{
  ConstRootIterator the_ele = mesh.beginRootElement();
  ConstRootIterator end_ele = mesh.endRootElement();
  RootIterator the_ele_1 = beginRootElement();
  for (;the_ele != end_ele;++ the_ele, ++ the_ele_1) {
    this->copyNonnegtiveSubtree(&*the_ele, &*the_ele_1);
  }
}



TEMPLATE
typename THIS::ir_mesh_t& THIS::operator=(const ir_mesh_t& mesh)
{
  clear();
  setGeometryTree(mesh.geometry_tree);
  copyNonnegtiveSubtree(mesh);
  return *this;
}



TEMPLATE
typename THIS::RootFirstIterator THIS::beginRootFirstElement()
{
  typename std::list<element_t *>::iterator element = root_element.begin();
  return RootFirstIterator(this, element, *element);
}



TEMPLATE
typename THIS::RootFirstIterator THIS::endRootFirstElement()
{
  typename std::list<element_t *>::iterator element = root_element.end();
  return RootFirstIterator(this, element, NULL);
}



TEMPLATE
typename THIS::ActiveIterator THIS::beginActiveElement()
{
  RootFirstIterator it = this->beginRootFirstElement();
  while (it->value > 0) ++ it;
  return ActiveIterator(it);
}



TEMPLATE
typename THIS::ActiveIterator THIS::endActiveElement()
{
  typename std::list<element_t *>::iterator element = root_element.end();
  return ActiveIterator(this, element, NULL);
}

TEMPLATE
void THIS::refineElement(element_t& ele) {
  ele.refine();
  ele.value = 1;
  for (int k = 0;k < ele.n_child;k ++) {
    ele.child[k]->value = 0;
  }
}

TEMPLATE
void THIS::renumerateElement()
{
  Assert (regular_mesh != NULL, ExcInternalError());
  int i, j, k, l, m, n;

  std::cerr << "Renumerating element of the mesh ..." << std::endl;
  int n_ele = regular_mesh->n_geometry(DIM);
  std::list<int> element_index;
  std::vector<std::list<int>::iterator> element_index_iterator(n_ele);
  for (i = 0;i < n_ele;i ++) {
    element_index_iterator[i] = 
      element_index.insert(element_index.end(), 
                           i);
  }

  std::vector<std::list<std::pair<int,std::list<int>::iterator> > >
    element_to_node(regular_mesh->n_point());
  for (i = 0;i < n_ele;i ++) {
    GeometryBM& ele = regular_mesh->geometry(DIM,i);
    std::list<int>::iterator& the_it = element_index_iterator[i];
    for (j = 0;j < ele.n_vertex();j ++) {
      element_to_node[ele.vertex(j)].push_back
        (std::pair<int,std::list<int>::iterator>(i, the_it));
    }
  }

  std::vector<int> value(regular_mesh->n_geometry(DIM), 0);
  std::vector<int> new_index(regular_mesh->n_geometry(DIM));
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
        if (regular_mesh->geometry(DIM,idx).n_vertex() == n_used_vtx) {
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

    GeometryBM& ele = regular_mesh->geometry(DIM,j);
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

  std::vector<GeometryBM> tmp_ele(regular_mesh->geometry(DIM));
  std::vector<int> old_index(n_ele);
#ifdef __SERIALIZATION__
  std::vector<HGeometryBase*> hgp(regular_mesh->h_geometry_ptr[DIM]);
#endif
  for (i = 0;i < n_ele;i ++) {
    GeometryBM& ele = regular_mesh->geometry(DIM,i);
    ele = tmp_ele[new_index[i]];
    ele.index() = i;
    old_index[new_index[i]] = i;
#ifdef __SERIALIZATION__
    regular_mesh->h_geometry_ptr[DIM][i] = hgp[new_index[i]];
#endif
  }
  ActiveIterator 
    the_ele = beginActiveElement(),
    end_ele = endActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    the_ele->index = old_index[the_ele->index];
  }
  std::cerr << " OK!" << std::endl;
}

#undef THIS /// IrregularMesh<DIM,DOW>

/////////////////////////////////////////////////////////////////////////

TEMPLATE std::ostream& 
operator<<(std::ostream& os, const HGeometry<DIM,DOW>& geometry)
{
  for (int i = 0;i < geometry.n_boundary;i ++)
    os << *(geometry.boundary[i]);
  return os;
}

TEMPLATE
std::ostream& operator<<(std::ostream& os, const HElement<DIM, DOW>& element)
{
  if (element.value == 1) {
    for (int i = 0;i < element.n_child;i ++) {
      os << *(element.child[i]);
    }
  }
  else if (element.value == 0) {
    os << *(element.h_element);
  }
  else {
    Assert(false, ExcInternalError()); // what happeden? something must be wrong!
  }
  return os;
}

TEMPLATE
std::ostream& operator<<(std::ostream& os, IrregularMesh<DIM, DOW>& mesh)
{
  //	IrregularMesh<DIM, DOW>::RootElementIterator the_ele = mesh.beginRootElement();
  //	IrregularMesh<DIM, DOW>::RootElementIterator end_ele = mesh.endRootElement();
  //	for (;the_ele != end_ele;++ the_ele)
  //		os << *the_ele;
  ///////////////////////////////////////////////////////////////////////
  //	ActiveElementIterator<DIM, DOW> the_ele = mesh.beginActiveElement();
  //	ActiveElementIterator<DIM, DOW> end_ele = mesh.endActiveElement();
  //	for (;the_ele != end_ele;++ the_ele)
  //		os << *(the_ele->h_element);
  return os;
}

TEMPLATE
ElementIterator<DIM, DOW>::ElementIterator() :
  mesh(NULL), element(NULL) {}

/////////////////////////////////////////////////////////////////////////

#define THIS ElementIterator<DIM,DOW>

TEMPLATE
THIS::ElementIterator(IrregularMesh<DIM, DOW> * m, 
                      root_t& r, 
                      HElement<DIM, DOW> * e) :
mesh(m), root_element(r), element(e) {}

TEMPLATE
THIS::ElementIterator(const ElementIterator<DIM, DOW>& it) :
mesh(it.mesh), root_element(it.root_element), element(it.element) {}

TEMPLATE
THIS::~ElementIterator()
{}

TEMPLATE
ElementIterator<DIM, DOW>& THIS::operator=(const ElementIterator<DIM, DOW>& it)
{
  mesh = it.mesh;
  root_element = it.root_element;
  element = it.element;
  return *this;
}

#undef THIS /// ElementIterator<DIM,DOW>

/////////////////////////////////////////////////////////////////////////

TEMPLATE
bool operator==(const ElementIterator<DIM, DOW>& it0,
                ElementIterator<DIM, DOW>& it1)
{
  return (it0.mesh == it1.mesh &&
          it0.root_element == it1.root_element &&
          it0.element == it1.element);
}

TEMPLATE
bool operator!=(const ElementIterator<DIM, DOW>& it0,
                ElementIterator<DIM, DOW>& it1)
{
  return (it0.mesh != it1.mesh ||
	  it0.root_element != it1.root_element ||
	  it0.element != it1.element);
}

/////////////////////////////////////////////////////////////////////////

/**
 * 在支持树结构中存在哑几何体的情况下，遍历器将会跳过所有的哑几何体。
 */
TEMPLATE
RootFirstElementIterator<DIM, DOW>& 
  RootFirstElementIterator<DIM, DOW>::operator++()
{
  if (element == NULL) return *this; // do nothing
  else if (element->value == 1) { /// 非叶子节点会进入其后代
    element = element->child[0];
  } else {
    assert (element->value == 0);

    HElement<DIM, DOW> * child = element;
    HElement<DIM, DOW> * parent = child->parent;
    while (parent != NULL) {
      if (child != parent->child[n_child - 1]) break;
      child = parent;
      parent = child->parent;
    };
    if (parent == NULL) {
      ++ root_element;
      if (root_element == mesh->rootElement().end()) {
	element = NULL;
      } else {
	element = *(root_element);
      }
    } else {
      int i = 0;
      for (;child != parent->child[i];++ i) {} /// child 是第 i 个孩子
      element = parent->child[++ i];
    }
  }
  return *this;
}

TEMPLATE
ActiveElementIterator<DIM, DOW>& 
  ActiveElementIterator<DIM, DOW>::operator++()
{
  do {
    RootFirstElementIterator<DIM, DOW>::operator++();
    if (this->element == NULL) break;
  } while (this->element->value > 0);
  return *this;
}

////////////////////////////////////////////////////////////////////

#define THIS IrregularMeshPair<DIM, DOW>

TEMPLATE
THIS::IrregularMeshPair(ir_mesh_t& m0, ir_mesh_t& m1) :
mesh0(&m0), mesh1(&m1)
{
  Assert(mesh0->geometry_tree == mesh1->geometry_tree, ExcInternalError());
}

TEMPLATE
THIS::IrregularMeshPair(ir_mesh_t * m0, ir_mesh_t * m1) :
mesh0(m0), mesh1(m1)
{
  Assert(mesh0->geometry_tree == mesh1->geometry_tree, ExcInternalError());
}

TEMPLATE
THIS::~IrregularMeshPair()
{}

TEMPLATE
typename THIS::iterator THIS::beginActiveElementPair()
{
  RootFirstElementIterator<DIM, DOW> 
    it0 = mesh0->beginRootFirstElement(),
    it1 = mesh1->beginRootFirstElement();
	
  while ((*it0).value == 1 && (*it1).value == 1) {
    ++ it0; ++ it1;
  }
  typename iterator::State state;
  if ((*it0).value == 0 && (*it1).value == 0) {
    state = iterator::EQUAL;
  }
  else if ((*it0).value == 0) { // then (*it1).value == 1
    while ((*it1).value > 0) ++ it1;
    state = iterator::GREAT_THAN;
  }
  else { // then (*it0).value == 1 && (*it1).value == 0
    while ((*it0).value > 0) ++ it0;
    state = iterator::LESS_THAN;
  }
  return iterator(this, state, it0, it1);
}

TEMPLATE
typename THIS::iterator THIS::endActiveElementPair()
{
  RootFirstElementIterator<DIM, DOW> 
    it0 = mesh0->endRootFirstElement(),
    it1 = mesh1->endRootFirstElement();
  return iterator(this, iterator::EQUAL, it0, it1);
}

#undef THIS /// IrregularMeshPair<DIM, DOW>

////////////////////////////////////////////////////////////////////

#define THIS ActiveElementPairIterator<DIM, DOW>

TEMPLATE
THIS::ActiveElementPairIterator(const THIS& i) :
mesh_pair(i.mesh_pair), st(i.st),
  iterator0(i.iterator0), iterator1(i.iterator1) {}

TEMPLATE
THIS& THIS::operator=(const THIS& i)
{
  mesh_pair = i.mesh_pair;
  st = i.st;
  iterator0 = i.iterator0;
  iterator1 = i.iterator1;
	
  return *this;
}

TEMPLATE
THIS& THIS::operator++()
{
  if (iterator0.element == NULL && iterator1.element == NULL) {
    Assert(st == EQUAL, ExcInternalError());
  } else if (st == EQUAL) {
    ++ iterator0;
    ++ iterator1;
    if (iterator0.element == NULL || iterator1.element == NULL) {
      Assert (iterator0.element == NULL && iterator1.element == NULL, ExcInternalError());
      return *this;
    }
    while (iterator0->value > 0 && iterator1->value > 0) {
      ++ iterator0;
      ++ iterator1;
      if (iterator0.element == NULL || iterator1.element == NULL) {
        Assert(iterator0.element == NULL && iterator1.element == NULL, ExcInternalError());
        return *this;
      }
    };
    if (iterator0->value == 0 && iterator1->value == 0) {
      st = EQUAL;
    } else if (iterator0->value == 0) { // then iterator1->value > 0
      while (iterator1->value > 0) ++ iterator1;
      st = GREAT_THAN;
    } else { // iterator0->value > 0 && iterator1->value == 0
      while (iterator0->value > 0) ++ iterator0;
      st = LESS_THAN;
    }
  } else if (st == GREAT_THAN) {
    RootFirstElementIterator<DIM, DOW> next0 = iterator0;
    ++ next0;
    ++ iterator1;
    if (iterator1.element == NULL) {
      Assert(next0.element == NULL, ExcInternalError());
      iterator0 = next0;
      st = EQUAL;
    } else if (next0.element == NULL) {
      while (iterator1->value > 0) ++ iterator1;
    } else if (next0->h_element == iterator1->h_element) {
      iterator0 = next0;
      while (iterator0->value > 0 && iterator1->value > 0) {
        ++ iterator0;
        ++ iterator1;
      };
      if (iterator0->value == 0 && iterator1->value == 0) {
        st = EQUAL;
      } else if (iterator0->value == 0) { // then iterator1->value > 0
        while (iterator1->value > 0) ++ iterator1;
        st = GREAT_THAN;
      } else { // iterator0->value > 0 && iterator1->value == 0
        while (iterator0->value > 0) ++ iterator0;
        st = LESS_THAN;
      }
    } else {
      while (iterator1->value > 0) ++ iterator1;
    }
  } else { // st == LESS_THAN
    RootFirstElementIterator<DIM, DOW> next1 = iterator1;
    ++ next1;
    ++ iterator0;
    if (iterator0.element == NULL) {
      Assert(next1.element == NULL, ExcInternalError());
      iterator1 = next1;
      st = EQUAL;
    } else if (next1.element == NULL) {
      while (iterator0->value > 0) ++ iterator0;
    } else if (next1->h_element == iterator0->h_element) {
      iterator1 = next1;
      while (iterator0->value > 0 && iterator1->value > 0) {
        ++ iterator0;
        ++ iterator1;
      };
      if (iterator0->value == 0 && iterator1->value == 0) {
        st = EQUAL;
      } else if (iterator0->value == 0) { // then iterator1->value > 0
        while (iterator1->value > 0) ++ iterator1;
        st = GREAT_THAN;
      } else { // iterator0->value > 0 && iterator1->value == 0
        while (iterator0->value > 0) ++ iterator0;
        st = LESS_THAN;
      }
    } else {
      while (iterator0->value > 0) ++ iterator0;
    }
  }
  return *this;
}

TEMPLATE
bool operator==(const THIS& i0,
                THIS& i1)
{
  return (i0.mesh_pair == i1.mesh_pair &&
	  i0.st == i1.st &&
	  i0.iterator0 == i1.iterator0 &&
	  i0.iterator1 == i1.iterator1);
}

TEMPLATE
bool operator!=(const THIS& i0,
                THIS& i1)
{
  return (i0.mesh_pair != i1.mesh_pair ||
          i0.st != i1.st ||
          i0.iterator0 != i1.iterator0 ||
          i0.iterator1 != i1.iterator1);
}

#undef THIS /// ActiveElementPairIterator<DIM, DOW>

////////////////////////////////////////////////////////////////////////////////////////////////

#define THIS MeshAdaptor<DIM, DOW>

TEMPLATE
THIS::MeshAdaptor() :
from_mesh(NULL),
  to_mesh(NULL),
  ind(NULL),
  convergence_order(1),
  refine_step(1),
  refine_threshold(1.33333),
  coarse_threshold(0.75),
  _is_refine_only(false)
{}



TEMPLATE
THIS::MeshAdaptor(ir_mesh_t& f) :
from_mesh(&f), 
  to_mesh(&f),
  ind(NULL),
  convergence_order(1),
  refine_step(1),
  refine_threshold(1.33333),
  coarse_threshold(0.75),
  _is_refine_only(false) 
{}



TEMPLATE
THIS::MeshAdaptor(ir_mesh_t& f, ir_mesh_t& t) :
from_mesh(&f), 
  to_mesh(&t),
  ind(NULL),
  convergence_order(1),
  refine_step(1),
  refine_threshold(1.33333),
  coarse_threshold(0.75),
  _is_refine_only(false) 
{}



TEMPLATE
THIS::~MeshAdaptor()
{}



TEMPLATE
void THIS::globalRefine(unsigned int s)
{
  typedef typename ir_mesh_t::ActiveIterator iterator;
  std::cerr << "Global refine the mesh ..." << std::endl;
  for (unsigned int i = 0;i < s;i ++) {
    std::cerr << "\r\tround " << i+1 << " ..." << std::flush;
    iterator
      the_active_ele = to_mesh->beginActiveElement(),
      end_active_ele = to_mesh->endActiveElement();
    for (;the_active_ele != end_active_ele;) {
      iterator the_ele = the_active_ele;
      ++ the_active_ele;
      the_ele->refine();
      the_ele->value = 1;
      for (int k = 0;k < the_ele->n_child;k ++) {
	the_ele->child[k]->value = 0;
      }
    }
  }
  std::cerr << std::endl;
}



TEMPLATE
void THIS::randomRefine(double percent)
{
  typedef typename ir_mesh_t::ActiveIterator iterator;
  std::cerr << "Randomly refine the mesh ..." << std::endl;
  iterator
    the_active_ele = to_mesh->beginActiveElement(),
    end_active_ele = to_mesh->endActiveElement();
  for (;the_active_ele != end_active_ele;) {
    iterator the_ele = the_active_ele;
    ++ the_active_ele;
    int r = rand();
    if (100.0*r >= RAND_MAX*percent) continue;

    the_ele->refine();
    the_ele->value = 1;
    for (int k = 0;k < the_ele->n_child;k ++) {
      the_ele->child[k]->value = 0;
    }
  }
  std::cerr << std::endl;
}



TEMPLATE
void THIS::adapt()
{
  prepareToMesh();
  collectIndicator();
  implementAdaption();
}

TEMPLATE
void THIS::collectIndicator(HElement<DIM,DOW>& ele,
                            double convergenceCoefficient) 
{
  if (ele.value == 0) {
    ele.indicator = indicator(ele.index);
  } else {
    int n_chd = 0;
    ele.indicator = 0.0;
    for (int i = 0;i < ele.n_child;++ i) {
      HElement<DIM,DOW> * p_chd = ele.child[i];
      collectIndicator(*p_chd, convergenceCoefficient);
      ele.indicator += p_chd->indicator;
      n_chd += 1;
    }
    ele.indicator *= convergenceCoefficient*ele.n_child/n_chd;
  }
}

TEMPLATE
void THIS::collectIndicator()
{
  double c = pow(2.0, convergenceOrder());
  typename ir_mesh_t::RootIterator
    the_root = to_mesh->beginRootElement(),
    end_root = to_mesh->endRootElement();
  for (;the_root != end_root;++ the_root) {
    this->collectIndicator(*the_root, c);
  }
}



TEMPLATE
void THIS::prepareToMesh()
{
  if (from_mesh == to_mesh) return;
  *to_mesh = *from_mesh;
}

TEMPLATE
void THIS::adaptElement(HElement<DIM,DOW>& ele,
                        double convergenceCoefficient,
                        int refine_level) {
  if ((! is_refine_only()) && is_indicator_underflow(ele.indicator)) {
    ele.value = 0; /// 将此单元设为叶子节点
  } else if (ele.value == 0 && /// 叶子节点
             refine_level <= refine_step && /// 自适应步数还不到
             is_indicator_overflow(ele.indicator, convergenceCoefficient)) { /// 且指示子过大
    ele.refine(); /// 对单元进行加密
    ele.value = 1; /// 将单元设为非叶子节点
    for (int k = 0;k < ele.n_child;k ++) { /// 然后对孩子递归做自适应
      ele.child[k]->value = 0; /// 设置孩子为叶子节点
      /// 猜测孩子的指示子的值
      ele.child[k]->indicator = ele.indicator/convergenceCoefficient;
      /// 对孩子做自适应
      adaptElement(*ele.child[k], convergenceCoefficient, refine_level + 1);
    }
  } else if (ele.value == 1) { /// 对非叶子节点，对孩子做加密
    for (int i = 0;i < ele.n_child;i ++) {
      adaptElement(*ele.child[i], convergenceCoefficient, refine_level);
    }
  }
}

TEMPLATE
void THIS::implementAdaption()
{
  std::cerr << "Implementing mesh adaption ..." << std::flush;
  double convergenceCoefficient = pow(2.0, DIM + convergenceOrder());
  typename ir_mesh_t::RootIterator
    the_root = to_mesh->beginRootElement(),
    end_root = to_mesh->endRootElement();
  for (;the_root != end_root;++ the_root) {
    this->adaptElement(*the_root, convergenceCoefficient, 0);
  }
  std::cerr << " OK!" << std::endl;
}

#undef THIS /// MeshAdaptor<DIM,DOW>

////////////////////////////////////////////////////////////////////

#define THIS RegularMesh<DIM,DOW>

TEMPLATE
void THIS::renumerateElementHSFC(void (*f)(const double *, double *))
{
  std::cerr << "Renumerating element of the mesh using Hibert space filling curve ..." 
            << std::flush;
  int n_ele = this->n_geometry(DIM);
  std::vector<std::vector<double> > x(DOW, 
                                      std::vector<double>(n_ele, 0.0));
  for (int i = 0;i < n_ele;++ i) {
    GeometryBM& ele = this->geometry(DIM, i);
    int n_vtx = ele.n_vertex();
    Point<DOW> pnt;
    for (int j = 0;j < n_vtx;++ j) {
      int vtx_idx = this->geometry(0, ele.vertex(j)).vertex(0);
      pnt += this->point(vtx_idx);
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

  std::vector<GeometryBM> tmp_ele(this->geometry(DIM));
  std::vector<int> old_index(n_ele);
#ifdef __SERIALIZATION__
  std::vector<HGeometryBase*> hgp(h_geometry_ptr[DIM]);
#endif
  for (int i = 0;i < n_ele;i ++) {
    GeometryBM& ele = this->geometry(DIM,i);
    ele = tmp_ele[new_index[i]];
    ele.index() = i;
    old_index[new_index[i]] = i;
#ifdef __SERIALIZATION__
    h_geometry_ptr[DIM][i] = hgp[new_index[i]];
#endif
  }
  typename ir_mesh_t::ActiveIterator 
    the_ele = irregularMesh().beginActiveElement(),
    end_ele = irregularMesh().endActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    the_ele->index = old_index[the_ele->index];
  }
  std::cerr << " OK!" << std::endl;
}

TEMPLATE
void THIS::writeEasyMesh(const std::string& filename) const
{
  Assert (DIM == 2, ExcInternalError());
}



TEMPLATE
void THIS::writeTecplotData(const std::string& filename) const
{}

#undef TEMPLATE

AFEPACK_CLOSE_NAMESPACE

#endif // _HGeometry_templates_h_

/**
 * end of file
 * 
 */

