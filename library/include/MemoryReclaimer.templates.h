/**
 * @file   MemoryReclaimer.templates.h
 * @author Robert Lie
 * @date   Sun Apr 29 15:50:15 2007
 * 
 * @brief  内存回收的实现
 * 
 * 
 */

#define TEMPLATE template <int DIM, int DOW>
#define THIS MemoryReclaimer<DIM,DOW>

/**
 * 几何体 g 使用 lab 进行标识。
 */
TEMPLATE
template <int DIM1> inline void 
THIS::labelHGeometry(HGeometry<DIM1,DOW>& g, int lab)
{
  for (int i = 0;i < g.n_vertex;++ i) {
    labelHGeometry(*(g.vertex[i]), lab);
  }
  for (int i = 0;i < g.n_boundary;++ i) {
    labelHGeometry(*(g.boundary[i]), lab);
  }
  g.index = lab;
}

/**
 * 对几何体 g 进行递归标识：如果它有后代，先标识其后代，然后标识自己。
 * 
 */
TEMPLATE
template <int DIM1> inline void 
THIS::labelHGeometryRecursively(HGeometry<DIM1,DOW>& g, int lab)
{
  for (int i = 0;i < g.n_boundary;++ i) {
    labelHGeometryRecursively(*(g.boundary[i]), lab);
  }
  if (g.isRefined()) {
    for (int i = 0;i < g.n_child;++ i) {
      labelHGeometryRecursively(*(g.child[i]), lab);
    }
  }

  labelHGeometry(g, lab);
}

/**
 * 对点几何体 g 使用 lab 进行标识。
 * 
 */
TEMPLATE inline void 
THIS::labelHGeometry(HGeometry<0,DOW>& g, int lab)
{
  g.index = lab;
}

TEMPLATE inline void 
THIS::labelHGeometryRecursively(HGeometry<0,DOW>& g, int lab)
{
  g.index = lab;
}

/**
 * 对IrregularMesh中未使用的内存做回收。我们只是简单地将每个叶子节点的
 * 所有后代删除即可。
 * 
 */
TEMPLATE void 
THIS::reclaimIrregularMesh(typename THIS::ir_mesh_t& m) 
{
  ActiveElementIterator<DIM, DOW> 
    the_ele = m.beginActiveElement(),
    end_ele = m.endActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    if (the_ele->isRefined()) {
      for (int i = 0;i < the_ele->n_child;++ i) {
        m.deleteTree(the_ele->child[i]);
        the_ele->child[i] = NULL;
      }
    }
  }
}

/**
 * 对整个HGeometryTree中的所有几何体都使用 -1 进行标识。
 * 
 */
TEMPLATE
void THIS::initialTreeLabel()
{
  typename HGeometryTree<DIM,DOW>::RootIterator 
    the_ele = h_tree->beginRootElement(),
    end_ele = h_tree->endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    labelHGeometryRecursively(*the_ele, -1);
  }
}

/**
 * 对于IrregularMesh中索引到的单元及其组件几何体使用 1 进行标识。
 * 
 */
TEMPLATE
void THIS::labelIrregularMesh(typename THIS::ir_mesh_t& m)
{
  RootFirstElementIterator<DIM, DOW>
    the_ele = m.beginRootFirstElement(),
    end_ele = m.endRootFirstElement();
  for (;the_ele != end_ele;++ the_ele) {
    labelHGeometry(*(the_ele->h_element), 1);
  }

  typename THIS::ir_mesh_t::mesh_t& mesh = m.regularMesh();
#define LABEL_USED_HGEOMETRY(D)                                         \
  if (dim >= D) {                                                       \
    int n_geo = mesh.h_geometry()[D].size();                            \
    for (int i = 0;i < n_geo;++ i) {                                    \
      HGeometry<D,dow> * p_geo = mesh.template h_geometry<D>(i);        \
      labelHGeometry(*p_geo, 1);                                        \
    }                                                                   \
  }

  LABEL_USED_HGEOMETRY(0);
  LABEL_USED_HGEOMETRY(1);
  LABEL_USED_HGEOMETRY(2);
  LABEL_USED_HGEOMETRY(3);

#undef LABEL_USED_HGEOMETRY
}

/**
 * 对几何体 g 的自身、组件和后代进行重标识：如果其原标识为 1，则我们什
 * 么都不做，返回 1；如果其原标识为 -1，说明这个几何体是第一次被索引到，
 * 我们将其标识修改为 -2，并返回 -1；如果其原标识为 -2，说明这个几和体
 * 已经不是第一次被索引到了，我们将这个索引拆除，然后返回 -2。
 * 
 */
TEMPLATE
template <int DIM1> inline int 
THIS::relabelHGeometryRecursively(HGeometry<DIM1,DOW>& g)
{
  for (int i = 0;i < g.n_vertex;++ i) {
    if (g.vertex[i] == NULL) continue;
    if (relabelHGeometryRecursively(*(g.vertex[i])) == -2) {
      g.vertex[i] = NULL;
    }
  }
  for (int i = 0;i < g.n_boundary;++ i) {
    if (g.boundary[i] == NULL) continue;
    if (relabelHGeometryRecursively(*(g.boundary[i])) == -2) {
      g.boundary[i] = NULL;
    }
  }
  if (g.isRefined()) {
    for (int i = 0;i < g.n_child;++ i) {
      if (g.child[i] == NULL) continue;
      if (relabelHGeometryRecursively(*(g.child[i])) == -2) {
        g.child[i] = NULL;
      }
    }
  }
  Assert ((g.index == -1) || 
          (g.index ==  1) ||
          (g.index == -2), ExcInternalError());
  if (g.index == -1) {
    g.index = -2;
    return -1;
  }
  else return g.index;
}


/**
 * 对点几何体 g 进行重新标识。
 * 
 */
TEMPLATE inline int 
THIS::relabelHGeometryRecursively(HGeometry<0,DOW>& g)
{
  Assert ((g.index == -1) || 
          (g.index ==  1) ||
          (g.index == -2), ExcInternalError());
  if (g.index == -1) {
    g.index = -2;
    return -1;
  }
  else return g.index;
}

/**
 * 对几何体 g 自身、组件已经其后代进行回收。
 * 
 */
TEMPLATE
template <int DIM1> inline int 
THIS::reclaimHGeometryRecursively(HGeometry<DIM1,DOW>& g)
{
  for (int i = 0;i < g.n_vertex;++ i) {
    if (g.vertex[i] != NULL) {
      if (reclaimHGeometryRecursively(*(g.vertex[i])) == -1) {
        g.vertex[i] = NULL;
      }
    }
  }
  for (int i = 0;i < g.n_boundary;++ i) {
    if (g.boundary[i] != NULL) {
      if (reclaimHGeometryRecursively(*(g.boundary[i])) == -1) {
        g.boundary[i] = NULL;
      }
    }
  }
  for (int i = 0;i < g.n_child;++ i) {
    if (g.child[i] == NULL) continue;
    if (reclaimHGeometryRecursively(*(g.child[i])) == -1) {
      g.child[i] = NULL;
    }
  }
  Assert ((g.index == 1 || g.index == -2), ExcInternalError());
  if (g.index == -2) {
    this->reclaimHGeometry(&g, DIM1);
    return -1;
  }
  else {
    return 1;
  }
}


/**
 * 对点几何体 g 进行回收。
 * 
 */
TEMPLATE inline int 
THIS::reclaimHGeometryRecursively(HGeometry<0,DOW>& g)
{
  Assert ((g.index == 1 || g.index == -2), ExcInternalError());
  if (g.index == -2) {
    this->reclaimHGeometry(&g, 0);
    return -1;
  }
  else {
    return 1;
  }
}


/**
 * 回收 HGeometryTree 中不用的内存。
 * 
 */
TEMPLATE
void THIS::reclaimTreeMemory()
{
  typename HGeometryTree<DIM,DOW>::RootIterator
    the_ele = h_tree->beginRootElement(),
    end_ele = h_tree->endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    relabelHGeometryRecursively(*the_ele);
  }

  the_ele = h_tree->beginRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    reclaimHGeometryRecursively(*the_ele);
  }
}

TEMPLATE
void THIS::reclaim()
{
  /// 先对各个IrregularMesh的不使用内存进行回收
  typename std::list<typename THIS::ir_mesh_t *>::iterator 
    the_ir_mesh = ir_mesh.begin(),
    end_ir_mesh = ir_mesh.end();
  for (;the_ir_mesh != end_ir_mesh;++ the_ir_mesh) {
    reclaimIrregularMesh(**the_ir_mesh);
  }
  
  /// 将所有组件标志为 -1
  initialTreeLabel();

  /// 对使用中的组件标识为 1
  the_ir_mesh = ir_mesh.begin();
  for (;the_ir_mesh != end_ir_mesh;++ the_ir_mesh) {
    labelIrregularMesh(**the_ir_mesh);
  }

  /// 对内存进行回收
  reclaimTreeMemory();
}

#undef THIS
#undef TEMPLATE

/**
 * end of file
 * 
 */
