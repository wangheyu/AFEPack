/**
 * @file   MPI_MemoryReclaimer.h
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Wed Jan 12 09:47:48 2011
 * 
 * @brief  此处是专门为使用 ULoadBalance 做了负载平衡以后产生的网格做
 *         内存回收的时候使用的。为了保险起见，我们对所有的共享几何体
 *         均未进行回收。
 * 
 */

#ifndef __MPI_MemoryReclaimer_h__
#define __MPI_MemoryReclaimer_h__

#include "MPI_HGeometry.h"

/*!

  此处所使用的算法和串行的 MemoryReclaimer 是一样的，所不同之处是对于共
  享几何体将不会进行回收。共享几何体在使用 ULoadBalance 做了负载平衡以
  后，情形变得很复杂，我们很难实现共享几何体的拷贝都同时回收，因此对所
  有的共享几何体都不进行回收。为了不回收共享几何体，但是又保持树结构的
  要求

      "一个几何体在某个分区上，则其组件和兄弟都必然在此分区上有拷贝"

  我们先通过下列步骤将所有不回收的几何体上都做一个标记，在后面的步骤
  中，我们检查此标记，以确保这些几何体不会被回收。步骤为：

  1. 将所有自身或者其组件为共享几何体的单元几何体进行标记，并将这个标
     记加到其兄弟身上；

  2. 将所有标记过的单元几何体的所有组件都做上标记；

  对所有加上了标记的几何体，通过在relabelHGeometryRecursively 函数中，
  将几何体的 index 改为 -2 的时候，对几何体是否是共享体或者是几何体的构
  建是共享几何体做一次判断，如果是，就不将 index 修改为-2。如果 index
  不是-2，几何体就不会被回收。共享的或者其构建是共享的几何体的回收将在
  负载平衡的操作之后，共享几何体变为不是这样的几何体以后才进行。因此此
  处的内存回收不是非常完全，但是由于共享几何体相对数量较少，总体上应该
  还是可以接受的。

*/

AFEPACK_OPEN_NAMESPACE

namespace MPI {
  template <class FOREST>
    class MemoryReclaimer {
  public:
    enum { dim = FOREST::dim, dow = FOREST::dow };
    typedef FOREST tree_t;
    typedef BirdView<tree_t> ir_mesh_t;
  private:
    tree_t * h_tree;
    std::list<ir_mesh_t *> ir_mesh;

  public:
  MemoryReclaimer() : h_tree(NULL) {}
  MemoryReclaimer(tree_t& _h_tree) : h_tree(&_h_tree) {}
    virtual ~MemoryReclaimer() {}

  public:
    /**
     * 设置进行内存回收操作的HGeometryTree。
     * 
     */
    void setGeometryTree(tree_t& _h_tree) {
      h_tree = &_h_tree;
    }
    /**
     * 加进一个参加操作的IrregularMesh。
     * 
     */
    void addIrregularMesh(ir_mesh_t& _ir_mesh) {
      if (&(_ir_mesh.getForest()) != h_tree) {
        std::cout << "warning: the irregular mesh added is not based on the geometry tree used."
                  << std::endl;
      }
      ir_mesh.push_back(&_ir_mesh);
    }
    void clear() {
      h_tree = NULL;
      ir_mesh.clear();
    }
    void reclaim() {
      /// 先对各个IrregularMesh的不使用内存进行回收
      typename std::list<ir_mesh_t *>::iterator 
        the_ir_mesh = ir_mesh.begin(),
        end_ir_mesh = ir_mesh.end();
      for (;the_ir_mesh != end_ir_mesh;++ the_ir_mesh) {
        reclaimIrregularMesh(**the_ir_mesh);
      }
  
      /// 对和共享几何体有关的几何体做上标记
      new_property_id(_pid_is_shared);
      markSharedObject();

      /// 将所有组件标志为 -1
      initialTreeLabel();

      /// 对使用中的组件标识为 1
      the_ir_mesh = ir_mesh.begin();
      for (;the_ir_mesh != end_ir_mesh;++ the_ir_mesh) {
        labelIrregularMesh(**the_ir_mesh);
      }

      /// 对内存进行回收
      reclaimTreeMemory();
      free_property_id(_pid_is_shared);
    }

  private:
    void markSharedObject() {
      typename HGeometryTree<dim,dow>::RootIterator 
        the_ele = h_tree->beginRootElement(),
        end_ele = h_tree->endRootElement();
      for (;the_ele != end_ele;++ the_ele) {
        HGeometry<dim,dow> * p_geo = &(*the_ele);
        do {
          if (p_geo->parent != NULL) {
            p_geo = p_geo->parent;
          } else break;
        } while (true);
        markSharedElement(*p_geo);
      }

      the_ele = h_tree->beginRootElement();
      for (;the_ele != end_ele;++ the_ele) {
        HGeometry<dim,dow> * p_geo = &(*the_ele);
        do {
          if (p_geo->parent != NULL) {
            p_geo = p_geo->parent;
          } else break;
        } while (true);
        markSharedGeometry(*p_geo);
      }
    }      

    void reclaimIrregularMesh(ir_mesh_t& m) {
      ActiveElementIterator<dim,dow>
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

    void initialTreeLabel() {
      typename HGeometryTree<dim,dow>::RootIterator 
        the_ele = h_tree->beginRootElement(),
        end_ele = h_tree->endRootElement();
      for (;the_ele != end_ele;++ the_ele) {
        labelHGeometryRecursively(*the_ele, -1);
      }
    }

    void labelIrregularMesh(ir_mesh_t& m) {
      RootFirstElementIterator<dim,dow>
        the_ele = m.beginRootFirstElement(),
        end_ele = m.endRootFirstElement();
      for (;the_ele != end_ele;++ the_ele) {
        labelHGeometry(*(the_ele->h_element), 1);
      }

      typename ir_mesh_t::mesh_t& mesh = m.regularMesh();
#define LABEL_USED_HGEOMETRY(D)                                         \
      if (dim >= D) {                                                   \
        int n_geo = mesh.h_geometry()[D].size();                        \
        for (int i = 0;i < n_geo;++ i) {                                \
          HGeometry<D,dow> * p_geo = mesh.template h_geometry<D>(i);    \
          labelHGeometry(*p_geo, 1);                                    \
        }                                                               \
      }

      LABEL_USED_HGEOMETRY(0);
      LABEL_USED_HGEOMETRY(1);
      LABEL_USED_HGEOMETRY(2);
      LABEL_USED_HGEOMETRY(3);

#undef LABEL_USED_HGEOMETRY
    }

    void reclaimTreeMemory() {
      typename HGeometryTree<dim,dow>::RootIterator
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

  private:
    property_id_t<> _pid_is_shared;

    /**
     * 判断一个几何体是否是共享的，或者具有共享组件。
     */
    template <class GEO> bool
      is_shared_geometry(const GEO& g) const {
      bool result = false;
      result |= (h_tree->get_shared_info(g) != NULL);
      if (! result) {
        for (int i = 0;i < g.n_vertex;++ i) {
          result |= is_shared_geometry(*(g.vertex[i]));
        }
      }
      if (result) {
        return result;
      } else {
        for (int i = 0;i < g.n_boundary;++ i) {
          result |= is_shared_geometry(*(g.boundary[i]));
        }
      }
      return result;
    }

    /**
     * 对单元 ele 递归搜索设置标记。如果自身需要标记，则将其兄弟也标记
     * 上，然后对其后代递归进行。
     */
    void markSharedElement(HGeometry<dim,dow>& ele) const {
      if ((ele.get_property(_pid_is_shared) == NULL) &&
          (is_shared_geometry(ele))) {
        ele.new_property(_pid_is_shared);

        if (ele.parent != NULL) {
          for (int i = 0;i < ele.parent->n_child;++ i) {
            HGeometry<dim,dow> * p_sibling = ele.parent->child[i];
            if (p_sibling == &ele) continue;
            if (p_sibling->get_property(_pid_is_shared) == NULL) {
              p_sibling->new_property(_pid_is_shared);
            }
          }
        }
      }
      
      if (ele.isRefined()) {
        for (int i = 0;i < ele.n_child;++ i) {
          markSharedElement(*ele.child[i]);
        }
      }
    }

    /**
     * 如果几何体已经被标记，则将其所有顶点和边界标记，并对边界递归进
     * 行标记。然后对自身的后代递归进行。
     */
    template <class GEO> void
      markSharedGeometry(GEO& g) const {
      if (g.get_property(_pid_is_shared) != NULL) {
        for (int i = 0;i < g.n_vertex;++ i) {
          if (g.vertex[i]->get_property(_pid_is_shared) == NULL) {
            g.vertex[i]->new_property(_pid_is_shared);
          }
        }
        for (int i = 0;i < g.n_boundary;++ i) {
          if (g.boundary[i]->get_property(_pid_is_shared) == NULL) {
            g.boundary[i]->new_property(_pid_is_shared);
          }
          markSharedGeometry(*g.boundary[i]);
        }
      }
      if (g.isRefined()) {
        for (int i = 0;i < g.n_child;++ i) {
          markSharedGeometry(*g.child[i]);
        }
      }
    }      

    template <class GEO> void 
      labelHGeometry(GEO& g, int lab) const {
      for (int i = 0;i < g.n_vertex;++ i) {
        labelHGeometry(*(g.vertex[i]), lab);
      }
      for (int i = 0;i < g.n_boundary;++ i) {
        labelHGeometry(*(g.boundary[i]), lab);
      }
      g.index = lab;
    }

    template <class GEO> void 
      labelHGeometryRecursively(GEO& g, int lab) const {
      labelHGeometry(g, lab);

      if (g.isRefined()) {
        for (int i = 0;i < g.n_child;++ i) {
          labelHGeometryRecursively(*(g.child[i]), lab);
        }
      }
    }

    template <class GEO> int 
      relabelHGeometryRecursively(GEO& g) const {
      if (g.get_property(_pid_is_shared) != NULL) {
        g.index = 1;
      }

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
      if (g.index == -1) {
        g.index = -2;
        return -1;
      } else {
        return g.index;
      }
    }

    template <class GEO> int 
      reclaimHGeometryRecursively(GEO& g) const {
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

      if ((g.index == -2) && 
          (g.get_property(_pid_is_shared) == NULL)) {
        this->reclaimHGeometry(&g);
        return -1;
      } else {
        return 1;
      }
    }

    template <class GEO> void 
      reclaimHGeometry(GEO * p_geo) const {
      assert ((p_geo->get_property(_pid_is_shared) == NULL));
      delete p_geo;
    }
  };

}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_MemoryReclaimer_h__

/**
 * end of file
 * 
 */
