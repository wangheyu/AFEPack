/**
 * @file   MemoryReclaimer.h
 * @author Robert Lie
 * @date   Sun Apr 29 15:48:52 2007
 * 
 * @brief  
 * 
 * 
 */

#ifndef _MemoryReclaimer_h_
#define _MemoryReclaimer_h_

#include "HGeometry.h"

/*!

  在HGeometryTree和IrregularMesh中，当一个HGeometry和HElement曾经被细
  分了以后，我们将其后代存储在内存中，在这些对象暂时不使用的时候，我们
  并没有立刻将这些内存释放掉，其主要的考虑是为了在这些对象再次需要的时
  候，可以节省处理的时间。但是如果在特定的时候，如果期望将这些暂时不用
  的内存收回，可以本类来实现，典型的代码行如下：

  <pre>

    HGeometryTree<DIM,DOW> h_tree;
    ... ...
    IrregularMesh<DIM,DOW> ir_mesh_0;
    IrregularMesh<DIM,DOW> ir_mesh_1;
    ... ...
    MemoryReclaimer<DIM,DOW> mr(h_tree);
    mr.addIrregularMesh(ir_mesh_0);
    mr.addIrregularMesh(ir_mesh_1);
    mr.reclaim();

  </pre>

  在调用reclaim函数之后，不使用的内存将会被回收。

  需要特别注意的是：在做此操作时，所有建立在h_tree上的、在随后还会使用
  的IrregularMesh，都需要使用addIrregularMesh 函数加入到这个回收的过程
  中，否则可能出现潜在的错误。

 */
template <int DIM, int DOW=DIM>
  class MemoryReclaimer {
 public:
 enum {dim = DIM, dow = DOW};

 typedef HGeometryTree<DIM,DOW> tree_t;
 typedef IrregularMesh<DIM,DOW> ir_mesh_t;
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
   if (&(_ir_mesh.geometryTree()) != h_tree) {
     std::cout << "warning: the irregular mesh added is not based on the geometry tree used."
     << std::endl;
   }
   ir_mesh.push_back(&_ir_mesh);
 }
 void clear() {
   h_tree = NULL;
   ir_mesh.clear();
 }
 /*!

   HGeometryTree中间的所有几何体相互索引，成为一个非常复杂的网络结构。
   为了回收内存，我们使用下面的步骤来实现：

   1. 对每个IrregularMesh中不用的内存进行回收。这个操作是非常简单的，
   因为我们只需要把每个叶子节点的子树通通删除就够了，实现在函数
   reclaimIrregularMesh这个函数中。

   2. 将整个HGeometryTree中的所有几何体使用index=-1做标识，这个操作在
   initialTreeLabel中进行。

   3. 将所有的IrregularMesh中索引过的单元用到的几何体都重新标识为1，
   这个操作在 labelIrregularMesh 中做到。

   4. 在整个HGeometryTree中的几何体进行遍历：如果一个几何体的标识为-1，
   则该几何体应该被删除。在这个几何体第一次被索引到的时候，我们就将
   其标识修改为-2，并返回-1表示是第一次索引到它。如果标识为1，那么
   这个几何体正在被某个IrregularMesh使用，我们不做什么。如果一个几
   何体标识为-2，那么这个几何体已经是至少第二次被索引到了，那么就切
   断对其进行这次索引的联系，将对其进行索引的那个指针设为 NULL。经
   过这个操作，HGeometryTree中原本是个网状结构的数据，那个需要被删
   除的部分，已经被拆成为了树状结构；

   5. 在整个HGeometryTree中的几何体进行遍历，如果一个几何体的标识为-2，
   我们就将其删除。4和5两步都是在 reclaimTreeMemory 中实现的；

 */
 void reclaim();

 private:
 void reclaimIrregularMesh(ir_mesh_t&);
 void initialTreeLabel();
 void labelIrregularMesh(ir_mesh_t&);
 void reclaimTreeMemory();

 template <int DIM1> void labelHGeometry(HGeometry<DIM1,DOW>&, int lab);
 template <int DIM1> void labelHGeometryRecursively(HGeometry<DIM1,DOW>& g, int lab);
 template <int DIM1> int relabelHGeometryRecursively(HGeometry<DIM1,DOW>& g);
 template <int DIM1> int reclaimHGeometryRecursively(HGeometry<DIM1,DOW>& g);

 void labelHGeometry(HGeometry<0,DOW>&, int lab);
 void labelHGeometryRecursively(HGeometry<0,DOW>& g, int lab);
 int relabelHGeometryRecursively(HGeometry<0,DOW>& g);
 int reclaimHGeometryRecursively(HGeometry<0,DOW>& g);

 virtual void 
 reclaimHGeometry(void * p_geo, int dim) const {
   switch (dim) {
   case 0: delete ((HGeometry<0,DOW> *)(p_geo)); break;
   case 1: delete ((HGeometry<1,DOW> *)(p_geo)); break;
   case 2: delete ((HGeometry<2,DOW> *)(p_geo)); break;
   case 3: delete ((HGeometry<3,DOW> *)(p_geo)); break;
   }
 }

 protected:
 tree_t * get_tree_ptr() const { return h_tree; }
};

#include "MemoryReclaimer.templates.h"

#endif // _MemoryReclaimer_h_

/**
 * end of file
 * 
 */
