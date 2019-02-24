/**
 * @file   Serialization.h
 * @author Ruo Li <rli@aztec>
 * @date   Wed Sep 23 15:25:54 2009
 * 
 * @brief 在区域分解并行实现中，为了完成动态负载平衡中的数据迁移，我们
 *        将必要的数据结构进行序列化，然后在不同的进程间进行数据传输，
 *        然后将数据恢复。
 *
 * 为了完成这个部分的工作，我们使用 boost 的序列化库来实现。对于下面的
 * 一些主要数据结构，我们手工实现了序列化的程序。这些类包括：
 *
 *   - HGeometryTree : 几何遗传树
 *   - IrregularMesh : 非正则网格
 *   - FEMFunction : 有限元函数
 * 
 */

#ifdef __SERIALIZATION__
#ifndef __Serialization_h__
#define __Serialization_h__

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

#include "HGeometry.h"

namespace boost {
  namespace serialization {

    /**
     * 序列化类Point<DOW>。
     */
    template <class AR, int DOW> void
      serialize(AR& ar, 
                AFEPack::Point<DOW>& pnt, 
                const u_int version) {
      for (int i = 0;i < DOW;++ i) {
        ar & pnt[i];
      }
    }

    /**
     * 序列化类BinaryBuffer<T>。
     * 
     */
    template <class AR, class T> void
      serialize(AR& ar, 
                BinaryBuffer<T>& buf, 
                const u_int version) {
      ar & boost::serialization::base_object<typename BinaryBuffer<T>::_Base>(buf);;
    }

    /**
     * 序列化类HGeometry<DIM,DOW>。
     */
    template <class AR, int DIM, int DOW> void 
      serialize(AR& ar, 
                HGeometry<DIM,DOW>& hg, 
                const u_int version) {
      ar & hg.buffer;

      ar & hg.index;
      ar & hg.vertex;
      ar & hg.boundary;
      ar & hg.parent;
      ar & hg.child;
      ar & hg.bmark;
    }

    /**
     * 系列化类HGeometry<0,DOW>。
     */
    template <class AR, int DOW> void
      serialize(AR& ar, 
                HGeometry<0,DOW>& hg, 
                const u_int version) {
      ar & hg.buffer;

      ar & boost::serialization::base_object<AFEPack::Point<DOW> >(hg);
      ar & hg.index;
      ar & hg.bmark;
    }

    /**
     * 序列化类HGeometryTree<DIM,DOW>。
     * 
     * @param ar 档案
     * @param ht 几何遗传树
     * @param version 版本
     */
    template <class AR, int DIM, int DOW> void
      serialize(AR& ar, 
                HGeometryTree<DIM,DOW>& ht, 
                const u_int version) {
      ar & ht.rootElement();
      ar & ht.is_locked();
    }

    /**
     * 从档案中载入HGeometryTree<DIM,DOW>的部分单元。
     * 
     * @param ar 输入型档案
     * @param tree 几何遗传树
     * @param version 版本
     */
    template <class AR, int DIM, int DOW> void
      load(AR& ar,
           HGeometryTree<DIM,DOW>& tree,
           const u_int version = 0) {
      u_int n_root_ele;
      ar & n_root_ele;
      for (u_int i = 0;i < n_root_ele;++ i) {
        HGeometry<DIM,DOW> * ele = new HGeometry<DIM,DOW>();
        ar & *ele;
        tree.rootElement().push_back(ele);
      }
    }


    /**
     * 向档案中写入HGeometryTree<DIM,DOW>的部分单元。
     * 
     * @param ar 输出型档案
     * @param tree 几何遗传树
     * @param start 开始序号
     * @param end 结束序号
     * @param version 版本
     */
    template <class AR, int DIM, int DOW> void
      save(AR& ar,
           HGeometryTree<DIM,DOW>& tree,
           int start, 
           int end,
           const u_int version = 0) {
      u_int n_root_ele = end - start;
      typename HGeometryTree<DIM,DOW>::RootIterator
        the_root_ele = tree.beginRootElement();
      std::advance(the_root_ele, start);

      ar & n_root_ele;
      for (u_int i = 0;i < n_root_ele;++ i) {
        ar & *the_root_ele;
        ++ the_root_ele;
      }
    }

    /**
     * 序列化类HElement<DIM,DOW>。
     */
    template <class AR, int DIM, int DOW> void 
      serialize(AR& ar, 
                HElement<DIM,DOW>& he, 
                const u_int version) {
      ar & he.index;
      ar & he.indicator;
      ar & he.value;
      ar & he.h_element;
      ar & he.parent;
      ar & he.child;
    }

    /**
     * 序列化类IrregularMesh<DIM,DOW>。
     * 
     * @param ar 档案
     * @param im 非正则网格
     * @param version 版本
     */
    template <class AR, int DIM, int DOW> void
      serialize(AR& ar, 
                IrregularMesh<DIM,DOW>& im, 
                const u_int version) {
      ar & im.rootElement();
      ar & &(im.geometryTree());
    }

    /**
     * 从档案中载入IrregularMesh<DIM,DOW>的部分单元。
     * 
     * @param ar 输入型档案
     * @param im 非正则网格
     * @param version 版本
     */
    template <class AR, int DIM, int DOW> void
      load(AR& ar,
           IrregularMesh<DIM,DOW>& im,
           const u_int version = 0) {
      u_int n_root_ele;
      ar & n_root_ele;
      for (u_int i = 0;i < n_root_ele;++ i) {
        HElement<DIM,DOW> * ele = new HElement<DIM,DOW>();
        ar & *ele;
        im.rootElement().push_back(ele);
      }
    }


    /**
     * 向档案中写入IrregularMesh<DIM,DOW>的部分单元。
     * 
     * @param ar 输出型档案
     * @param im 非正则网格
     * @param start 开始序号
     * @param end 结束序号
     * @param version 版本
     */
    template <class AR, int DIM, int DOW> void
      save(AR& ar,
           IrregularMesh<DIM,DOW>& im,
           int start, 
           int end,
           const u_int version = 0) {
      u_int n_root_ele = end - start;
      typename IrregularMesh<DIM,DOW>::RootIterator
        the_root_ele = im.beginRootElement();
      std::advance(the_root_ele, start);

      ar & n_root_ele;
      for (u_int i = 0;i < n_root_ele;++ i) {
        ar & *the_root_ele;
        ++ the_root_ele;
      }
    }

  }
}
#endif // __Serialization_h__
#endif // __SERIALIZATION__
/**
 * end of file
 * 
 */
