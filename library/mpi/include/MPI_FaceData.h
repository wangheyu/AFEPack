/**
 * @file   MPI_FaceData.h
 * @author Ruo Li <rli@aztec>
 * @date   Thu Oct 22 20:48:38 2009
 * 
 * @brief  
 * 
 * 
 */

#ifndef __MPI_FaceData_h__
#define __MPI_FaceData_h__

#include "MPI_HGeometry.h"
#include "../../include/DGFEMSpace.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {
  namespace FaceData {

    /**
     * Syncer 主要用来提供一个方便的对面上的双边数据进行交换的接口，这
     * 是有限体积类型方法中最主要的数据交换方式。其特点在于，每个需要交
     * 换的数据是双边的，即只需要发送给一个邻居分区和从一个邻居分区获取。
     * 如果数据事实上不是双边的，操作结果是未知的。
     *
     * FOREST 模板参数用来指定所依赖的几何遗传树。PACKER 模板参数是完成
     * 交换的最主要的实现部分，我们要求其提供需要进行交换的数据类型
     *
     * typename PACKER::data_t;
     *
     * 和提供用户方面接口的
     *
     * DATA * PACKER::new_property(const OBJ&, const property_id_t<DATA>&) const 
     * DATA * PACKER::get_property(const OBJ&, const property_id_t<DATA>&) cosnt 
     *
     * 成员函数，其中 OBJ 是模板。PACKER 的设计主要由库完成，对于用户这
     * 应该是困难的工作。为了方便，我们提供另外一个包装方案Packer，使得
     * 对于 RegularMesh 和 DGFEMSpace 这两种实用情况的得以实现。
     *
     * 下面是对这个类进行使用的例子代码：
     *
     *   /// 需要交换的数据类型为 std::vector<double>
     *   typedef MPI::FaceData::FESpace<std::vector<double>,fe_space_t> data_packer_t;
     *   typedef MPI::FaceData::Syncer<forest_t,typename data_packer_t::packer_t> syncer_t;
     *
     *   syncer_t syncer(htree, data_packer_t::get_packer(fem_space));
     *   for (the_dgele 是单元之间的边界) {
     *     if (the_dgele 是一个需要交换数据的分区之间的边) {
     *       /// 取出需要交换的数据缓冲区
     *       std::vector<double> * p_buf = syncer.get_send_buffer(*the_dgele);
     *       ... ... /// 对需要交换的数据进行填充
     *     }
     *   }
     *   syncer.sync(); /// 交换数据
     *   for (the_dgele 是单元之间的边界) {
     *     if (the_dgele 是一个需要交换数据的分区之间的边) {
     *       /// 取出需要交换的数据缓冲区
     *       std::vector<double> * p_buf = syncer.get_recv_buffer(*the_dgele);
     *       ... ... /// 对交换后的数据进行解码
     *     }
     *   }
     *     
     */
    template <class FOREST, class PACKER>
      class Syncer {
    public:
      enum { dim = FOREST::dim, dow = FOREST::dow };
      typedef FOREST forest_t;
      typedef PACKER packer_t;
      typedef typename packer_t::data_t data_t;
      typedef HGeometry<dim-1,dow> geometry_t;

    private:
      typedef Syncer<forest_t,packer_t> this_t;

      const forest_t * _forest;
      packer_t _packer;

      property_id_t<data_t> _pid_out;
      property_id_t<data_t> _pid_in;

      Transmit_map<geometry_t> _transmit_map;

    public:
      Syncer() { new_property(); }
      Syncer(const forest_t& forest) : _forest(&forest) {
        new_property();
      }
      Syncer(const forest_t& forest, const packer_t& packer) :
      _forest(&forest), _packer(packer) {
        new_property();
      }
      Syncer(const this_t& syncer) :
      _forest(syncer._forest), _packer(syncer._packer) {
        new_property();
      }

    private:
      void new_property() {
        new_property_id(_pid_in);
        new_property_id(_pid_out);
      }
      void free_property() {
        free_property_id(_pid_in);
        free_property_id(_pid_out);
      }

    public:
      void set_forest(const forest_t& forest) { _forest = &forest; }
      const forest_t& forest() const { return *_forest; }

      void set_packer(const packer_t& packer) { 
        _packer = packer; 
        free_data_buffer();
      }
      const typename packer_t::packer_t * packer() const {
        return _packer.packer();
      }

      /**
       * 获取依附于对象上的数据发送缓冲区。OBJ 可以是面几何体或者
       * DGElement。如果缓冲区没有分配，则会进行分配。
       */
      template <class OBJ> data_t *
        get_send_buffer(const OBJ& obj) const {
        data_t * p_data = _packer.get_property(obj, _pid_out);
        if (p_data == NULL) p_data = _packer.new_property(obj, _pid_out);
        return p_data;
      }
      /**
       * 获取依附于对象上的数据接收缓冲区。OBJ 可以是面几何体或者
       * DGElement。
       */
      template <class OBJ> data_t *
        get_recv_buffer(const OBJ& obj) const {
        data_t * p_data = _packer.get_property(obj, _pid_in);
        return p_data;
      }

      /**
       * 显式释放数据缓冲区。
       */
      void free_data_buffer() {
        free_property();
        new_property();
      }

      /**
       * 更新数据交换图表。
       */
      void update_transmit_map() {
        _transmit_map.build(_forest->template get_shared_list<dim-1>(),
                            *this, &this_t::is_pack_data);
      }

      /**
       * 完成数据的交换。在数据交换前，应将发送数据填入数据发送缓冲区，
       * 在交换后可从数据接收缓冲区中获取解收到的数据。在本类的对象自
       * 身析构后，数据缓冲区将自行释放。数据缓冲区是可以重用的，但是
       * 也可以显式调用 free_data_buffer 将所有数据缓冲区释放。
       */
      void sync(bool is_update_transmit_map = true) {
        if (is_update_transmit_map) {
          update_transmit_map();
        }

        sync_data(_forest->communicator(), 
                  _transmit_map, 
                  *this, 
                  &this_t::pack_data, 
                  &this_t::unpack_data);
      }

      /**
       * 此处数据应该是由使用方填入，但并不是所有的面几何体均有此数据。
       * 此处判断仅仅对有数据发送的几何体进行操作。
       */
      bool is_pack_data(geometry_t * geo) const {
        assert (_forest->get_shared_info(*geo) != NULL);
        return (geo->get_property(_pid_out) != NULL);
      }

      void pack_data(geometry_t * geo,
                     int remote_rank,
                     Migration::ostream<>& os) {
        assert (_forest->get_shared_info(*geo) != NULL);
        data_t * p_data = geo->get_property(_pid_out);
        assert (p_data != NULL);
        os << *p_data;
      }
      void unpack_data(geometry_t * geo,
                       int remote_rank,
                       Migration::istream<>& is) {
        assert (_forest->get_shared_info(*geo) != NULL);
        data_t * p_data = geo->get_property(_pid_in);
        if (p_data == NULL) p_data = geo->new_property(_pid_in);
        is >> *p_data;
      }
    };

    /**
     * SyncerPtr 和 Syncer 的不同之处是实际交换数据的空间是用户自己分配
     * 的，SyncerPtr 自己仅仅记录了用户的数据空间的指针。在部分情况下，
     * 这将会减少一次数据拷贝，使得执行效率得以提升。
     * 
     * 请注意下面的例子代码和 Syncer 的不同：
     *
     *   /// 需要交换的数据类型为 std::vector<double>，此处需要在类型中加上指针
     *   typedef MPI::FaceData::FESpace<std::vector<double> *,fe_space_t> data_packer_t;
     *   typedef MPI::FaceData::SyncerPtr<forest_t,typename data_packer_t::packer_t> syncer_t;
     *
     *   syncer_t syncer(htree, data_packer_t::get_packer(fem_space));
     *   for (the_dgele 是单元之间的边界) {
     *     if (the_dgele 是一个需要交换数据的分区之间的边) {
     *       /// 用户需要自行准备发送数据的缓冲区
     *       std::vector<double> * p_send_buf = new std::vector<double)();
     *       ... ... /// 对需要发送的数据进行填充
     *       syncer.attach_send_buffer(*the_dgele, p_send_buf);
     *
     *       /// 与此同时，用户需要自行准备接收数据缓冲区
     *       std::vector<double> * p_recv_buf = new std::vector<double)();
     *       syncer.attach_send_buffer(*the_dgele, p_recv_buf);
     *     }
     *   }
     *   syncer.sync(); /// 交换数据
     *
     */
    template <class FOREST, class PACKER>
      class SyncerPtr {
    public:
      enum { dim = FOREST::dim, dow = FOREST::dow };
      typedef FOREST forest_t;
      typedef PACKER packer_t;
      typedef typename packer_t::data_t data_t; /// 是一个指针类型
      typedef HGeometry<dim-1,dow> geometry_t;

    private:
      typedef SyncerPtr<forest_t,packer_t> this_t;

      packer_t _packer;
      const forest_t * _forest;

      property_id_t<data_t> _pid_out;
      property_id_t<data_t> _pid_in;

      Transmit_map<geometry_t> _transmit_map;

    public:
      SyncerPtr() { new_property(); }
      SyncerPtr(const forest_t& forest) : _forest(&forest) {
        new_property();
      }
      SyncerPtr(const forest_t& forest, const packer_t& packer) :
      _forest(&forest), _packer(packer) {
        new_property();
      }
      SyncerPtr(const this_t& syncer) :
      _forest(syncer._forest), _packer(syncer._packer) {
        new_property();
      }

    private:
      void new_property() {
        new_property_id(_pid_in);
        new_property_id(_pid_out);
      }
      void free_property() {
        free_property_id(_pid_in);
        free_property_id(_pid_out);
      }

    public:
      void set_forest(const forest_t& forest) { _forest = &forest; }
      const forest_t& forest() const { return *_forest; }

      void set_packer(const packer_t& packer) { 
        free_data_buffer();
        _packer = packer; 
      }
      const typename packer_t::packer_t * packer() const {
        return _packer.packer();
      }

      /**
       * 获取依附于对象上的数据发送缓冲区。OBJ 可以是面几何体或者
       * DGElement。如果缓冲区没有分配，则会进行分配。
       */
      template <class OBJ> void
        attach_send_buffer(const OBJ& obj, const data_t data) const {
        data_t * p_data = _packer.get_property(obj, _pid_out);
        if (p_data == NULL) p_data = _packer.new_property(obj, _pid_out);
        (*p_data) = data;
      }
      /**
       * 获取依附于对象上的数据接收缓冲区。OBJ 可以是面几何体或者
       * DGElement。如果缓冲区没有分配，则会进行分配。
       */
      template <class OBJ> void
        attach_recv_buffer(const OBJ& obj, const data_t data) const {
        data_t * p_data = _packer.get_property(obj, _pid_in);
        if (p_data == NULL) p_data = _packer.new_property(obj, _pid_in);
        (*p_data) = data;
      }

      /**
       * 显式释放数据缓冲区。
       */
      void free_data_buffer() {
        free_property();
        new_property();
      }

      void update_transmit_map() {
        _transmit_map.build(_forest->template get_shared_list<dim-1>(),
                            *this, &this_t::is_pack_data);
      }

      /**
       * 完成数据的交换。在数据交换前，应将发送数据填入数据发送缓冲区，
       * 在交换后可从数据接收缓冲区中获取解收到的数据。在本类的对象自
       * 身析构后，数据缓冲区将自行释放。数据缓冲区是可以重用的，但是
       * 也可以显式调用 free_data_buffer 将所有数据缓冲区释放。
       */
      void sync(bool is_update_transmit_map = true) {
        if (is_update_transmit_map) {
          update_transmit_map();
        }

        sync_data(_forest->communicator(), 
                  _transmit_map, 
                  *this, 
                  &this_t::pack_data, 
                  &this_t::unpack_data);
      }

      /**
       * 此处数据应该是由使用方填入，但并不是所有的面几何体均有此数据。
       * 此处判断仅仅对有数据发送的几何体进行操作。
       */
      bool is_pack_data(geometry_t * geo) const {
        return (geo->get_property(_pid_out) != NULL);
      }

      void pack_data(geometry_t * geo,
                     int remote_rank,
                     Migration::ostream<>& os) {
        const data_t * p_data = geo->get_property(_pid_out);
        assert (p_data != NULL);
        assert (geo->get_property(_pid_in) != NULL);
        os << **p_data;
      }
      void unpack_data(geometry_t * geo,
                       int remote_rank,
                       Migration::istream<>& is) {
        const data_t * p_data = geo->get_property(_pid_in);
        assert (p_data != NULL);
        assert (geo->get_property(_pid_out) != NULL);
        is >> **p_data;
      }
    };

    namespace details {
      template <class OBJ, class DATA>
        struct _dummy_packer {
          DATA * _dummy_func(const OBJ&, const property_id_t<DATA>&) { return NULL; }
        };

      /**
       * 以成员函数做函数指针构造的双边数据打包器。
       */
      template <class OBJ, class DATA, 
        class PACKER = _dummy_packer<OBJ,DATA> >
        class Packer {
      public:
      typedef OBJ object_t;
      typedef DATA data_t;
      typedef PACKER packer_t;
      typedef property_id_t<data_t> pid_t;
      typedef data_t * (packer_t::*fun_ptr_t)(const object_t&, const pid_t&) const;

      private:
      typedef Packer<object_t,data_t,packer_t> this_t;
      const packer_t * _packer;
      fun_ptr_t _new_property;
      fun_ptr_t _get_property;

      public:
      Packer() {}
      Packer(const packer_t& pac, fun_ptr_t _np, fun_ptr_t _gp) :
      _packer(&pac), _new_property(_np), _get_property(_gp) {}
      Packer(const this_t& fdp) : _packer(fdp._packer),
      _new_property(fdp._new_property), _get_property(fdp._get_property) {}
      this_t& operator=(const this_t& fdp) {
        _packer = fdp._packer;
        _new_property = fdp._new_property;
        _get_property = fdp._get_property;
        return *this;
      }

      public:
      const packer_t * packer() const { return _packer; }
      data_t * get_property(const object_t& obj, const pid_t& pid) const {
        return (_packer->*_get_property)(obj, pid);
      };
      data_t * new_property(const object_t& obj, const pid_t& pid) const {
        return (_packer->*_new_property)(obj, pid);
      };
      };

      /**
       * 以非成员函数做函数指针构造的双边数据打包器。
       */
      template <class OBJ, class DATA>
        class Packer<OBJ,DATA,_dummy_packer<OBJ,DATA> > {
      public:
        typedef OBJ object_t;
        typedef DATA data_t;
        typedef _dummy_packer<object_t,DATA> packer_t;
        typedef property_id_t<data_t> pid_t;
        typedef data_t * (*fun_ptr_t)(const object_t&, const pid_t&);
      private:
        typedef Packer<object_t,data_t,packer_t> this_t;
        fun_ptr_t _new_property;
        fun_ptr_t _get_property;
      public:
        Packer() {}
      Packer(fun_ptr_t _np, fun_ptr_t _gp) :
        _new_property(_np), _get_property(_gp) {}
      Packer(const this_t& fdp) :
        _new_property(fdp._new_property), _get_property(fdp._get_property) {}
        this_t& operator=(const this_t& fdp) {
          _new_property = fdp._new_property;
          _get_property = fdp._get_property;
          return *this;
        }
      public:
        const packer_t * packer() const { return NULL; }
        data_t * get_property(const object_t& obj, const pid_t& pid) const {
          return (*_get_property)(obj, pid);
        };
        data_t * new_property(const object_t& obj, const pid_t& pid) const {
          return (*_new_property)(obj, pid);
        };
      };

    } // namespace details

    /**
     * 用于网格的双边数据打包器。提供类型 packer_t 和打包器构造函数
     * get_packer。网格必须是 RegularMesh。
     */
    template <class DATA, class MESH>
      struct Mesh {
        typedef MESH mesh_t;
        typedef typename mesh_t::ir_mesh_t ir_mesh_t;
        typedef details::Packer<GeometryBM, DATA, mesh_t> packer_t;
        typedef HGeometry<mesh_t::dim-1,mesh_t::dow> h_geometry_t;
        static packer_t get_packer(const mesh_t& mesh) {
          return packer_t(mesh, 
                          &mesh_t::template new_property<DATA,mesh_t::dim-1>,
                          &mesh_t::template get_property<DATA,mesh_t::dim-1>);
        }
        static packer_t get_packer(const ir_mesh_t& mesh) {
          return packer_t(mesh.regularMesh(), 
                          &mesh_t::template new_property<DATA,mesh_t::dim-1>,
                          &mesh_t::template get_property<DATA,mesh_t::dim-1>);
        }
      };

    /**
     * 用于有限元空间的双边数据打包器。提供类型 packer_t 和打包器构造函数
     * get_packer。有限元空间必须是 DGFEMSpace。
     */
    template <class DATA, class SPACE>
      struct FESpace {
        typedef SPACE fe_space_t;
        typedef typename fe_space_t::dg_element_t dg_element_t;
        typedef details::Packer<dg_element_t, DATA, fe_space_t> packer_t;
        typedef HGeometry<fe_space_t::tdim1, fe_space_t::dow> h_geometry_t;
        static packer_t get_packer(const fe_space_t& sp) {
          return packer_t(sp,
                          &fe_space_t::template new_property<DATA>,
                          &fe_space_t::template get_property<DATA>);
        }
      };

  } // namespace FaceData
} // namespace MPI

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_FaceData_h__

/**
 * end of file
 * 
 */
