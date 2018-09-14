/**
 * @file   MPI_Migration.h
 * @author Ruo Li <rli@aztec>
 * @date   Mon Oct 19 15:42:05 2009
 * 
 * @brief  
 * 
 * 
 */

#ifndef __MPI_Migration_h__
#define __MPI_Migration_h__

#include <mpi.h>
#include "../../include/Miscellaneous.h"
#include "../../include/Migration.details.h"
#include "../../include/HGeometry.h"
#include "../../include/mpi/MPI_HGeometry.h"
#include "../../include/mpi/MPI_SyncProp.h"
#include "../../include/FEMSpace.h"

AFEPACK_OPEN_NAMESPACE

namespace Migration {
  namespace details {
    //@{
    /**
     * 将几何体上的数据缓冲去除。
     */
    template <int DOW> void
      clear_geo_data_buffer(HGeometry<0,DOW>& geo) {
      geo.buffer.clear();
    }
    template <int DOW> void
      clear_geo_data_buffer(HGeometry<1,DOW>& geo) {
      geo.buffer.clear();
      for (u_int i = 0;i < geo.n_vertex;++ i) {
        clear_geo_data_buffer(*geo.vertex[i]);
      }
      if (geo.isRefined()) {
        for (u_int i = 0;i < geo.n_child;++ i) {
          clear_geo_data_buffer(*geo.child[i]);
        }
      }
    }
    template <class GEO> void
      clear_geo_data_buffer(GEO& geo) {
      geo.buffer.clear();
      for (u_int i = 0;i < geo.n_boundary;++ i) {
        clear_geo_data_buffer(*geo.boundary[i]);
      }
      if (geo.isRefined()) {
        for (u_int i = 0;i < geo.n_child;++ i) {
          clear_geo_data_buffer<GEO>(*geo.child[i]);
        }
      }
    }
    //@}

    //@{
    /**
     * 共享几何体上的数据缓冲进行同步。
     */
    struct buffer_trait {
      buffer_t * p_buf;

      friend Migration::ostream<>&
      operator<<(Migration::ostream<>& os,
                 const buffer_trait& bt) {
        const buffer_t& buf = *(bt.p_buf);
        u_int n_data = buf.size();
        os << n_data;
        buffer_t::const_iterator
          the_data = buf.begin(),
          end_data = buf.end();
        for (;the_data != end_data;++ the_data) {
          os << the_data->first << the_data->second;
        }
        return os;
      }
      friend Migration::istream<>&
      operator>>(Migration::istream<>& is,
                 buffer_trait& bt) {
        buffer_t& buf = *(bt.p_buf);
        data_buffer_t spool;
        u_int n_data;
        is >> n_data;
        for (u_int i = 0;i < n_data;++ i) {
          data_id_t id;
          is >> id;
          if (buf.find(id) == buf.end()) {
            is >> buf[id]; /// 如果自身没有该数据则接收
          } else { /// 否则丢弃掉接收的数据
            is >> spool;
          }
        }
        return is;
      }
    };

    template <int D, class HTREE> void
      sync_data_buffer_prepare(const HTREE& htree,
                               const property_id_t<buffer_trait>& pid) {
      enum {dow = HTREE::dow};
      typedef HGeometry<D,dow> GEO;
      const MPI::Shared_ptr_list<GEO>& so = htree.template get_shared_list<D>();
      typename MPI::Shared_ptr_list<GEO>::const_iterator
        the_entry = so.begin(), end_entry = so.end();
      for (;the_entry != end_entry;++ the_entry) {
        GEO * p_geo = the_entry->local_pointer();
        if (p_geo->get_property(pid) == NULL) {
          buffer_trait * p_bt = p_geo->new_property(pid);
          p_bt->p_buf = &(p_geo->buffer);
        }
      }
    }

    template <class HTREE> void
      sync_data_buffer_prepare(const HTREE& htree,
                               const property_id_t<buffer_trait>& pid,
                               int D) {
      switch (D) {
      case 0: sync_data_buffer_prepare<0,HTREE>(htree, pid); break;
      case 1: sync_data_buffer_prepare<1,HTREE>(htree, pid); break;
      case 2: sync_data_buffer_prepare<2,HTREE>(htree, pid); break;
      case 3: sync_data_buffer_prepare<3,HTREE>(htree, pid); break;
      }
    }

    template <class HTREE> void 
      sync_data_buffer(const HTREE& htree) {
      enum {dim = HTREE::dim};
      property_id_t<buffer_trait> pid;
      new_property_id(pid);

      MPI::PropSyncer<HTREE,buffer_trait> syncer(htree, pid);
      for (int i = 0;i <= dim;++ i) {
        sync_data_buffer_prepare(htree, pid, i);
        syncer.sync(i);
      }
      free_property_id(pid);
    }

    //@}


    //@{
    /**
     * 对几何体上的一个性质数据进行输出和输入。此功能可能将几何体上分配
     * 的性质 ID 进行输入和输出，但是其能力是非常有限的，如果性质数据本
     * 身含有指针对象结构的网络，其行为是不正确的。要求对于性质数据本身
     * 定义了对 AFEPack 内部流的输入输出算子。
     */
    template <class HGEO, class T>
      void export_property(const HGEO& geo,
                           const data_id_t& data_id,
                           const property_id_t<T>& pid,
                           const property_id_t<>& flag) {
      if (geo.has_property(flag)) return; /// 已经输出过了
      if (! geo.has_property(pid)) return; /// 并无要输出的性质

      geo.new_property(flag); /// 做上标记并进行输出

      T * p_prp = geo.get_property(pid);
      Migration::ostream<> os;
      os.set_buffer(get_buffer(geo, data_id, true));
      os << *p_prp;

      for (u_int i = 0;i < geo.n_vertex;++ i) {
        export_property(*geo.vertex[i], data_id, pid, flag);
      }

      for (u_int i = 0;i < geo.n_boundary;++ i) {
        export_property(*geo.boundary[i], data_id, pid, flag);
      }

      if (geo.isRefined()) {
        for (u_int i = 0;i < geo.n_child;++ i) {
          export_property(*geo.child[i], data_id, pid, flag);
        }
      }
    }

    template <class HGEO, class T>
      void import_property(const HGEO& geo,
                           const data_id_t& data_id,
                           const property_id_t<T>& pid,
                           const property_id_t<>& flag) {
      if (geo.has_property(flag)) return; /// 已经输入过了

      BinaryBuffer<>& buf = get_buffer(geo, data_id, false);
      if (&buf == NULL) return; /// 并无要输入的性质

      geo.new_property(flag); /// 做上标记并进行输出

      T * p_prp = geo.get_property(pid);
      if (p_prp == NULL) {
        p_prp = geo.new_property(pid);
      }

      Migration::istream<> is;
      is.set_buffer(buf);
      is >> *p_prp;

      for (u_int i = 0;i < geo.n_vertex;++ i) {
        import_property(*geo.vertex[i], data_id, pid, flag);
      }

      for (u_int i = 0;i < geo.n_boundary;++ i) {
        import_property(*geo.boundary[i], data_id, pid, flag);
      }

      if (geo.isRefined()) {
        for (u_int i = 0;i < geo.n_child;++ i) {
          import_property(*geo.child[i], data_id, pid, flag);
        }
      }
    }

    template <class HTREE, class T>
      void export_property(const HTREE& tree,
                           const data_id_t& data_id,
                           const property_id_t<T>& pid) {
      property_id_t<> flag;
      new_property_id(flag);

      enum {dim = HTREE::dim, dow = HTREE::dow};
      typename HTREE::ConstRootIterator 
        the_ele = tree.beginRootElement(),
        end_ele = tree.endRootElement();
      for (;the_ele != end_ele;++ the_ele) {
        const HGeometry<dim,dow> * p_geo = &(*the_ele);
        do {
          if (p_geo->parent != NULL) {
            p_geo = p_geo->parent;
          } else break;
        } while (true);
        export_property(*p_geo, data_id, pid, flag);
      }
    }

    template <class HTREE, class T>
      void import_property(const HTREE& tree,
                           const data_id_t& data_id,
                           const property_id_t<T>& pid) {
      property_id_t<> flag;
      new_property_id(flag);

      enum {dim = HTREE::dim, dow = HTREE::dow};
      typename HTREE::ConstRootIterator 
        the_ele = tree.beginRootElement(),
        end_ele = tree.endRootElement();
      for (;the_ele != end_ele;++ the_ele) {
        const HGeometry<dim,dow> * p_geo = &(*the_ele);
        do {
          if (p_geo->parent != NULL) {
            p_geo = p_geo->parent;
          } else break;
        } while (true);
        import_property(*p_geo, data_id, pid, flag);
      }
    }

    //@}
  }

  /**
   * 将一棵树上的所有数据缓冲清除，将内存释放出来。在数据迁移过后，应
   * 该进行此操作以释放内存；在新一轮数据迁移进行以前，如果数据缓冲没
   * 有清除，则必须进行此操作，否则将会出现数据混乱。
   */
  template <class HTREE> void
    clear_data_buffer(HTREE& tree) {
    enum {dim = HTREE::dim, dow = HTREE::dow};
    typename HTREE::RootIterator 
      the_ele = tree.beginRootElement(),
      end_ele = tree.endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      HGeometry<dim,dow> * p_geo = &(*the_ele);
      do {
        if (p_geo->parent != NULL) {
          p_geo = p_geo->parent;
        } else break;
      } while (true);
      details::clear_geo_data_buffer(*p_geo);
    }
  }

  void initialize(MPI_Comm comm); /// 对迁移环境进行初始化
  data_id_t register_data_name(const data_name_t& dn); /// 登记数据名称
  void load_config(const std::string& filename); /// 载入配置文件
  void save_config(const std::string& filename); /// 输出配置文件

  void ensured_open_fstream(const std::string& filename,
                            std::ifstream& is);
  void ensured_open_filtering_stream(const std::string& filename,
                                     filtering_istream& is);

  /** 
   * 输出一个有限元函数。对于附着于同一个几何体上的自由度，我们使用自由
   * 度的识别协议旗标和插值点来进行识别。
   * 
   * @param fun 有限元函数
   * @param data_id 数据 ID
   */  
  template <class FUNC>
    void export_fe_func(const FUNC& fun,
                        const data_id_t& data_id) {
    typedef FUNC fe_func_t;
    typedef typename fe_func_t::fe_space_t fe_space_t;
    typedef typename fe_space_t::element_t element_t;
    typedef typename fe_space_t::dof_info_t dof_info_t;
    typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
    typedef typename mesh_t::point_t point_t;
    const fe_space_t& sp = fun.femSpace();
    const mesh_t& mesh = dynamic_cast<const mesh_t&>(sp.mesh());
    Migration::ostream<> os;
    u_int n_dof = sp.n_dof();
    std::vector<double> h(n_dof, 0.0);
    typename fe_space_t::ConstElementIterator
      the_ele = sp.beginElement(),
      end_ele = sp.endElement();
    for (;the_ele != end_ele;++ the_ele) {
      const GeometryBM& geo = the_ele->geometry();
      const point_t& p0 = mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0));
      const point_t& p1 = mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0));
      double ele_h = distance(p0, p1);
      const std::vector<int>& ele_dof = the_ele->dof();
      int n_ele_dof = ele_dof.size();
      for (int i = 0;i < n_ele_dof;++ i) {
	double& dof_h = h[ele_dof[i]];
	dof_h = std::max(ele_h, dof_h);
      }
    }
    for (u_int i = 0;i < n_dof;++ i) {
      const dof_info_t& di = sp.dofInfo(i);
      const DOFIndex& dof_idx = sp.dofIndex(i);
      get_export_stream(mesh,
                        data_id,
                        dof_idx.dimension,
                        dof_idx.geometry_index,
                        os);
      os << di.identity
         << di.interp_point 
	 << h[i]
         << fun(i);
    }
  }
  /** 
   * 输出一个有限元函数的向量。对于附着于同一个几何体上的自由度，我们使用自由
   * 度的识别协议旗标和插值点来进行识别。
   *
   * 特别注意：此向量中的所有有限元函数需要建立在同一个有限元空间中！函
   * 数中并为对此进行检查，用户需自行保持这一点。
   *
   * 此函数由邓剑提供。
   * 
   * @param fun 有限元函数的向量
   * @param data_id 数据 ID
   */  
  template <class FUNC>
    void export_fe_func(const std::vector<FUNC>& fun,
                        const data_id_t& data_id) {
    typedef FUNC fe_func_t;
    typedef typename fe_func_t::fe_space_t fe_space_t;
    typedef typename fe_space_t::dof_info_t dof_info_t;
    typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
    typedef typename mesh_t::point_t point_t;
    const fe_space_t& sp = fun[0].femSpace();
    const mesh_t& mesh = dynamic_cast<const mesh_t&>(sp.mesh());
    Migration::ostream<> os;
    u_int n_dof = sp.n_dof();
    std::vector<double> h(n_dof, 0.0);
    typename fe_space_t::ConstElementIterator
      the_ele = sp.beginElement(),
      end_ele = sp.endElement();
    for (;the_ele != end_ele;++ the_ele) {
      const GeometryBM& geo = the_ele->geometry();
      const point_t& p0 = mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0));
      const point_t& p1 = mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0));
      double ele_h = distance(p0, p1);
      const std::vector<int>& ele_dof = the_ele->dof();
      int n_ele_dof = ele_dof.size();
      for (int i = 0;i < n_ele_dof;++ i) {
	double& dof_h = h[ele_dof[i]];
	dof_h = std::max(ele_h, dof_h);
      }
    }
    for (u_int i = 0;i < n_dof;++ i) {
      const dof_info_t& di = sp.dofInfo(i);
      const DOFIndex& dof_idx = sp.dofIndex(i);
      get_export_stream(mesh,
                        data_id,
                        dof_idx.dimension,
                        dof_idx.geometry_index,
                        os);
      os << di.identity
         << di.interp_point
	 << h[i];
      for (u_int j = 0;j < fun.size();++ j)
        os << fun[j](i);
    }
  }

  /** 
   * 载入一个有限元函数。对于附着于同一个几何体上的自由度，我们使用自由
   * 度的识别协议旗标和插值点的局部指标来进行识别。
   * 
   * @param fun 有限元函数
   * @param data_id 数据 ID
   */  
  template <class FUNC>
    void import_fe_func(FUNC& fun,
                        const data_id_t& data_id) {
    typedef FUNC fe_func_t;
    typedef typename fe_func_t::fe_space_t fe_space_t;
    typedef typename fe_space_t::dof_info_t dof_info_t;
    typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
    const fe_space_t& sp = fun.femSpace();
    const mesh_t& mesh = dynamic_cast<const mesh_t&>(sp.mesh());
    Migration::istream<> is;
    u_int n_dof = sp.n_dof();
    for (u_int i = 0;i < n_dof;++ i) {
      const dof_info_t& di = sp.dofInfo(i);
      const DOFIndex& dof_idx = sp.dofIndex(i);
      get_import_stream(mesh,
                        data_id,
                        dof_idx.dimension,
                        dof_idx.geometry_index,
                        is);
      do {
        dof_info_t dof_info;
	double h;
        is >> dof_info.identity
           >> dof_info.interp_point 
	   >> h
           >> fun(i);
        if (!(dof_info.identity == di.identity)) continue;
        if (distance(dof_info.interp_point, 
                     di.interp_point) <= 1.0e-04*h) break;
      } while (1);
    }
  }

  /** 
   * 载入一个有限元函数的向量。对于附着于同一个几何体上的自由度，我们使用自由
   * 度的识别协议旗标和插值点的局部指标来进行识别。
   *
   * 特别注意：此向量中的所有有限元函数需要建立在同一个有限元空间中！函
   * 数中并为对此进行检查，用户需自行保持这一点。
   *
   * 此函数由邓剑提供。
   * 
   * @param fun 有限元函数的向量
   * @param data_id 数据 ID
   */  
  template <class FUNC>
    void import_fe_func(std::vector<FUNC> & fun,
                        const data_id_t& data_id) {
    typedef FUNC fe_func_t;
    typedef typename fe_func_t::fe_space_t fe_space_t;
    typedef typename fe_space_t::dof_info_t dof_info_t;
    typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
    const fe_space_t& sp = fun[0].femSpace();
    const mesh_t& mesh = dynamic_cast<const mesh_t&>(sp.mesh());
    Migration::istream<> is;
    u_int n_dof = sp.n_dof();
    for (u_int i = 0;i < n_dof;++ i) {
      const dof_info_t& di = sp.dofInfo(i);
      const DOFIndex& dof_idx = sp.dofIndex(i);
      get_import_stream(mesh,
                        data_id,
                        dof_idx.dimension,
                        dof_idx.geometry_index,
                        is);
      do {
        dof_info_t dof_info;
	double h;
        is >> dof_info.identity
           >> dof_info.interp_point
	   >> h;
        for (u_int j = 0;j < fun.size();++ j)
          is >> fun[j](i);
        if (!(dof_info.identity == di.identity)) continue;
        if (distance(dof_info.interp_point, 
                     di.interp_point) <= 1.0e-04*h) break;
      } while (1);
    }
  }

  template <class HTREE, class T>
    void export_property(const HTREE& htree,
                         const data_id_t& data_id,
                         const property_id_t<T>& pid) {
    details::export_property(htree, data_id, pid);
  }

  template <class HTREE, class T>
    void import_property(const HTREE& htree,
                         const data_id_t& data_id,
                         const property_id_t<T>& pid) {
    details::import_property(htree, data_id, pid);
  }

  template <class HTREE>
    void sync_data_buffer(const HTREE& htree) {
    details::sync_data_buffer(htree);
  }

}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_Migration_h__

/**
 * end of file
 * 
 */
