/**
 * @file   MPI_SyncDOF.h
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Mon Jan  3 12:29:14 2011
 * 
 * @brief  对函数的自由度进行同步。
 * 
 * 
 */

#ifndef __MPI_SyncDOF_h__
#define __MPI_SyncDOF_h__

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  template <class FOREST, class FESPACE>
    class DOFSyncer {
  private:
    typedef FOREST forest_t;
    typedef FESPACE fe_space_t;

    const forest_t * _p_forest;
    const fe_space_t * _p_fe_sp;

  public:
    DOFSyncer(const forest_t& forest,
              const fe_space_t& fe_sp) :
    _p_forest(&forest), _p_fe_sp(&fe_sp) {}

  private:
    template <class FEFUN>
    bool encode_geometry_data(const FEFUN& fun,
                              int i,
                              const property_id_t<BinaryBuffer<> >& pid) const {
      const DOFIndex& di = _p_fe_sp->dofIndex(i);
      const int& geo_dim = di.dimension;
      bool result = false;
      switch (geo_dim) {
      case 0: result = encode_geometry_data<0,FEFUN>(fun, i, di, pid); break;
      case 1: result = encode_geometry_data<1,FEFUN>(fun, i, di, pid); break;
      case 2: result = encode_geometry_data<2,FEFUN>(fun, i, di, pid); break;
      case 3: result = encode_geometry_data<3,FEFUN>(fun, i, di, pid); break;
      }
      return result;
    }

    template <int D, class FEFUN>
      bool encode_geometry_data(const FEFUN& fun,
                                int i,
                                const DOFIndex& di,
                                const property_id_t<BinaryBuffer<> >& pid) const {
      typedef RegularMesh<forest_t::dim,forest_t::dow> mesh_t;
      typedef typename fe_space_t::dof_info_t dof_info_t;
      const int& geo_idx = di.geometry_index;
      const mesh_t& mesh = dynamic_cast<const mesh_t&>(_p_fe_sp->mesh());
      const HGeometry<D,forest_t::dow> * p_geo = mesh.template h_geometry<D>(geo_idx);
      if (_p_forest->get_shared_info(*p_geo) == NULL) return false;

      BinaryBuffer<> * p_buf = p_geo->get_property(pid);
      if (p_buf == NULL) {
        p_buf = p_geo->new_property(pid);
      }

      Migration::ostream<> os(*p_buf);
      const dof_info_t& dof_info = _p_fe_sp->dofInfo(i);
      os << dof_info.identity << dof_info.interp_point << fun(i);
      return true;
    }

    template <class FEFUN>
    void decode_geometry_data(FEFUN& fun,
                              int i,
                              const property_id_t<BinaryBuffer<> >& pid) const {
      const DOFIndex& di = _p_fe_sp->dofIndex(i);
      const int& geo_dim = di.dimension;
      switch (geo_dim) {
      case 0: decode_geometry_data<0,FEFUN>(fun, i, di, pid); break;
      case 1: decode_geometry_data<1,FEFUN>(fun, i, di, pid); break;
      case 2: decode_geometry_data<2,FEFUN>(fun, i, di, pid); break;
      case 3: decode_geometry_data<3,FEFUN>(fun, i, di, pid); break;
      }
    }

    template <int D, class FEFUN>
      void decode_geometry_data(FEFUN& fun,
                                int i,
                                const DOFIndex& di,
                                const property_id_t<BinaryBuffer<> >& pid) const {
      typedef RegularMesh<forest_t::dim,forest_t::dow> mesh_t;
      typedef typename fe_space_t::dof_info_t dof_info_t;
      const int& geo_idx = di.geometry_index;
      const mesh_t& mesh = dynamic_cast<const mesh_t&>(_p_fe_sp->mesh());
      const HGeometry<D,forest_t::dow> * p_geo = mesh.template h_geometry<D>(geo_idx);
      if (_p_forest->get_shared_info(*p_geo) == NULL) return;

      BinaryBuffer<> * p_buf = p_geo->get_property(pid);
      if (p_buf == NULL) return;

      const dof_info_t& dof_info = _p_fe_sp->dofInfo(i);
      double d = dof_info.interp_point.length();
      Migration::istream<> is(*p_buf);
      do {
        dof_info_t dof_info_1;
        is >> dof_info_1.identity >> dof_info_1.interp_point >> fun(i);
        if (!(dof_info.identity == dof_info_1.identity)) continue;
        double d1 = dof_info_1.interp_point.length() + d;
        if (distance(dof_info.interp_point,
                     dof_info_1.interp_point) <= 1.0e-08*d1) break;
      } while (true);
    }

  public:
    /**
     * 对有限元函数的自由度进行分区间的同步，其中 mask 是对需要同步的自
     * 由度的遮罩：在发送的时候，如果 mask[i] 为真，则发送第 i 个自由
     * 度；在接收的时候，如果 mask[i] 为假，并且对面发送了第 i 个自由
     * 度，才会接收第 i 个自由度。
     */
    template <class FEFUN, class MASK>
      void sync(FEFUN& fun, const MASK& mask) const {
      property_id_t<BinaryBuffer<> > pid;
      new_property_id(pid);

      std::vector<int> dim_mask(forest_t::dim + 1, 0);
      u_int n_dof = _p_fe_sp->n_dof();
      for (u_int i = 0;i < n_dof;++ i) {
        if (mask[i]) {
          if (encode_geometry_data(fun, i, pid)) {
            const DOFIndex& di = _p_fe_sp->dofIndex(i);
            const int& geo_dim = di.dimension;
            dim_mask[geo_dim] = 1;
          }
        }
      }

      PropSyncer<forest_t,BinaryBuffer<> > syncer(*_p_forest, pid);
      for (u_int i = 0;i <= forest_t::dim;++ i) {
        int is_sync;
        MPI_Allreduce(&dim_mask[i], &is_sync, 1, MPI_INT, 
                      MPI_SUM, _p_forest->communicator());
        if (is_sync > 0) syncer.sync(i);
      }

      for (u_int i = 0;i < n_dof;++ i) {
        if (! mask[i]) {
          decode_geometry_data(fun, i, pid);
        }
      }
      free_property_id(pid);
    }
  };
}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_SyncDOF_h__

/**
 * end of file
 * 
 */
