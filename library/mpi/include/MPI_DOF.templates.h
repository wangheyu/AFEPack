/**
 * @file   MPI_DOF.templates.h
 * @author Ruo Li <rli@aztec>
 * @date   Mon Oct 26 17:22:11 2009
 * 
 * @brief  
 * 
 * 
 */

#define TEMPLATE template <class FOREST, class FESPACE>
#define THIS GlobalIndex<FOREST,FESPACE>

TEMPLATE int
THIS::get_dof_minimal_rank(u_int i) const {
  typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
  const mesh_t& mesh = dynamic_cast<const mesh_t&>(_fe_sp->mesh());

  const DOFIndex& dof_idx = _fe_sp->dofIndex(i);
  const int& dim = dof_idx.dimension;
  const int& gid = dof_idx.geometry_index;

  int * p_min_rank = mesh.template get_property<int>(dim, gid, 
						     _pid_min_rank);
  if (p_min_rank == NULL) {
    return _forest->rank();
  } else {
    return *p_min_rank;
  }
}

TEMPLATE bool
THIS::is_dof_on_minimal_geometry(u_int i) const {
  return _idopg[i];
}

/// 对所有共享自由度的编号在全局进行同步
TEMPLATE void
THIS::sync_idx() {
  u_int n_dof = _fe_sp->n_dof();
  typedef typename fe_space_t::value_t value_t;
  typedef typename fe_space_t::dof_info_t dof_info_t;
  typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
  const mesh_t& mesh = dynamic_cast<const mesh_t&>(_fe_sp->mesh());

  /// 首先准备发送数据
  new_property_id(_pid_global_idx);
  BinaryBuffer<> * p_buf;
  for (u_int i = 0;i < n_dof;++ i) {
    if ((*this)(i) == UNUSED_IDX) continue;

    const dof_info_t& di = _fe_sp->dofInfo(i);
    const DOFIndex& dof_idx = _fe_sp->dofIndex(i);
    const int& dim = dof_idx.dimension;
    const int& gid = dof_idx.geometry_index;
    if (_forest->is_geometry_shared(mesh, dim, gid)) {
      p_buf = mesh.get_property(dim, gid, _pid_global_idx);
      if (p_buf == NULL) {
        p_buf = mesh.new_property(dim, gid, _pid_global_idx);
      }
      Migration::ostream<> os(*p_buf);
      os << di.interp_point << di.identity << (*this)(i);
    }
  }

  /// 对数据进行交换
  MPI_Comm comm = _forest->communicator();
#define MPI_DOF_SYNC_DATA_DIM(D)                                        \
  if (forest_t::dim >= D) {                                             \
    sync_data(comm, _forest->template get_shared_list<D>(), *this,      \
              &this_t::template pack_global_idx<D>,                     \
              &this_t::template unpack_global_idx<D>);                  \
  }

  MPI_DOF_SYNC_DATA_DIM(0);
  MPI_DOF_SYNC_DATA_DIM(1);
  MPI_DOF_SYNC_DATA_DIM(2);
  MPI_DOF_SYNC_DATA_DIM(3);

#undef MPI_DOF_SYNC_DATA_DIM

  /// 对接收到的数据进行解码
  const typename forest_t::matcher_t& matcher = _forest->matcher();

  std::vector<bool> flag(n_dof, false);
  typename fe_space_t::ConstElementIterator
    the_ele = _fe_sp->beginElement(),
    end_ele = _fe_sp->endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    u_int n_ele_dof = ele_dof.size();
    for (u_int j = 0;j < n_ele_dof;++ j) {
      const int& i = ele_dof[j];
      if (flag[i] || (*this)(i) != UNUSED_IDX) continue;

      const dof_info_t& di = _fe_sp->dofInfo(i);
      const DOFIndex& dof_idx = _fe_sp->dofIndex(i);
      const int& dim = dof_idx.dimension;
      const int& gid = dof_idx.geometry_index;
      if (_forest->is_geometry_shared(mesh, dim, gid)) {
        flag[i] = true;

        const GeometryBM& geo = the_ele->geometry();
        double h = distance(mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0)),
                            mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0)));

        p_buf = mesh.get_property(dim, gid, _pid_global_idx);
        assert (p_buf != NULL);
        Migration::istream<> is(*p_buf);
        do {
          dof_info_t dof_info;
          is >> dof_info.interp_point >> dof_info.identity >> (*this)(i);
          if (!(dof_info.identity == di.identity)) continue;
          if (matcher.value(dof_info.interp_point, di.interp_point, 1.0e-04*h) >= 0) break;
        } while (1);
      }
    }
  }

  free_property_id(_pid_global_idx);
}

namespace __details {
  struct min_rank {
    int operator()(int a, int b) const {
      return std::min(a, b);
    }
  };
};

TEMPLATE void
THIS::set_dof_minimal_rank() {
  free_property_id(_pid_min_rank);
  new_property_id(_pid_min_rank);

  typedef typename fe_space_t::value_t value_t;
  typedef typename fe_space_t::dof_info_t dof_info_t;
  typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
  const mesh_t& mesh = dynamic_cast<const mesh_t&>(_fe_sp->mesh());
  
  int rank_has_dim_dof[4] = {0, 0, 0, 0};
  u_int n_dof = _fe_sp->n_dof();
  for (u_int i = 0;i < n_dof;++ i) {
    const DOFIndex& dof_idx = _fe_sp->dofIndex(i);
    const int& dim = dof_idx.dimension;
    const int& gid = dof_idx.geometry_index;
    if (_forest->is_geometry_shared(mesh, dim, gid)) {
      rank_has_dim_dof[dim] = 1;

      int * p_min_rank = mesh.template get_property<int>(dim, gid,
							 _pid_min_rank);
      if (p_min_rank == NULL) {
	p_min_rank = mesh.template new_property<int>(dim, gid,
						     _pid_min_rank);
	(*p_min_rank) = _forest->rank();
      }
    }
  }

  /// 交换最小秩
  int has_dim_dof[4] = {0, 0, 0, 0};
  MPI_Allreduce(rank_has_dim_dof, has_dim_dof, 4, MPI_INT, MPI_MAX, _forest->communicator());
  PropSyncer<forest_t,int> syncer(*_forest, _pid_min_rank);
  for (u_int dim = 0;dim <= 3;++ dim) {
    if (has_dim_dof[dim] == 1) {
      syncer.template sync(dim, __details::min_rank());
    }
  }

  /// 清点在本进程上进行编号的自由度个数
  base_t::resize(n_dof, UNUSED_IDX); /// 为数组分配空间
  _idopg.resize(n_dof, false);
  int idx = 0;
  for (u_int i = 0;i < n_dof;++ i) {
    int min_rank = get_dof_minimal_rank(i);
    if (min_rank == _forest->rank()) {
      idx += 1;
      (*this)(i) = 0;
      _idopg[i] = true;
    }
  }
  _n_primary_dof = idx;
}

TEMPLATE void
THIS::build() {
  u_int n_dof = _fe_sp->n_dof();
  typedef typename fe_space_t::value_t value_t;
  typedef typename fe_space_t::dof_info_t dof_info_t;
  typedef RegularMesh<fe_space_t::dim,fe_space_t::dow> mesh_t;
  const mesh_t& mesh = dynamic_cast<const mesh_t&>(_fe_sp->mesh());

  /// 设置所有自由度的最小秩
  set_dof_minimal_rank();

  /// 获取本进程上全局编号的首个指标
  int idx = _n_primary_dof;
  int global_idx;
  MPI_Comm comm = _forest->communicator();
  u_int n_rank = _forest->n_rank();
  MPI_Scan(&idx, &global_idx, 1, MPI_INT, MPI_SUM, comm);
  _n_global_dof = global_idx; /// 获取总的自由度个数
  MPI_Bcast(&_n_global_dof, 1, MPI_UNSIGNED, n_rank - 1, comm);
  if (_forest->rank() == 0) {
    std::cerr << "Total #DOF: " << _n_global_dof << std::endl;
  }

  _first_idx = global_idx - idx; /// 本进程上的首指标
  /// 对需要在本进程上进行编号的自由度进行编号
  std::vector<int>::iterator
    the_idx = base_t::begin(),
    end_idx = base_t::end();
  for (idx = _first_idx;the_idx != end_idx;++ the_idx) {
    if (*the_idx == UNUSED_IDX) continue;
    *the_idx = idx ++;
  }

  sync_idx(); /// 进程间同步
}

TEMPLATE void
THIS::build_primary_index(int * idx_ptr) const {
  std::vector<int>::const_iterator
    the_idx = base_t::begin(),
    end_idx = base_t::end();
  for (u_int i = 0;the_idx != end_idx;++ the_idx) {
    if (is_dof_on_minimal_geometry(i ++)) {
      *idx_ptr = *the_idx;
      ++ idx_ptr;
    }
  }
}

TEMPLATE 
template <class MAP> void
THIS::build_epetra_map(MAP& map) const {
  int * idx_ptr = map.MyGlobalElements();
  std::vector<int>::const_iterator
    the_idx = base_t::begin(),
    end_idx = base_t::end();
  for (u_int i = 0;the_idx != end_idx;++ the_idx) {
    if (is_dof_on_minimal_geometry(i ++)) {
      *idx_ptr = *the_idx;
      ++ idx_ptr;
    }
  }
}

TEMPLATE 
template <class MAP, class INVMAP> void
THIS::build_epetra_map(MAP& map, INVMAP& inv_map) const {
  int * idx_ptr = map.MyGlobalElements();
  std::vector<int>::const_iterator
    the_idx = base_t::begin(),
    end_idx = base_t::end();
  for (u_int i = 0, j = 0;the_idx != end_idx;++ the_idx) {
    if (is_dof_on_minimal_geometry(i)) {
      *idx_ptr = *the_idx;
      inv_map[j ++] = i ++;
      ++ idx_ptr;
    }
  }
}

TEMPLATE
template <class LC, class GC> void
  THIS::translate(const LC& lc, GC& gc) const {
  u_int n = lc.size();
  gc.resize(n);
  for (u_int i = 0;i < n;++ i) {
    gc[i] = (*this)(lc[i]);
  }
}

TEMPLATE
template <int GDIM> void
THIS::pack_global_idx(HGeometry<GDIM,fe_space_t::dow> * geo,
                      int remote_rank,
                      Migration::ostream<>& os) {
  BinaryBuffer<> * p_buf = geo->get_property(_pid_global_idx);
  u_int n = 0;
  if (p_buf != NULL) {
    n = p_buf->size();
    os << n;
    os.encode_binary(&(*p_buf)[0], n);
  } else {
    os << n;
  }
}

TEMPLATE
template <int GDIM> void
THIS::unpack_global_idx(HGeometry<GDIM,fe_space_t::dow> * geo,
                        int remote_rank,
                        Migration::istream<>& is) {
  u_int n = 0;
  is >> n;
  if (n > 0) {
    BinaryBuffer<> * p_buf;;
    assert ((p_buf  = geo->get_property(_pid_global_idx)) == NULL);
    p_buf = geo->new_property(_pid_global_idx);
    p_buf->resize(n);
    is.decode_binary(&(*p_buf)[0], n);
  }
}

TEMPLATE
template <class INVEC, class OUTVEC>
  void THIS::scatter_hypre_vector(INVEC& iv,
                                  OUTVEC& ov) const {
  /**
   * 此函数的实现分为下面三步：
   *
   * 1. 分析本类所需要的远程元素的值为哪些元素和应该来自于哪些进程；
   *
   * 2. 向远程进程发送请求，要求对方向自身发送相应的值；
   *
   * 3. 获取远程进程发来的数据包，填充到本进程向量的相应位置上；
   * 
   */

  /// 分析需要接收的远程元素的位置和来源进程
  const int& n_rank = _forest->n_rank();
  int n_ldof = n_local_dof();
  std::vector<std::list<std::pair<int,int> > > rdof(n_rank);
  for (int i = 0;i < n_ldof;++ i) {
    if (is_dof_on_minimal_geometry(i)) continue;
    int remote_rank = get_dof_minimal_rank(i);
    rdof[remote_rank].push_back(std::pair<int,int>(i, (*this)(i)));
  }

  /// 准备发送给远程进程数据请求
  std::vector<int> target(n_rank);
  std::vector<BinaryBuffer<> > obuf(n_rank), ibuf(n_rank);
  for (int i = 0;i < n_rank;++ i) {
    target[i] = i;
    Migration::ostream<> os(obuf[i]);

    u_int n_item = rdof[i].size();
    os << n_item;
    typename std::list<std::pair<int,int> >::iterator
      the_dof = rdof[i].begin(),
      end_dof = rdof[i].end();
    for (;the_dof != end_dof;++ the_dof) {
      os << the_dof->first << the_dof->second;
    }
  }

  /// 发送和接收远程数据请求
  MPI_Barrier(_forest->communicator());
  sendrecv_data(_forest->communicator(), n_rank, obuf.begin(),
                ibuf.begin(), target.begin());

  /// 分析接收到的数据，并准备好发送的数据
  obuf.clear();
  obuf.resize(n_rank);
  for (int i = 0;i < n_rank;++ i) {
    Migration::istream<> is(ibuf[i]);
    Migration::ostream<> os(obuf[i]);

    u_int n_item;
    is >> n_item;
    os << n_item;
    if (n_item > 0) {
      std::vector<int> ldof(n_item);
      std::vector<int> gdof(n_item);
      std::vector<double> value(n_item);
      for (u_int j = 0;j < n_item;++ j) {
        is >> ldof[j] >> gdof[j];
      }
      HYPRE_IJVectorGetValues(iv, n_item, &gdof[0], &value[0]);
      for (u_int j = 0;j < n_item;++ j) {
        os << ldof[j] << value[j];
      }
    }
  }

  /// 发送和接收远程发来的向量数据
  ibuf.clear();
  ibuf.resize(n_rank);
  MPI_Barrier(_forest->communicator());
  sendrecv_data(_forest->communicator(), n_rank, obuf.begin(),
                ibuf.begin(), target.begin());

  /// 对接收的向量数据进行分析和填充
  for (int i = 0;i < n_rank;++ i) {
    Migration::istream<> is(ibuf[i]);

    u_int n_item;
    is >> n_item;

    int ldof;
    double value;

    for (u_int j = 0;j < n_item;++ j) {
      is >> ldof >> value;
      ov(ldof) = value;
    }
  }

  /**
   * 填充本来就在本地的数据，首先建立自身身上的全局指标到局部指标的映
   * 射表。
   */
  for (u_int i = 0, j = 0;i < n_ldof;++ i) {
    if (is_dof_on_minimal_geometry(i)) {
      HYPRE_IJVectorGetValues(iv, 1, &(*this)(i), &ov(i));
    }
  }
}

#undef THIS
#undef TEMPLATE 

/**
 * end of file
 * 
 */
