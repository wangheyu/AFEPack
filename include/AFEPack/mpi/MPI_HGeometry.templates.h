/**
 * @file   MPI_HGeometry.templates.h
 * @author Ruo Li <rli@aztec>
 * @date   Fri Sep 11 16:58:34 2009
 * 
 * @brief  
 * 
 * 
 */

#define TEMPLATE template <int DIM, int DOW, class MATCHER>
#define THIS HGeometryForest<DIM,DOW,MATCHER>

TEMPLATE void 
THIS::readMesh(const std::string& filename) {

  /**
   * 取秩和进程个数。
   */
  const int& rank = this->_rank;
  const int& n_rank = this->_n_rank;

  /**
   * 数据文件的名字缺省定为
   *
   *     filename + rank + ".mesh"
   *
   * 的格式。这个是从现在已经有的程序的习惯中继承过来的。现在的网格文
   * 件名称和原本AFEPack的网格文件的格式是完全相同的。
   */
  char filename_buf[256];
  if (n_rank > 1) {
    sprintf(filename_buf, "%s%d.mesh", filename.c_str(), rank);
  } else {
    sprintf(filename_buf, "%s.mesh", filename.c_str());
  }
  std::cerr << "Reading in mesh data file " 
            << filename_buf 
            << " as geometry tree root ..." 
            << std::endl;
  std::ifstream is(filename_buf);

  /// 读入顶点的坐标
  u_int n_point;
  is >> n_point;
  std::cerr << "\t# points: " << n_point << std::endl;
  std::vector<Point<this_t::dow> > point(n_point);
  for (u_int i = 0;i < n_point;i ++) is >> point[i];
  is >> n_point;
  std::vector<HGeometry<0,this_t::dow> *> geo_0d(n_point, (HGeometry<0,this_t::dow> *)NULL);
  for (u_int i = 0;i < n_point;i ++) {
    u_int j, k;
    is >> j; geo_0d[j] = new HGeometry<0,this_t::dow>();
    is >> k >> k; *dynamic_cast<Point<this_t::dow> *>(geo_0d[j]) = point[k];
    is >> k >> k >> geo_0d[j]->bmark;
  }
  point.clear();
	
#define GDIM 1
  u_int n_geo_1d;
  std::vector<HGeometry<GDIM,this_t::dow>*> geo_1d;
  if (DIM >= 1) {/// 读入一维几何体的信息
    is >> n_geo_1d;
    std::cerr << "\t# 1D-geometry: " << n_geo_1d << std::endl;
    geo_1d.resize(n_geo_1d, NULL);
    for (u_int i = 0;i < n_geo_1d;i ++) {
      u_int j, k, l;
      is >> j >> k; geo_1d[j] = new HGeometry<GDIM,this_t::dow>();
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
  std::vector<HGeometry<GDIM,this_t::dow>*> geo_2d;
  if (DIM >= 2) {/// 读入二维几何体的信息
    is >> n_geo_2d;
    std::cerr << "\t# 2D-geometry: " << n_geo_2d << std::endl;
    geo_2d.resize(n_geo_2d, NULL);
    for (u_int i = 0;i < n_geo_2d;i ++) {
      u_int j, k, l;
      is >> j >> k; geo_2d[j] = new HGeometry<GDIM,this_t::dow>();
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
  std::vector<HGeometry<GDIM,this_t::dow>*> geo_3d;
  if (DIM >= 3) { /// 读入三维几何体的信息
    is >> n_geo_3d;
    std::cerr << "\t# 3D-geometry: " << n_geo_3d << std::endl;
    geo_3d.resize(n_geo_3d, NULL);
    for (u_int i = 0;i < n_geo_3d;i ++) {
      u_int j, k, l;
      is >> j >> k; geo_3d[j] = new HGeometry<GDIM,this_t::dow>();
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
      this->rootElement().push_back((HGeometry<DIM,this_t::dow> *)geo_1d[i]);
    }
  } else if (DIM == 2) {
    for (u_int i = 0;i < n_geo_2d;++ i) {
      this->rootElement().push_back((HGeometry<DIM,this_t::dow> *)geo_2d[i]);
    }
  } else if (DIM == 3) {
    for (u_int i = 0;i < n_geo_3d;++ i) {
      this->rootElement().push_back((HGeometry<DIM,this_t::dow> *)geo_3d[i]);
    }
  }

  /**
   * 下面开始读入共享信息并构造最基本的背景网格的共享结构。
   *
   * 共享信息文件名为 filename + rank + ".share"
   */
  if (n_rank == 1) return; /// 只有一个进程时不需要 .share 文件

  sprintf(filename_buf, "%s%d.share", filename.c_str(), rank);
  is.open(filename_buf); /// 打开文件
  std::cerr << "Reading in shared data file " 
            << filename_buf
            << " as shared information ..." 
            << std::endl;

  /// 先准备好读入的数据的缓冲区和流
  std::vector<bool> proc_flag(n_rank, false);
  std::vector<BinaryBuffer<> > proc_buf_in(n_rank), proc_buf_out(n_rank);
  std::vector<Migration::ostream<> > proc_os(n_rank);
  for (u_int i = 0;i < n_rank;++ i) {
    proc_os[i].set_buffer(proc_buf_out[i]);
  }

  /// 读入节点信息
  int n_item;
  std::vector<std::pair<std::vector<std::pair<int,int> >, int> > 
    shared_geo_0d, shared_geo_1d, shared_geo_2d, shared_geo_3d;        
  if (DIM >= 0) {
    std::vector<int> n_rank_shared_0d(n_rank, 0);
    int n_shared_geo_0d;
    is >> n_shared_geo_0d >> n_shared_geo_0d;
    std::cerr << "\t# 0d-shared geometry: " << n_shared_geo_0d << std::endl;
    shared_geo_0d.resize(n_shared_geo_0d);
    for (int i = 0;i < n_shared_geo_0d;++ i) {
      is >> n_item; /// 共享份数
      shared_geo_0d[i].first.resize(n_item);
      for (int j = 0;j < n_item;++ j) {
        is >> shared_geo_0d[i].first[j].first   /// 几何体所在进程的秩
           >> shared_geo_0d[i].first[j].second; /// 几何体的指标
        n_rank_shared_0d[shared_geo_0d[i].first[j].first] += 1;
      }
      is >> shared_geo_0d[i].second; /// 边界标识
    }

    for (int i = 0;i < n_rank;++ i) {
      proc_os[i] << n_rank_shared_0d[i];
      if (n_rank_shared_0d[i] > 0) {
        proc_flag[i] = true;
      }
    }
    for (int i = 0;i < n_shared_geo_0d;++ i) {
      n_item = shared_geo_0d[i].first.size();
      int local_idx;
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_0d[i].first[j].first;
        int idx = shared_geo_0d[i].first[j].second;
        if (rnk == rank) {
          local_idx = idx;
          break;
        }
      }
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_0d[i].first[j].first;
        int idx = shared_geo_0d[i].first[j].second;
        if (rnk != rank) {
          proc_os[rnk] << idx << geo_0d[local_idx];
        }
      }
    }
  }

  if (DIM >= 1) {/// 一维共享几何体
    std::vector<int> n_rank_shared_1d(n_rank, 0);
    int n_shared_geo_1d;
    is >> n_shared_geo_1d >> n_shared_geo_1d;
    std::cerr << "\t# 1d-shared geometry: " << n_shared_geo_1d << std::endl;
    shared_geo_1d.resize(n_shared_geo_1d);        
    for (int i = 0;i < n_shared_geo_1d;++ i) {
      is >> n_item; /// 共享份数
      shared_geo_1d[i].first.resize(n_item);
      for (int j = 0;j < n_item;++ j) {
        is >> shared_geo_1d[i].first[j].first   /// 几何体所在进程的秩
           >> shared_geo_1d[i].first[j].second; /// 几何体的指标
        n_rank_shared_1d[shared_geo_1d[i].first[j].first] += 1;
      }
      is >> shared_geo_1d[i].second; /// 边界标识
    }

    for (int i = 0;i < n_rank;++ i) {
      proc_os[i] << n_rank_shared_1d[i];
      if (n_rank_shared_1d[i] > 0) {
        proc_flag[i] = true;
      }
    }
    for (int i = 0;i < n_shared_geo_1d;++ i) {
      n_item = shared_geo_1d[i].first.size();
      int local_idx;
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_1d[i].first[j].first;
        int idx = shared_geo_1d[i].first[j].second;
        if (rnk == rank) {
          local_idx = idx;
          break;
        }
      }
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_1d[i].first[j].first;
        int idx = shared_geo_1d[i].first[j].second;
        if (rnk != rank) {
          proc_os[rnk] << idx << geo_1d[local_idx];
        }
      }
    }
  }

  if (DIM >= 2) {/// 二维共享几何体
    std::vector<int> n_rank_shared_2d(n_rank, 0);
    int n_shared_geo_2d;
    is >> n_shared_geo_2d >> n_shared_geo_2d;
    std::cerr << "\t# 2d-shared geometry: " << n_shared_geo_2d << std::endl;
    shared_geo_2d.resize(n_shared_geo_2d);        
    for (int i = 0;i < n_shared_geo_2d;++ i) {
      is >> n_item; /// 共享份数
      shared_geo_2d[i].first.resize(n_item);
      for (int j = 0;j < n_item;++ j) {
        is >> shared_geo_2d[i].first[j].first   /// 几何体所在进程的秩
           >> shared_geo_2d[i].first[j].second; /// 几何体的指标
        n_rank_shared_2d[shared_geo_2d[i].first[j].first] += 1;
      }
      is >> shared_geo_2d[i].second; /// 边界标识
    }

    for (int i = 0;i < n_rank;++ i) {
      proc_os[i] << n_rank_shared_2d[i];
      if (n_rank_shared_2d[i] > 0) {
        proc_flag[i] = true;
      }
    }
    for (int i = 0;i < n_shared_geo_2d;++ i) {
      n_item = shared_geo_2d[i].first.size();
      int local_idx;
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_2d[i].first[j].first;
        int idx = shared_geo_2d[i].first[j].second;
        if (rnk == rank) {
          local_idx = idx;
          break;
        }
      }
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_2d[i].first[j].first;
        int idx = shared_geo_2d[i].first[j].second;
        if (rnk != rank) {
          proc_os[rnk] << idx << geo_2d[local_idx];
        }
      }
    }
  }

  if (DIM >= 3) {/// 三维共享几何体
    std::vector<int> n_rank_shared_3d(n_rank, 0);
    int n_shared_geo_3d;
    is >> n_shared_geo_3d >> n_shared_geo_3d;
    std::cerr << "\t# 3d-shared geometry: " << n_shared_geo_3d << std::endl;
    shared_geo_3d.resize(n_shared_geo_3d);        
    for (int i = 0;i < n_shared_geo_3d;++ i) {
      is >> n_item; /// 共享份数
      shared_geo_3d[i].first.resize(n_item);
      for (int j = 0;j < n_item;++ j) {
        is >> shared_geo_3d[i].first[j].first   /// 几何体所在进程的秩
           >> shared_geo_3d[i].first[j].second; /// 几何体的指标
        n_rank_shared_3d[shared_geo_3d[i].first[j].first] += 1;
      }
      is >> shared_geo_3d[i].second; /// 边界标识
    }

    for (int i = 0;i < n_rank;++ i) {
      proc_os[i] << n_rank_shared_3d[i];
      if (n_rank_shared_3d[i] > 0) {
        proc_flag[i] = true;
      }
    }
    for (int i = 0;i < n_shared_geo_3d;++ i) {
      n_item = shared_geo_3d[i].first.size();
      int local_idx;
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_3d[i].first[j].first;
        int idx = shared_geo_3d[i].first[j].second;
        if (rnk == rank) {
          local_idx = idx;
          break;
        }
      }
      for (int j = 0;j < n_item;++ j) {
        int rnk = shared_geo_3d[i].first[j].first;
        int idx = shared_geo_3d[i].first[j].second;
        if (rnk != rank) {
          proc_os[rnk] << idx << geo_3d[local_idx];
        }
      }
    }
  }

  /**
   * 下面发送和接收信息。
   */
  std::cerr << "\ttransfering data ..." << std::endl;
  int tag = 97, n_req = 0;
  proc_flag[rank] = false; /// 往本进程自身不用发送信息
  MPI_Request request[2*n_rank];
  MPI_Status status[2*n_rank];
  for (int i = 0;i < n_rank;++ i) {
    if (proc_flag[i] == false) continue;
    MPI_Isend(proc_buf_out[i].start_address(), proc_buf_out[i].size(), MPI_CHAR,
              i, tag, _comm, &request[n_req ++]);
  }
  for (int i = 0;i < n_rank;++ i) {
    MPI_Status mpi_status;
    if (proc_flag[i] == false) continue;
    MPI_Probe(i, tag, _comm, &mpi_status);
    int msg_size;
    MPI_Get_count(&mpi_status, MPI_CHAR, &msg_size);
    proc_buf_in[i].resize(msg_size);
    MPI_Irecv(proc_buf_in[i].start_address(), msg_size, MPI_CHAR,
              i, tag, _comm, &request[n_req ++]);
  }
  MPI_Waitall(n_req, request, status);

  /**
   * 下面对收到的数据进行解码。我们已知数据的发送都是成对进行的，即从
   * 进程 i 向进程 j 有数据发送，则从进程 j 向进程 i 也必定有数据发送。
   */
  std::cerr << "\tdecoding data received ..." << std::endl;
  for (int i = 0;i < n_rank;++ i) { /// 对进程进行循环
    if (proc_flag[i] == false) continue;
    Migration::istream<> proc_is(proc_buf_in[i]);
    int n_shared_item, local_idx;

    if (DIM >= 0) { /// 解码零维信息
#define GDIM 0
#define GEO HGeometry<GDIM,this_t::dow>
#define OBJ Shared_object<GEO>
      proc_is >> n_shared_item;
      for (int j = 0;j < n_shared_item;++ j) {
        GEO * remote_obj, * local_obj;
        proc_is >> local_idx >> remote_obj;
        local_obj = geo_0d[local_idx];
        OBJ * p_info = get_shared_info(*local_obj);
        if (p_info == NULL) p_info = new_shared_info(*local_obj);
        p_info->add_clone(i, remote_obj);
      }
#undef OBJ
#undef GEO
#undef GDIM
    }

    if (DIM >= 1) {/// 解码一维信息
#define GDIM 1
#define GEO HGeometry<GDIM,this_t::dow>
#define OBJ Shared_object<GEO>
      proc_is >> n_shared_item;
      for (int j = 0;j < n_shared_item;++ j) {
        GEO * remote_obj, * local_obj;
        proc_is >> local_idx >> remote_obj;
        local_obj = geo_1d[local_idx];
        OBJ * p_info = get_shared_info(*local_obj);
        if (p_info == NULL) p_info = new_shared_info(*local_obj);
        p_info->add_clone(i, remote_obj);
      }
#undef OBJ
#undef GEO
#undef GDIM
    }

    if (DIM >= 2) {/// 解码二维信息
#define GDIM 2
#define GEO HGeometry<GDIM,this_t::dow>
#define OBJ Shared_object<GEO>
      proc_is >> n_shared_item;
      for (int j = 0;j < n_shared_item;++ j) {
        GEO * remote_obj, * local_obj;
        proc_is >> local_idx >> remote_obj;
        local_obj = geo_2d[local_idx];
        OBJ * p_info = get_shared_info(*local_obj);
        if (p_info == NULL) p_info = new_shared_info(*local_obj);
        p_info->add_clone(i, remote_obj);
      }
#undef OBJ
#undef GEO
#undef GDIM
    }

    if (DIM >= 3) {/// 解码三维信息
#define GDIM 3
#define GEO HGeometry<GDIM,this_t::dow>
#define OBJ Shared_object<GEO>
      proc_is >> n_shared_item;
      for (int j = 0;j < n_shared_item;++ j) {
        GEO * remote_obj, * local_obj;
        proc_is >> local_idx >> remote_obj;
        local_obj = geo_3d[local_idx];
        OBJ * p_info = get_shared_info(*local_obj);
        if (p_info == NULL) p_info = new_shared_info(*local_obj);
        p_info->add_clone(i, remote_obj);
      }
#undef OBJ
#undef GEO
#undef GDIM
    }
  }

  set_shared_info_sent();
}

TEMPLATE
void THIS::renumerateRootElement(void (*f)(const double *, double *))
{
  std::cerr << "Renumerating element of the mesh using Hibert space filling curve ..." 
            << std::flush;
  int n_ele = this->n_rootElement();
  std::vector<HGeometry<DIM,this_t::dow>*> tmp_ele(n_ele);
  std::vector<std::vector<double> > x(DOW, 
                                      std::vector<double>(n_ele, 0.0));
  typename std::list<HGeometry<DIM,this_t::dow>*>::iterator
    the_ele = this->rootElement().begin();
  for (int i = 0;i < n_ele;++ i, ++ the_ele) {
    tmp_ele[i] = *the_ele;
    HGeometry<DIM,this_t::dow>& ele = **the_ele;
    const int& n_vtx = ele.n_vertex;
    Point<DOW> pnt;
    for (int j = 0;j < n_vtx;++ j) {
      pnt += *ele.vertex[j];
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

  the_ele = this->rootElement().begin();
  for (int i = 0;i < n_ele;++ i, ++ the_ele) {
    *the_ele = tmp_ele[new_index[i]];
  }
  std::cerr << " OK!" << std::endl;
}


TEMPLATE
void THIS::eraseRootElement(u_int level) {
  for (u_int i = 0;i < level;++ i) {
    this->eraseRootElementOneLevel();
  }
}

TEMPLATE
void THIS::eraseRootElementOneLevel() {
  /// 先收集背景单元的孩子
  std::list<HGeometry<DIM,this_t::dow>*> new_root;
  typename std::list<HGeometry<DIM,this_t::dow>*>::iterator
    the_ele = this->rootElement().begin(),
    end_ele = this->rootElement().end();
  for (;the_ele != end_ele;++ the_ele) {
    for (u_int i = 0;i < HGeometry<DIM,this_t::dow>::n_child;++ i) {
      HGeometry<DIM,this_t::dow> * p_geo = (*the_ele)->child[i];
      if (p_geo == NULL) {
        std::cerr << "Macro element is not refined. The root elements are not erased!" 
                  << std::endl;
        abort();
      }
      new_root.push_back(p_geo);
      this->nullParent(*p_geo);
    }
  }
  new_root.swap(this->rootElement());

  /**
   * 删除掉不要的对象，分为两步来进行：第一步将重复删除的指针都清零，
   * 使用 tryDeleteGeometry，第二步将剩余的指针所引用的对象删除，使用
   * deleteGeometry。在 tryDeleteGeometry 的过程中，我们顺手将几何体的
   * 共享信息在系统列表中的入口删除掉。
   */
  property_id_t<bool> pid;
  new_property_id(pid);
  end_ele = new_root.end();
  the_ele = new_root.begin();
  for (;the_ele != end_ele;++ the_ele) {
    this->tryDeleteGeometry(*the_ele, pid);
  }
  free_property_id(pid);
  the_ele = new_root.begin();
  for (;the_ele != end_ele;++ the_ele) {
    this->deleteGeometry(*the_ele);
  }
}

TEMPLATE
template <class HGEO>
void THIS::set_is_used_up(HGEO& geo) const {
  if (geo.parent != NULL) {
    set_is_used_up(*geo.parent);
  }
  u_int n_bound = geo.n_boundary;
  for (u_int i = 0;i < n_bound;++ i) {
    this->set_is_used_up(*geo.boundary[i]);
  }
  HTools tools;
  tools.setGeometryUsed(geo);
}

TEMPLATE
void THIS::set_is_used_up() {
  /**
   * 在使用了 ULoadBalance 以后，我们需要对所有使用过的几何体的祖先也
   * 设置上使用过的标志。
   */
  typename std::list<HGeometry<this_t::dim,this_t::dow>*>::iterator
    the_ele = this->rootElement().begin(),
    end_ele = this->rootElement().end();
  for (;the_ele != end_ele;++ the_ele) {
    if ((*the_ele)->parent != NULL) {
      this->set_is_used_up(*(*the_ele)->parent);
    }
  }

}

#undef THIS
#undef TEMPLATE

//////////////////////////////////////////////////////

#define TEMPLATE template <class FOREST>
#define THIS HGeometryMatcher<FOREST>

TEMPLATE
bool THIS::match_geometry() {
  typedef THIS this_t;

  std::cerr << "Rank " << _forest->rank()
            << ": Matching h-geometry among partitions ..." << std::endl;
  bool is_operated = false;
  do {
    _is_operated = false;

    if (this_t::dim >= 3)
      sync_data(_forest->communicator(), _forest->_shared_list_3d, *this, 
                &this_t::template pack_match_geometry<3>,
                &this_t::template unpack_match_geometry<3>,
                &this_t::template is_pack_match_geometry<3>);
    if (this_t::dim >= 2)
      sync_data(_forest->communicator(), _forest->_shared_list_2d, *this,
                &this_t::template pack_match_geometry<2>,
                &this_t::template unpack_match_geometry<2>,
                &this_t::template is_pack_match_geometry<2>);
    if (this_t::dim >= 1)
      sync_data(_forest->communicator(), _forest->_shared_list_1d, *this,
                &this_t::template pack_match_geometry<1>,
                &this_t::template unpack_match_geometry<1>,
                &this_t::template is_pack_match_geometry<1>);

    int src = _is_operated?1:0, dst = 0; /// 对是否做了操作进行同步
    MPI_Allreduce(&src, &dst, 1, MPI_INT, MPI_SUM, _forest->communicator());
    if ((_is_operated = (dst > 0))) {
      is_operated = true;
    }
  } while (_is_operated);
  _forest->set_shared_info_sent();

  return is_operated;
}

TEMPLATE
template <int GDIM> bool 
THIS::is_pack_match_geometry(HGeometry<GDIM,dow> * geo) const {
  u_int n_unsent_vtx = 0;
  u_int n_unsent_bnd = 0;
  if (GDIM == 1) { /**
                    * 如果维数GDIM为1，则打包没有发送过的顶点几何体的信息。
                    */

    /// 先清点需要发送的顶点个数
    u_int n_vertex = geo->n_vertex;
    for (u_int i = 0;i < n_vertex;++ i) {
      if (! is_shared_info_sent(*geo->vertex[i])) {
        n_unsent_vtx += 1;
      }
    }
  } else { /**
            * 否则，打包没有父亲的边界几何体的信息。
            */
    /// 先清点需要发送的边界个数
    u_int n_boundary = geo->n_boundary;
    for (u_int i = 0;i < n_boundary;++ i) {
      HGeometry<GDIM-1,dow> * bnd = geo->boundary[i];
      if (bnd->parent == NULL && /// 几何体没有父亲
          (! is_shared_info_sent(*bnd))) {
        n_unsent_bnd += 1;
      }
    }
  }

  u_int n_child = 0;
  if (geo->isRefined() && /// 如果几何体没有被细分，则不发送孩子信息
      (! is_shared_info_sent(*geo->child[0]))) {
    n_child = geo->n_child;
  }
  return ((n_unsent_vtx > 0) ||
          (n_unsent_bnd > 0) ||
          (n_child > 0));
}

TEMPLATE
template <int GDIM> void 
THIS::pack_match_geometry(HGeometry<GDIM,dow> * geo,
                          int remote_rank,
                          Migration::ostream<>& os) {
  Point<dow> grp; /// 存储几何体的参考点坐标

  if (GDIM == 1) { /**
                    * 如果维数GDIM为1，则打包没有发送过的顶点几何体的信息。
                    */

    /// 先清点需要发送的顶点个数
    u_int n_vertex = geo->n_vertex;
    u_int n_unsent_vtx = 0;
    for (u_int i = 0;i < n_vertex;++ i) {
      if (! is_shared_info_sent(*geo->vertex[i])) {
        n_unsent_vtx += 1;
      }
    }

    os << n_unsent_vtx; /// 输出没发送过的顶点个数
    for (u_int i = 0;i < n_vertex;++ i) {
      HGeometry<0,dow> * vtx = geo->vertex[i];
      if (! is_shared_info_sent(*vtx)) {
        grp = (*vtx); /// 取出顶点坐标
        os << vtx << grp;
      }
    }
  } else { /**
            * 否则，打包没有父亲的边界几何体的信息。
            */
    /// 先清点需要发送的边界个数
    u_int n_boundary = geo->n_boundary;
    u_int n_unsent_bnd = 0;
    for (u_int i = 0;i < n_boundary;++ i) {
      HGeometry<GDIM-1,dow> * bnd = geo->boundary[i];
      if (bnd->parent == NULL && /// 几何体没有父亲
          (! is_shared_info_sent(*bnd))) {
        n_unsent_bnd += 1;
      }
    }

    os << n_unsent_bnd; /// 输出需要发送的边界几何体的个数
    for (u_int i = 0;i < n_boundary;++ i) {
      HGeometry<GDIM-1,dow> * bnd = geo->boundary[i];
      if (bnd->parent == NULL && /// 几何体没有父亲
          (! is_shared_info_sent(*bnd))) {
        geometry_reference_point(*bnd, grp);
        os << bnd << grp;
      }
    }
  }

  /**
   * 打包发送孩子的信息。
   */
  if (geo->isRefined() && /// 如果几何体没有被细分，则不发送孩子信息
      (! is_shared_info_sent(*geo->child[0]))) {
    /**
     * 此几何体的第一个孩子信息已经发送过，则此几何体上的匹配应该已
     * 经完成。对一个几何体的孩子来说，其所有孩子的匹配必定是同时完
     * 成的。
     */

    u_int n_child = geo->n_child;
    os << n_child; /// 输出孩子的个数
    for (u_int i = 0;i < n_child;++ i) {
      HGeometry<GDIM,dow> * chd = geo->child[i]; /// 取出孩子
      geometry_reference_point(*chd, grp); /// 计算孩子参考点坐标
      os << chd << grp; /// 将孩子的信息输出到流中
    }
  } else { /// 否则输出一个0，表示没有孩子
    u_int n_child = 0;
    os << n_child;
  }
}

TEMPLATE
template <int GDIM> void 
THIS::unpack_match_geometry(HGeometry<GDIM,dow> * geo,
                            int remote_rank,
                            Migration::istream<>& is) {
  /// 计算几何体的典型尺寸
  double h = distance(*geo->vertex[0], *geo->vertex[1]);

  if (GDIM == 1) { /**
                    * 如果维数GDIM为1，则解码没有发送过的顶点几何体的信息。
                    */

    /// 首先将相应的信息从流中解出来
    u_int n_unsent_vtx = 0;
    is >> n_unsent_vtx;
    std::vector<HGeometry<0,dow>*> remote_vtx(n_unsent_vtx);
    std::vector<Point<dow> > pnt(n_unsent_vtx);
    for (u_int i = 0;i < n_unsent_vtx;++ i) {
      is >> remote_vtx[i] >> pnt[i];
    }

    /// 然后和本几何体的顶点进行匹配
    u_int n_vertex = geo->n_vertex;
    for (u_int i = 0;i < n_unsent_vtx;++ i) {
      for (u_int j = 0;j < n_vertex;++ j) {
        HGeometry<0,dow> * vtx = geo->vertex[j];
        int type = matcher().value(pnt[i], *vtx, 1.0e-4*h);
        if (type >= 0) {
          Shared_object<HGeometry<0,dow> > * p_info = get_shared_info(*vtx);
          if (p_info == NULL) p_info = new_shared_info(*vtx);
          _is_operated |= p_info->add_clone(remote_rank, type, remote_vtx[i]);
          break;
        }
      }
    }
  } else { /**
            * 否则，解码没有父亲的边界几何体的信息
            */

    /// 首先将相应的信息从流中解出来
    u_int n_unsent_bnd = 0;
    is >> n_unsent_bnd;
    std::vector<HGeometry<GDIM-1,dow>*> remote_bnd(n_unsent_bnd);
    std::vector<Point<dow> > remote_grp(n_unsent_bnd);
    for (u_int i = 0;i < n_unsent_bnd;++ i) {
      is >> remote_bnd[i] >> remote_grp[i];
    }

    /// 然后和本几何体的边界进行匹配
    Point<dow> grp;
    u_int n_boundary = geo->n_boundary;
    for (u_int i = 0;i < n_unsent_bnd;++ i) {
      for (u_int j = 0;j < n_boundary;++ j) {
        HGeometry<GDIM-1,dow> * bnd = geo->boundary[j];
        if (bnd->parent != NULL) continue;
        geometry_reference_point(*bnd, grp);
        int type = matcher().value(remote_grp[i], grp, 1.0e-4*h);
        if (type >= 0) {
          Shared_object<HGeometry<GDIM-1,dow> > * p_info = get_shared_info(*bnd);
          if (p_info == NULL) p_info = new_shared_info(*bnd);
          _is_operated |= p_info->add_clone(remote_rank, type, remote_bnd[i]);
          break;
        }
      }
    }
  }
      
  u_int n_child;
  is >> n_child;
  if (n_child != 0) { /// 如果对方发来了孩子的信息
    assert (n_child == geo->n_child); /// 确认孩子的数目正确
        
    /// 先将发来的孩子数据读掉
    std::vector<HGeometry<GDIM,dow>*> remote_chd(n_child);
    std::vector<Point<dow> > remote_grp(n_child);
    for (u_int i = 0;i < n_child;++ i) {
      is >> remote_chd[i] >> remote_grp[i];
    }

    /**
     * 如果几何体不是哑几何体，并且自身还没有加密，则对几何体做加密。哑
     * 几何体会一直保持其从归档中读入的样子，不会做任何加密。
     */
    if (! (_forest->is_dummy(*geo) || geo->isRefined())) {
      geo->refine();
    }

    /**
     * 如果几何体没有加密，则不匹配后代；否则对孩子用参考点做匹配。
     */
    if (geo->isRefined()) { 
      Point<dow> grp;
      u_int n_matched_child = 0;
      for (u_int i = 0;i < n_child;++ i) {
        HGeometry<GDIM,dow> * chd = geo->child[i];
        geometry_reference_point(*chd, grp);
        for (u_int j = 0;j < n_child;++ j) {
          if (remote_chd[j] == NULL) continue; /// 这是匹配过的位置

          /// 否则计算参考点的是相同的点，是就匹配上
          int type = matcher().value(grp, remote_grp[j], 1.0e-4*h);
          if (type >= 0) {
            Shared_object<HGeometry<GDIM,dow> > * p_info = get_shared_info(*chd); 
            if (p_info == NULL) p_info = new_shared_info(*chd);
            _is_operated |= p_info->add_clone(remote_rank, type, remote_chd[j]);
            remote_chd[j] = NULL;
            n_matched_child += 1;
            break;
          }
        }
      }
      assert (n_matched_child == n_child);
    }
  }
}

TEMPLATE
bool THIS::sync_is_used() {
  typedef THIS this_t;
  std::cerr << "Rank " << _forest->rank()
            << ": Sync is_used flag ..." << std::endl;

  _is_operated = false;
  if (this_t::dim >= 1)
    sync_data(_forest->communicator(), _forest->_shared_list_1d, *this, 
              &this_t::template pack_sync_is_used<1>,
              &this_t::template unpack_sync_is_used<1>);
  if (this_t::dim >= 2)
    sync_data(_forest->communicator(), _forest->_shared_list_2d, *this, 
              &this_t::template pack_sync_is_used<2>,
              &this_t::template unpack_sync_is_used<2>);
  if (this_t::dim >= 3)
    sync_data(_forest->communicator(), _forest->_shared_list_3d, *this, 
              &this_t::template pack_sync_is_used<3>,
              &this_t::template unpack_sync_is_used<3>);

  int src = _is_operated?1:0, dst = 0; /// 对是否做了操作进行同步
  MPI_Allreduce(&src, &dst, 1, MPI_INT, MPI_SUM, _forest->communicator());
  return (dst > 0);
}

TEMPLATE
template <int GDIM> void 
THIS::pack_sync_is_used(HGeometry<GDIM,dow> * geo,
                        int remote_rank,
                        Migration::ostream<>& os) {
  bool is_used = _tools.isGeometryUsed(*geo);
  os << is_used;
}

TEMPLATE
template <int GDIM> void 
THIS::unpack_sync_is_used(HGeometry<GDIM,dow> * geo,
                          int remote_rank,
                          Migration::istream<>& is) {
  bool is_used;
  is >> is_used;

  if (is_used) {
    if (! _tools.isGeometryUsed(*geo)) {
      _tools.setGeometryUsed(*geo);
      _is_operated = true;
    }
  }
}

#undef THIS

//////////////////////////////////////////////////////

#define THIS BirdView<FOREST>

TEMPLATE void 
THIS::semiregularize() {
  forest_t& forest = this->getForest();

  if (! forest.lock()) { // 对几何遗传树进行加锁
    std::cerr << "The hierarchy geometry tree is locked, aborting ...";
    abort();
  }

  if (forest.rank() == 0)
    std::cerr << "Semiregularizing the mesh ...  " << std::flush;
  const char * timer = "-/|\\";
  int n_element_refined = 0;
  this->prepareSemiregularize();

  int round = 0;
  bool is_operated;
  HGeometryMatcher<forest_t> matcher(forest);
  /**
   * 首先进行一次几何体匹配以及将 is_used 旗标进行同步
   */
  matcher.match_geometry();
  matcher.sync_is_used(); 
  do {
    bool flag;
    do {
      flag = false;
      this->semiregularizeHelper(flag, n_element_refined);
    } while (flag);

    if (forest.rank() == 0)
      std::cerr << "\b" << timer[round] << std::flush;
    round = (round + 1)%4;

    is_operated = false;
    if (1) { // forest.n_rank() > 1) {
      is_operated |= matcher.match_geometry();
      is_operated |= matcher.sync_is_used();
    }
  } while (is_operated);


  int n_ele_refined = 0;
  MPI_Reduce(&n_element_refined, &n_ele_refined, 1, MPI_INT, MPI_SUM, 0,
             forest.communicator());
  if (forest.rank() == 0) 
    std::cerr << "\bOK!\n"
              << "\t" << n_ele_refined 
              << " elements refined in semiregularization."
              << std::endl;
}

TEMPLATE
void THIS::eraseRootElement(u_int level) {
  for (u_int i = 0;i < level;++ i) {
    this->eraseRootElementOneLevel();
  }
}

TEMPLATE
void THIS::eraseRootElementOneLevel() {
  /// 先收集背景单元的孩子
  typedef typename base_t::container_t container_t;
  typedef typename base_t::element_t element_t;
  container_t new_root;
  typename container_t::iterator
    the_ele = this->rootElement().begin(),
    end_ele = this->rootElement().end();
  for (;the_ele != end_ele;++ the_ele) {
    for (u_int i = 0;i < element_t::n_child;++ i) {
      element_t * p_geo = (*the_ele)->child[i];
      if (p_geo == NULL) {
        std::cerr << "Macro element is not refined. The root elements CAN'T be erased!" 
                  << std::endl;
        abort();
      }
      new_root.push_back(p_geo);
      p_geo->parent = NULL;
    }
  }
  new_root.swap(this->rootElement());

  the_ele = new_root.begin();
  end_ele = new_root.end();
  for (;the_ele != end_ele;++ the_ele) {
    delete *the_ele;
  }
}

#undef THIS
#undef TEMPLATE

/**
 * end of file
 * 
 */

