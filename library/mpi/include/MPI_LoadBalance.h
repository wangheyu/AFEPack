/**
 * @file   MPI_LoadBalance.h
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct  6 18:15:44 2009
 * 
 * @brief  完成动态负载平衡的功能。
 *
 * 关于动态负载平衡的基本目标确定为：对背景网格进行重新分区，从而完成一
 * 定程度上的动态负载平衡。这个目标的实现分为下面的步骤：
 *
 * 1. 将每个背景单元上的负载进行计算；
 *
 * 2. 根据背景单元上的负载进行重新分区；
 *
 * 3. 对新的分区方式进行数据分析和准备；
 *
 * 4. 利用序列化将准备的数据存储，以此完成数据分拆；
 *
 * 5. 利用序列化将存储的数据重新在恰当的进程载入；
 *
 * 6. 对载入的数据进行分析、数据合并以及恢复；
 * 
 * 在真正的并行机上实际进行并行计算的时候，数据在存储的时候将会由网络
 * 文件系统(NFS)完成最终的硬件存储操作，这会导致非常严重的延迟。因此在
 * 这样的情况下，我们使用下面的方式来进行实用的负载平衡：
 *
 * 1. 使用 HLoadBalance 类将系列化的数据存储在计算节点的本地硬盘上；
 *
 * 2. 使用 lb_sync_local_data_dir 函数将存储的数据在进程间进行交换，使
 *    得每个进程都获得其进行载入时需要的数据文件；
 *
 * 3. 使用 HLoadBalance 类将数据从本地硬盘上载入；
 *
 * HLoadBalance 在存储的时候会将本进程的在负载平衡以前的数据存储下来，
 * 然后由 lb_sync_local_data_dir 进行分发，使得每个进程拥有负载平衡后
 * 应该载入的数据，这样就在使用网络文件系统的硬件上解决了延迟问题。由
 * 于存储在硬盘上的数据是一个系统的目录结构，对存储和载入数据的目录我
 * 们有下面两点要求：
 *
 * 1. 为了避免出现由于偶然性所带来的问题，需要将数据交换前后的目录名给
 *    一个不会碰巧被别人用的名字。与此同时，在数据网本节点的本地硬盘上
 *    进行存储的时候，请先将可能出问题的目录删除一遍，以避免冲突。
 *
 * 2. 使用中，要求整个通讯器上的存储和载入数据的目录名是一致的，这样可
 *    以使得在同一个计算节点上运行的进程不必进行数据的传输，提高程序运
 *    行的效率。否则，检测不到数据在本地有，将会进行数据传输，这样尽管
 *    运行不会出错，但是效率将会降低。
 *
 */

#ifndef __MPI_LoadBalance_h__
#define __MPI_LoadBalance_h__

#ifdef __MPI_ULoadBalance_h__
#error "MPI_LoadBalance.h 和 MPI_ULoadBalance.h 是冲突的，不能一起用。"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdio>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "../../include/Serialization.h"
#include "MPI.h"
#include "MPI_HGeometry_archive.h"
#include "MPI_Migration.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  template <class FOREST>
    class HLoadBalance {
  public:
    typedef FOREST forest_t;
    enum {dim = FOREST::dim, dow = FOREST::dow};
    typedef typename forest_t::matcher_t matcher_t;
  private:
  // 在计算每个几何体上的负载是使用的性质
  property_id_t<double> _pid_loading;

  // 这个性质的 first 为新分区中的秩， second 为对其进行存储的秩                    
  property_id_t<std::map<int, int> > _pid_rank_map;

  // 部分需要拆分和合并的几何体的全局标号
  property_id_t<unsigned long> _pid_global_idx;

  // 元素为(全局标号，对象指针)的对
  std::map<unsigned long, void *> _global_pointer_to_merge;

  // (进程，列表(全局标号，维数，对象指针)组)的映射
  typedef std::map<unsigned long, std::pair<int, void *> > tuple_map_t;
  std::map<int, tuple_map_t> _global_pointer_to_share;

  typedef BirdView<forest_t> birdview_t;
  typedef BirdViewSet<forest_t> birdview_set_t;
  typedef HLoadBalance<forest_t> this_t;

  forest_t * _forest; /// 树结构
  std::vector<u_int> _cut_point; /// 重新分区的拆分点
  std::vector<int> _new_rank; /// 新分区的秩

  public:
  HLoadBalance(forest_t& forest) : _forest(&forest) {}
  ~HLoadBalance() {}

  private:
  template <class GEO> std::map<int,int> *
  new_rank_map(GEO& geo) const {
    return geo.new_property(_pid_rank_map);
  }
  template <class GEO> unsigned long *
  new_global_idx(const GEO& geo) const {
    return geo.new_property(_pid_global_idx);
  }

  public:
  void clear() {
    free_property_id(_pid_loading);
    free_property_id(_pid_rank_map);
    free_property_id(_pid_global_idx);
    _global_pointer_to_merge.clear();
    _global_pointer_to_share.clear();
    _cut_point.clear();
    _new_rank.clear();
  }

  template <class GEO> std::map<int,int> *
  get_rank_map(GEO& geo) const {
    return geo.get_property(_pid_rank_map);
  }
  template <class GEO> unsigned long *
  get_global_idx(const GEO& geo) const {
    return geo.get_property(_pid_global_idx);
  }

  /**
   * 判断一个几何体在新分区上的存储是否由本进程完成。
   * 
   * @param map 几何体的 rank_map 性质指针
   * @param new_rank 几何体在新分区上的秩
   * 
   * @return 在本进程存储为真，否则为假 
   */
  bool is_save_on_this_rank(std::map<int,int> * map, 
                            int new_rank) {
    typename std::map<int,int>::iterator
    the_pair = map->find(new_rank);
    /// 如果几何体的 rank_map 性质中应该不会没有 new_rank 的入口
    assert (the_pair != map->end());

    /// 如果几何体对 new_rank 分区的存储在本进程上完成，返回真
    return (the_pair->second == _forest->rank()); /// 否则返回假
  }

  /** 
   * 一个几何体的多个备份归并到一个秩上时，完成其相应指针的归并。存储
   * 的时候，因为一个几何体的多个备份事实上只存储了一次，其它分区上的
   * 这个指针在载入的时候将会仅仅获得全局标号。因此，如果该几何体事实
   * 上存储在此归档中，我们会记录下该指针的位置；如果该几何体不存储于
   * 此归档中，我们直接取出该指针。
   * 
   * @param type 类型(0: 存储在此文档；1: 不存储在此文档中)
   * @param global_idx 全局标号
   * @param p_geo 几何体指针的引用
   */
  template <class GEO>
  void merge_global_pointer(int type,
                            unsigned long global_idx,
                            GEO *& p_geo) {
    std::map<unsigned long, void *>::iterator
    the_pair = _global_pointer_to_merge.find(global_idx);
    if (type == 2) {
      if (the_pair == _global_pointer_to_merge.end()) {
        _global_pointer_to_merge.insert(std::pair<unsigned long,void*>(global_idx, p_geo));
      } else {
        assert (the_pair->second == p_geo);
      }
    } else {
      assert (type == 4 && the_pair != _global_pointer_to_merge.end());
      p_geo = (GEO *)(the_pair->second);
    }
  }

  /** 
   * 记录一个指针将被共享的信息。其将会被共享给秩为 rank 的进程，全局指
   * 标为 global_idx，几何体的指针为 p_geo。
   * 
   * @param rank 去给共享的秩
   * @param global_idx 全局标号
   * @param p_geo 几何体的指针
   */
  template <class GEO>
  void share_global_pointer(int rank, 
                            unsigned long global_idx,
                            GEO *& p_geo) {
    _global_pointer_to_share[rank][global_idx] = 
    std::pair<int,void *>(p_geo->dimension, p_geo);
  }

  private:
  /** 
   * 在进程之间交换共享几何体的信息，建立共享信息表。
   */
  void share_global_pointer() {
    std::cerr << "Exchanging shared global pointers ..." << std::endl;

    typedef std::map<int, tuple_map_t> map_t;

    int n = 0;
    std::list<int> target;
    std::list<BinaryBuffer<> > out_buf, in_buf;
    typename map_t::iterator
    the_pair = _global_pointer_to_share.begin(),
    end_pair = _global_pointer_to_share.end();
    for (;the_pair != end_pair;++ the_pair, ++ n) {
      target.push_back(the_pair->first);

      out_buf.push_back(BinaryBuffer<>());
      Migration::ostream<> os(out_buf.back());
      
      in_buf.push_back(BinaryBuffer<>());

      u_int n_item = the_pair->second.size();
      os << n_item;
      std::cerr << "sending " << the_pair->second.size()
                << " items from " << _forest->rank()
                << " to " << the_pair->first << std::endl;
      typename tuple_map_t::iterator
        the_tuple = the_pair->second.begin(),
        end_tuple = the_pair->second.end();
      for (;the_tuple != end_tuple;++ the_tuple) {
        os << the_tuple->first /// global index
           << the_tuple->second.first /// dimension
           << the_tuple->second.second; /// pointer
      }
    }

    MPI_Barrier(_forest->communicator());
    sendrecv_data(_forest->communicator(), n, out_buf.begin(), 
                  in_buf.begin(), target.begin());

    typename std::list<int>::iterator 
    the_rank = target.begin(),
    end_rank = target.end();
    typename std::list<BinaryBuffer<> >::iterator 
    the_buf = in_buf.begin();
    for (;the_rank != end_rank;++ the_rank, ++ the_buf) {
      Migration::istream<> is(*the_buf);
      u_int n_item;
      is >> n_item;
      std::cerr << "received " << n_item
                << " items from " << *the_rank
                << " to " << _forest->rank() << std::endl;
      for (u_int i = 0;i < n_item;++ i) {
        unsigned long global_idx;
        int dimension;
        void * remote_obj;
        is >> global_idx >> dimension >> remote_obj;
        std::pair<int,void *>& the_pair = 
          _global_pointer_to_share[*the_rank][global_idx];
        assert (dimension == the_pair.first && 
                dimension >= 0 && dimension <= 3);

        void * local_obj = the_pair.second;
        if (dimension == 0) {
#define GDIM 0
#define GEO HGeometry<GDIM,dow>
#define OBJ Shared_object<GEO>
          GEO * p_geo = (GEO *)local_obj;
          OBJ * p_info = _forest->get_shared_info(*p_geo);
          if (p_info == NULL) p_info = _forest->new_shared_info(*p_geo);
          p_info->add_clone(*the_rank, (GEO *)remote_obj);
#undef OBJ
#undef GEO
#undef GDIM
        } else if (dimension == 1) {
#define GDIM 1
#define GEO HGeometry<GDIM,dow>
#define OBJ Shared_object<GEO>
          GEO * p_geo = (GEO *)local_obj;
          OBJ * p_info = _forest->get_shared_info(*p_geo);
          if (p_info == NULL) p_info = _forest->new_shared_info(*p_geo);
          p_info->add_clone(*the_rank, (GEO *)remote_obj);
#undef OBJ
#undef GEO
#undef GDIM
        } else if (dimension == 2) {
#define GDIM 2
#define GEO HGeometry<GDIM,dow>
#define OBJ Shared_object<GEO>
          GEO * p_geo = (GEO *)local_obj;
          OBJ * p_info = _forest->get_shared_info(*p_geo);
          if (p_info == NULL) p_info = _forest->new_shared_info(*p_geo);
          p_info->add_clone(*the_rank, (GEO *)remote_obj);
#undef OBJ
#undef GEO
#undef GDIM
        } else if (dimension == 3) {
#define GDIM 3
#define GEO HGeometry<GDIM,dow>
#define OBJ Shared_object<GEO>
          GEO * p_geo = (GEO *)local_obj;
          OBJ * p_info = _forest->get_shared_info(*p_geo);
          if (p_info == NULL) p_info = _forest->new_shared_info(*p_geo);
          p_info->add_clone(*the_rank, (GEO *)remote_obj);
#undef OBJ
#undef GEO
#undef GDIM
        }
      }
    }
  }

  //@{
  public:
  /**
   * config 计算在基于几何遗传树上的一个正则网格上的每个单元上的负载量，
   * 计算的方法由用户提供。此函数有两个接口，一个是类的成员函数形式提供
   * 计算方法，一个是由一个普通函数提供计算方法。
   */
  template <class LOADER>
  void config(birdview_t& ir_mesh,
              LOADER& loader,
              double (LOADER::*value)(GeometryBM&)) {

    if (&(ir_mesh.getForest()) != _forest) {
      std::cerr << "The mesh should be on the same forest of the HLoadBalance." 
                << std::endl;
      abort();
    }

    RegularMesh<dim,dow>& mesh = ir_mesh.regularMesh();
    new_property_id(_pid_loading); 
    u_int n_ele = mesh.n_geometry(dim);
    for (u_int i = 0;i < n_ele;++ i) {
      HGeometry<dim,dow> * p_geo = mesh.template h_geometry<dim>(i);
      double * p_loading = p_geo->get_property(_pid_loading);
      if (p_loading == NULL) {
        p_loading = p_geo->new_property(_pid_loading);
      }
      (*p_loading) = (loader.*value)(mesh.geometry(dim,i));
    }
  }

  void config(birdview_t& ir_mesh,
              double (*value)(GeometryBM&) = &_default_element_loading) {
    if (&(ir_mesh.getForest()) != _forest) {
      std::cerr << "The mesh should be on the same forest of the HLoadBalance." 
                << std::endl;
      abort();
    }

    RegularMesh<dim,dow>& mesh = ir_mesh.regularMesh();
    new_property_id(_pid_loading); 
    u_int n_ele = mesh.n_geometry(dim);
    for (u_int i = 0;i < n_ele;++ i) {
      const HGeometry<dim,dow> * p_geo = mesh.template h_geometry<dim>(i);
      double * p_loading = p_geo->get_property(_pid_loading);
      if (p_loading == NULL) {
        p_loading = p_geo->new_property(_pid_loading);
      }
      (*p_loading) = (*value)(mesh.geometry(dim,i));
    }
  }

  static double _default_element_loading(GeometryBM&) { return 1.0; }
  //@}

  //@{
  public:
  /**
   * 对区域进行重新拆分，获取使得负载平衡的拆分点位置。这个功能由
   * partition, lump_loading, get_loading 三个函数完成。
   */
  void partition(u_int n_new_rank = 0) {
    int rank =  _forest->rank();
    int n_rank = _forest->n_rank();

    double loading = this->lump_loading();
    double total_loading, partial_loading;
    MPI_Comm comm = _forest->communicator();
    MPI_Scan(&loading, &partial_loading, 1, MPI_DOUBLE, MPI_SUM, comm);

    if (rank == n_rank - 1) {
      total_loading = partial_loading;
    }
    if (n_new_rank == 0) n_new_rank = n_rank; /// 缺省取为旧的分区数

    MPI_Bcast(&total_loading, 1, MPI_UNSIGNED, n_rank - 1, comm);
    double mean_loading = n_new_rank;
    mean_loading = total_loading/mean_loading;

    _cut_point.clear();
    _new_rank.clear();

    u_int idx = 0;
    _cut_point.push_back(idx); /// put a "0" at beginning

    loading  = partial_loading - loading;
    u_int current_rank = (u_int)(loading/mean_loading);
    _new_rank.push_back(current_rank);

    typename forest_t::RootIterator
    the_ele = _forest->beginRootElement(),
    end_ele = _forest->endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      idx += 1;
      loading += this->get_loading(*the_ele);
      u_int the_rank = (u_int)(loading/mean_loading);
      the_rank = std::min(the_rank, n_new_rank - 1); /// 最后一个进程上会多几个单元
      if (the_rank > current_rank) {
        _cut_point.push_back(idx);
        _new_rank.push_back(++ current_rank);
      }
    }
    if (_cut_point.back() != idx) {
      _cut_point.push_back(idx);
    } else {
      _new_rank.pop_back();
    }

    free_property_id(_pid_loading); 
  }

  /**
   * 使用原有分区。此功能可以用来存储和重新载入数据。
   */
  void reuse_partition() {
    _cut_point.clear();
    _new_rank.clear();

    _cut_point.resize(2, 0);
    _cut_point[1] = _forest->n_rootElement();
    _new_rank.resize(1, _forest->rank());
  }

  private:
  /**
   * lump_loading 函数利用 get_loading 函数的对后代递归的性质将正则网格
   * 的单元上的负载集中到背景网格的单元上。
   */
  double lump_loading() const {
    double loading = 0.0;
    typename forest_t::RootIterator
    the_ele = _forest->beginRootElement(),
    end_ele = _forest->endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      loading += this->get_loading(*the_ele);
    }
    return loading;
  }

  /**
   * 计算几何体 geo 上的负载。此函数将计算所有后代的负载，并将其累加到
   * 自身的负载上。
   */
  template <class GEO> double
  get_loading(GEO& geo) const {
    double * p_loading = geo.get_property(_pid_loading);
    if (p_loading == NULL) {
      p_loading = geo.new_property(_pid_loading);
      for (u_int i = 0;i < geo.n_child;++ i) {
        (*p_loading) += get_loading(*geo.child[i]);
      }
    }
    return (*p_loading);
  }
  //@}

  //@{
  private:
  /**
   * 对整个网格的几何体进行标识，判断几何体在负载平衡完成后每个几何体所
   * 处的分区。这个功能由 set_new_rank, pack_collect_rank,
   * unpack_collect_rank, geometry_set_new_rank 四个函数完成。
   */
  void set_new_rank() {
    new_property_id(_pid_rank_map);

    /**
     * 在每个分区上计算每个背景网格中的几何体所在分区。
     */
    typename forest_t::RootIterator
    the_ele = _forest->beginRootElement();
    u_int n_part = _new_rank.size();
    for (u_int i = 0;i < n_part;++ i) {
      for (u_int j = _cut_point[i];j < _cut_point[i + 1];++ j, ++ the_ele) {
        geometry_set_new_rank(*the_ele, _new_rank[i]);
      }
    }

    if (_forest->n_rank() <= 1) return;
    /**
     * 将各个分区上所加的分区标识进行汇总。
     */
    Shared_type_filter::only<0> type_filter;
    MPI_Comm comm = _forest->communicator();
#define SYNC_DATA(D) \
    if (dim >= D) { \
      sync_data(comm, _forest->template get_shared_list<D>(), *this, \
                &this_t::template pack_set_new_rank<D>, \
                &this_t::template unpack_set_new_rank<D>, \
                type_filter); \
    }
    SYNC_DATA(0);
    SYNC_DATA(1);
    SYNC_DATA(2);
    SYNC_DATA(3);
#undef SYNC_DATA
  }

  public:
  template <int GDIM> void
  pack_set_new_rank(HGeometry<GDIM,dow> * geo,
                    int remote_rank,
                    Migration::ostream<>& os) {
    std::map<int,int> * p_map = get_rank_map(*geo);
    assert (p_map != NULL);

    os << p_map->size();
    typename std::map<int,int>::iterator
    the_rank = p_map->begin(),
    end_rank = p_map->end();
    for (;the_rank != end_rank;++ the_rank) {
      os << the_rank->first << the_rank->second;
    }
  }

  template <int GDIM> void
  unpack_set_new_rank(HGeometry<GDIM,dow> * geo,
                      int remote_rank,
                      Migration::istream<>& is) {
    /**
     * 设置对每个几何体进行输出的最小秩。输出的时候，几何体跟随单元进行
     * 输出，因此只有当一个几何体所属于的单元，在某个进程上输出的时候，
     * 该几何体才会被输出。对于非单元的几何体，我们就将对其进行输出的最
     * 小秩进程作为最后对其进行输出的进程。
     */
    std::map<int,int> * p_map = get_rank_map(*geo);
    assert (p_map != NULL);

    int new_rank, old_rank;
    std::size_t n;
    is >> n;
    for (std::size_t i = 0;i < n;++ i) {
      is >> new_rank >> old_rank;
      typename std::map<int,int>::iterator 
        the_pair = p_map->find(new_rank);
      if (the_pair == p_map->end()) {
        p_map->insert(std::pair<int,int>(new_rank, old_rank));
      } else {
        /// 取较小的输出秩
        if (the_pair->second > old_rank) {
          the_pair->second = old_rank;
        }
      }
    }
  }

  private:
  /**
   * 对几何体标记为将位于新分区的第 new_rank 个进程上，并且由本进程进行
   * 存储。此函数对其所有的顶点、边界和后代做递归。
   */
  template <class GEO> void
  geometry_set_new_rank(GEO& geo, int new_rank) const {
    std::map<int,int> * p_map = get_rank_map(geo);
    if (p_map == NULL) p_map = new_rank_map(geo);
    p_map->insert(std::pair<int,int>(new_rank, _forest->rank())); /// 将存储秩设为当前的秩
    
    /// 然后对顶点、边界和后代做递归    
    for (u_int i = 0;i < geo.n_vertex;++ i) {
      geometry_set_new_rank(*geo.vertex[i], new_rank);
    }
    for (u_int i = 0;i < geo.n_boundary;++ i) {
      geometry_set_new_rank(*geo.boundary[i], new_rank);
    }
    if (geo.isRefined()) {
      for (u_int i = 0;i < geo.n_child;++ i) {
        geometry_set_new_rank(*geo.child[i], new_rank);
      }
    }
  }
  //@}

  //@{
  private:
  /**
   * 对现在处于共享状态的所有几何体进行自洽标号，使得几何体在存储后恢复
   * 的时候能够被正确合并。几何体的合并是一个非常困难的操作，目前设想的
   * 实现方式如下：对于每一个在整个分布式内存中具有多个拷贝的对象，所在
   * 的m个进程的秩为(r1, ..., rm)，假设在重分区以后，这个对象将会在n个
   * 进程上被恢复，这n个进程的秩为(R1, ..., Rn)。其中，在Rk上恢复出来的
   * 这个对象，将会被进行输出的m个进程中的一部分进程输出，假定为(r1,
   * ..., rl)这l个进程。在实现过程中，我们使得(r1, ..., rl)这l个进程事
   * 实上只是在进程号最小的进程上输出一个备份，并对自身记录一个标号，其
   * 它l-1进程对此对象的指针进行输出的时候，则输出这个标号，这样在恢复
   * 的时候就获得唯一的备份指针。
   *
   * 实现这个全局自洽标号的功能由 global_index, pack_global_index,
   * unpack_global_index, geometry_global_index 四个函数完成。
   */
  void global_index() {
    /**
     * 先清点需要加全局标号的几何体的个数。
     */
    new_property_id(_pid_global_idx);
    unsigned long idx = 0;
    typename forest_t::RootIterator
    the_ele = _forest->beginRootElement(),
    end_ele = _forest->endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      geometry_global_index(*the_ele, idx);
    }
    free_property_id(_pid_global_idx); /// 对从0开始的标号进行清除

    /**
     * 计算本进程上全局标号的开始值。
     */
    MPI_Comm comm = _forest->communicator();
    unsigned long global_idx = idx;
    MPI_Scan(&idx, &global_idx, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
    idx = global_idx - idx;

    /**
     * 从 idx 开始重新对几何体进行标号，完成自洽标号的指定。
     */
    new_property_id(_pid_global_idx);
    the_ele = _forest->beginRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      geometry_global_index(*the_ele, idx);
    }
    
    if (_forest->n_rank() <= 1) return;
    /**
     * 将各个分区上具有多个拷贝的几何体的变号进行同步，我们取那个最小
     * 的标号。
     */
    Shared_type_filter::only<0> type_filter;
#define SYNC_DATA(D) \
    if (dim >= D) { \
      sync_data(comm, _forest->template get_shared_list<D>(), *this, \
                &this_t::template pack_global_index<D>, \
                &this_t::template unpack_global_index<D>, \
                type_filter); \
    }
    SYNC_DATA(0);
    SYNC_DATA(1);
    SYNC_DATA(2);
    SYNC_DATA(3);
#undef SYNC_DATA
  }

  public:
  template <int GDIM> void
  pack_global_index(HGeometry<GDIM,dow> * geo,
                    int remote_rank,
                    Migration::ostream<>& os) {
    unsigned long * p_idx = get_global_idx(*geo);
    assert (p_idx != NULL);

    os << *p_idx; /// 传输自身的原始全局标号
  }

  /**
   * 对于在目前分区状态下已经有共享的几何体，其最终的全局标号取为其所
   * 有备份上的预设全局标号的最小值。这使得所有具有全局标号的几何体的
   * 全局标号并不成为一个连续的序列，但是是自洽的。
   */
  template <int GDIM> void
  unpack_global_index(HGeometry<GDIM,dow> * geo,
                      int remote_rank,
                      Migration::istream<>& is) {
    unsigned long * p_idx = get_global_idx(*geo);
    assert (p_idx != NULL);

    unsigned long remote_idx;
    is >> remote_idx;
    if (*p_idx > remote_idx) {
      *p_idx = remote_idx; /// 将标号取为所有拷贝中标号最小的
    }
  }

  private:
  /**
   * 对已经共享的几何体和将会共享的几何体进行全局标号。
   */
  template <class GEO> void
  geometry_global_index(GEO& geo, 
                        unsigned long& idx) const {
    unsigned long * p_idx = get_global_idx(geo);
    if (p_idx != NULL) return; /// 跳过已经处理过的几何体

    std::map<int,int> * p_map = get_rank_map(geo);
    if (_forest->get_shared_info(geo) != NULL || /// 原本是共享的几何体
        p_map->size() > 1) { /// 或者新秩个数大于 1 的情形都编号
      p_idx = new_global_idx(geo);
      *p_idx = idx ++;
    }

    /// 然后对顶点、边界和后代做递归    
    for (u_int i = 0;i < geo.n_vertex;++ i) {
      geometry_global_index(*geo.vertex[i], idx);
    }
    for (u_int i = 0;i < geo.n_boundary;++ i) {
      geometry_global_index(*geo.boundary[i], idx);
    }
    if (geo.isRefined()) {
      for (u_int i = 0;i < geo.n_child;++ i) {
        geometry_global_index(*geo.child[i], idx);
      }
    }
  }
  //@}

  //@{
  public:
  void write_config_file(const std::string& dirname) {
    char filename[1024];
    sprintf(filename, "%s/.config", dirname.c_str());
    std::ofstream os(filename);
    os << _forest->n_rank() << "\t# number of old rank" << std::endl;
    os.close();
    Migration::save_config(dirname);
  }

  private:
  void export_birdview(birdview_t& ir_mesh, 
                       Migration::data_id_t data_id, 
                       u_int idx) {
    typename birdview_t::ActiveIterator
    the_ele = ir_mesh.beginActiveElement(),
    end_ele = ir_mesh.endActiveElement();
    for (;the_ele != end_ele;++ the_ele) {
      HGeometry<dim,dow> * p_geo = the_ele->h_element;
      Migration::ostream<> os(p_geo->buffer[data_id]);
      os << idx;
    }
  }

  public:
  void save_data(const std::string& dirname, 
                 birdview_set_t& bvs) {
    set_new_rank();
    global_index();

    char command[1024], *filename = command;

    Migration::data_id_t id = Migration::name_to_id("__internal_birdviewflag__"); 
    if (! Migration::is_valid(id)) {
      id = Migration::register_data_name("__internal_birdviewflag__");
    }
    typename birdview_set_t::iterator
    the_bv = bvs.begin(), end_bv = bvs.end();
    for (u_int idx = 0;the_bv != end_bv;++ the_bv) {
      export_birdview(*the_bv, id, idx ++);
    }

    typedef HGeometry_oarchive<this_t> archive_t;
    typename forest_t::RootIterator
    the_ele = _forest->beginRootElement();
    u_int n_part = _new_rank.size();
    for (u_int i = 0;i < n_part;++ i) {
      sprintf(command, "mkdir -p %s && mkdir -p %s/%d", 
              dirname.c_str(), dirname.c_str(), _new_rank[i]);
      system(command);
      sprintf(filename, "%s/%d/%d.dat", dirname.c_str(), 
              _new_rank[i], _forest->rank());
      std::ofstream os(filename, std::ios::binary);

      {
        archive_t oa(*this, _new_rank[i], os);

        u_int n_ele = _cut_point[i + 1] - _cut_point[i];
        oa & boost::serialization::make_binary_object(&n_ele, sizeof(u_int));
        std::cerr << "saving data in " << filename 
                  << ", # macro element = " << n_ele << std::endl;

        for (u_int j = _cut_point[i];j < _cut_point[i + 1];++ j) {
          HGeometry<dim,dow> * p_ele = *(the_ele ++);
          oa & p_ele;
        }
      }

      os.close();
    }

    if (_forest->rank() == 0) write_config_file(dirname);
    MPI_Barrier(_forest->communicator()); /// 最后做一次同步
  }
  //@}

  //@{
  public:
  int load_config_file(const std::string& dirname) {
    Migration::load_config(dirname);

    char filename[1024];
    sprintf(filename, "%s/.config", dirname.c_str());
    filtering_istream is;
    Migration::ensured_open_filtering_stream(filename, is);

    int old_rank;
    is >> old_rank;
    return old_rank;
  }

  private:
  void reconstruct_birdview(HElement<dim,dow>& ele,
                            Migration::data_id_t data_id,
                            u_int idx) const {
    HGeometry<dim,dow> * p_geo = ele.h_element;
    typename Migration::buffer_t::iterator it = p_geo->buffer.find(data_id);
    bool is_active = true;
    ele.value = 0;
    if (it == p_geo->buffer.end()) {
      is_active = false;
    } else {
      BinaryBuffer<>& buf = it->second;
      int n = buf.size()/sizeof(u_int), i;
      u_int idx1;
      Migration::istream<> is(buf);
      for (i = 0;i < n;++ i) {
        is >> idx1;
        if (idx1 == idx) break;
      }
      if (i == n) is_active = false;
    }
    if (! is_active) {
      ele.value = 1;
      ele.refine();
      for (int i = 0;i < ele.n_child;++ i) {
        reconstruct_birdview(*ele.child[i], data_id, idx);
      }
    }
  }

  void reconstruct_birdview(birdview_t& ir_mesh,
                            Migration::data_id_t data_id,
                            u_int idx) const {
    typename birdview_t::RootIterator
    the_ele = ir_mesh.beginRootElement(),
    end_ele = ir_mesh.endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      reconstruct_birdview(*the_ele, data_id, idx);
    }
  }

  /**
   * 对几何体的 shared_info_sent 旗标进行设置。如果一个几何体自身是被
   * 加密了的，则我们设置此旗标，这是因为在 share_global_pointer 的时
   * 候，此信息已经进行了进程间的同步。
   */
  template <class GEO>
    void set_shared_info_sent(GEO& geo) const {
    if ((_forest->get_shared_info(geo) != NULL) &
        (! _forest->is_shared_info_sent(geo))) {
      _forest->set_shared_info_sent(geo);
    }

    for (u_int i = 0;i < geo.n_vertex;++ i) {
      if (! _forest->is_shared_info_sent(*geo.vertex[i])) {
        _forest->set_shared_info_sent(*geo.vertex[i]);
      }
    }
    for (u_int i = 0;i < geo.n_boundary;++ i) {
      this->set_shared_info_sent(*geo.boundary[i]);
    }
    if (geo.isRefined()) {
      for (u_int i = 0;i < geo.n_child;++ i) {
        this->set_shared_info_sent(*geo.child[i]);
      }
    }
  }

  /**
   * 对所有的根单元，设置其所有几何体的 shared_info_sent 旗标。
   */
  void set_shared_info_sent() const {
    typename forest_t::RootIterator
      the_ele = _forest->beginRootElement(),
      end_ele = _forest->endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      this->set_shared_info_sent(*the_ele);
    }
  }

  public:
  void load_data(const std::string& dirname,
                 birdview_set_t& bvs) {
    clear();

    int n_old_rank = load_config_file(dirname);

    /// 开始正式读数据以前同步一次，以保证要读入的数据已经准备好
    MPI_Barrier(_forest->communicator()); 

    /// 建立目录
    char filename[1024];

    typedef HGeometry_iarchive<this_t> archive_t;

    int rank = 0;
    do {
      std::ifstream is;
      do {
        sprintf(filename, "%s/%d/%d.dat", dirname.c_str(), 
                _forest->rank(), rank ++);
        if (rank > n_old_rank) {
          share_global_pointer();
          set_shared_info_sent();

          /// 重造半正则网格
          Migration::data_id_t id = Migration::name_to_id("__internal_birdviewflag__"); 
          typename birdview_set_t::iterator
            the_bv = bvs.begin(), end_bv = bvs.end();
          for (u_int idx = 0;the_bv != end_bv;++ the_bv) {
            the_bv->reinit(*_forest);
            reconstruct_birdview(*the_bv, id, idx ++);
          }
          return;
        }
        is.open(filename, std::ios::binary); /// 则打开文件为流
        if (is.good()) break; /// 打开正确，则开始读入
      } while (1);

      {
        u_int n_ele;
        std::cerr << "loading data from " << filename;
        archive_t ia(*this, _forest->rank(), is); /// 将流做成输入档案
        ia & boost::serialization::make_binary_object(&n_ele, sizeof(u_int));
        std::cerr << ", # macro element = " << n_ele << std::endl;

        for (u_int i = 0;i < n_ele;++ i) {
          HGeometry<dim,dow> * p_ele;
          ia & p_ele;
          _forest->rootElement().push_back(p_ele);
        }
      }

      is.close();
    } while (1);

  }

  };

  template <class FOREST>
    void load_forest(const std::string& dirname,
                     FOREST& forest,
                     bool is_has_orphans = false) {
    HLoadBalance<FOREST> hlb(forest);
    BirdViewSet<FOREST> bvs;
    hlb.load_data(dirname, bvs);
  }

  template <class FOREST>
    void load_mesh(const std::string& dirname,
                   FOREST& forest,
                   BirdView<FOREST>& mesh,
                   bool is_has_orphans = false) {
    HLoadBalance<FOREST> hlb(forest);
    BirdViewSet<FOREST> bvs;
    bvs.add(mesh);
    hlb.load_data(dirname, bvs);
  }

  template <class FOREST>
    void load_mesh(const std::string& dirname,
                   FOREST& forest,
                   u_int n_mesh, ...) {
    HLoadBalance<FOREST> hlb(forest);
    typedef BirdView<FOREST> * mesh_ptr_t;

    BirdViewSet<FOREST> bvs;
    va_list ap;
    va_start(ap, n_mesh);
    for (u_int i = 0;i < n_mesh;++ i) {
      mesh_ptr_t p_mesh = va_arg(ap, mesh_ptr_t);
      bvs.add(*p_mesh);
    }
    va_end(ap);

    hlb.load_data(dirname, bvs);
  }

  template <class FOREST>
    void load_mesh_set(const std::string& dirname,
                       FOREST& forest,
                       BirdViewSet<FOREST>& bvs,
                       bool is_has_orphans = false) {
    HLoadBalance<FOREST> hlb(forest);
    hlb.load_data(dirname, bvs);
  }

  /**
   * 假定在数据存储的时候，所有数据都存储在每个计算节点的局部存储上。本
   * 函数将所有这些文件都收集到网络文件系统上的同一个目录中。此实现中直
   * 接使用shell命令来完成，数据传输由网络文件系统来完成，事实上不可能
   * 做到实时，而且可能有很久的延时。其中 src_dir 应该是在计算节点的局
   * 部存储上， dst_dir 应该在全局存储上。
   */
  void lb_collect_local_data_dir(MPI_Comm comm,
                                 const std::string& src_dir,
                                 const std::string& dst_dir);
  /**
   * 在做负载平衡的时候，假定数据被存储在每个计算节点的局部存储上。我
   * 们使用下面的函数将进程之间的文件进行交换，使得每个进程都拿到其需
   * 要载入的数据文件。
   */
  void lb_sync_local_data_dir(MPI_Comm comm,
                              const std::string& src_dir,
                              const std::string& dst_dir);
}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_LoadBalance_h__

/**
 * end of file
 * 
 */
