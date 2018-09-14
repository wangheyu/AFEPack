/**
 * @file   MPI_ULoadBalance.h
 * @author Ruo Li <rli@math.pku.edu.cn>
 * @date   Sun Dec  5 21:00:16 2010
 * 
 * @brief  完成动态负载平衡的功能。此处的功能是为了使得分区能够不仅仅
 *         在最初的背景单元上进行，使用的技术途径是将几何遗传树的根单
 *         元几何体进行了更换。
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

#ifndef __MPI_ULoadBalance_h__
#define __MPI_ULoadBalance_h__

#ifdef __MPI_LoadBalance_h__
#error "MPI_LoadBalance.h 和 MPI_ULoadBalance.h 是冲突的，不能一起用。"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdio>

#include <boost/program_options.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "../../include/Serialization.h"
#include "MPI.h"
#include "MPI_UGeometry_archive.h"
#include "MPI_Migration.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  template <class FOREST>
    class HLoadBalance {
  public:
    typedef FOREST forest_t;
    enum {dim = FOREST::dim, dow = FOREST::dow};
    typedef typename forest_t::matcher_t matcher_t;
    
    struct rank_map_val {
      int old_rank; /// 对几何体进行存储的秩
      int priority; /**
                     * 存储的优先级：0 哑几何体；1 一般几何体；2 在新秩
                     * 上加密的几何体。值越大的优先级越优先存储。当两个
                     * 拷贝的优先级相同的时候，秩比较小的进程具有优先级。
                     */
    };
    typedef std::map<int,rank_map_val> rank_map_t;

    template <class GEO> bool 
      is_refined_on_new_rank(const GEO& geo, 
                             int new_rank) const {
      bool result = false;
      if (geo.isRefined() &&
          is_on_this_new_rank(*geo.child[0], new_rank)) {
        result = true;
      }
      return result;
    }

  private:
    // 在计算每个几何体上的负载是使用的性质
    property_id_t<double> _pid_loading;
    bool _is_configed; // 是否设置了负载
    property_id_t<> _pid_leaf_loading; // 是否是叶子

    // 标记一个共享几何体是否已经存储了的性质
    property_id_t<> _pid_is_saved;
  public:
    template <class GEO> void
      set_is_saved(const GEO& geo) const {
      if (geo.get_property(_pid_is_saved) == NULL) {
	geo.new_property(_pid_is_saved);
      }
    }

  private:
    /**
     * 这个性质的 first 为新分区中的秩， second 为对其进行存储的秩和在新
     * 分区上是否进行了加密。
     */
    property_id_t<rank_map_t> _pid_rank_map;

    // 部分需要拆分和合并的几何体的全局标号
    property_id_t<unsigned long> _pid_global_idx;

    // 元素为(全局标号，对象指针)的对
    std::map<unsigned long, void *> _global_pointer_to_merge;
    std::map<unsigned long, std::list<void *> > _global_pointer_to_merge_later;

    // (进程，列表(全局标号，维数，对象指针)组)的映射
    typedef std::map<unsigned long, std::pair<int, void *> > tuple_map_t;
    std::map<int, tuple_map_t> _global_pointer_to_share;

    typedef BirdView<forest_t> birdview_t;
    typedef BirdViewSet<forest_t> birdview_set_t;
    typedef HLoadBalance<forest_t> this_t;

    forest_t * _forest; /// 树结构
    std::vector<u_int> _cut_point; /// 重新分区的拆分点
    std::vector<int> _new_rank; /// 新分区的秩

    typename forest_t::container_t _nre; /// 存储时候的根单元

  public:
    HLoadBalance() : _is_configed(false) {}
      HLoadBalance(forest_t& forest) : _is_configed(false), _forest(&forest) {}
	~HLoadBalance() {}

	forest_t& forest() const { return *_forest; }
	int rank() const { return _forest->rank(); }

  private:
	template <class GEO> rank_map_t *
	  new_rank_map(const GEO& geo) const {
	  return geo.new_property(_pid_rank_map);
	}
	template <class GEO> unsigned long *
	  new_global_idx(const GEO& geo) const {
	  return geo.new_property(_pid_global_idx);
	}

  public:
	void clear() {
	  free_property_id(_pid_loading);
	  _is_configed = false;
	  free_property_id(_pid_leaf_loading);
	  free_property_id(_pid_rank_map);
	  free_property_id(_pid_global_idx);
	  _global_pointer_to_merge.clear();
	  _global_pointer_to_merge_later.clear();
	  _global_pointer_to_share.clear();
	  _cut_point.clear();
	  _new_rank.clear();
	}

	template <class GEO> rank_map_t *
	  get_rank_map(const GEO& geo) const {
	  return geo.get_property(_pid_rank_map);
	}
	template <class GEO> unsigned long *
	  get_global_idx(const GEO& geo) const {
	  return geo.get_property(_pid_global_idx);
	}

	//@{
	/**
	 * 打印出几何体的 rank_map 供调试时用。
	 */
	template <class GEO> void
	  print_geo_rank_map(const GEO& geo) const {
	  std::cerr << "[" << _forest->rank() << ":" << geo.dimension;
	  rank_map_t * p_map = get_rank_map(geo);
	  if (p_map != NULL) print_rank_map(p_map);
	  std::cerr << "]";
	}

	void print_rank_map(const rank_map_t * p_map) const {
	  typename rank_map_t::const_iterator
	    the_pair = p_map->begin(),
	    end_pair = p_map->end();
	  for (;the_pair != end_pair;++ the_pair) {
	    std::cerr << "(" << the_pair->first << ":"
		      << the_pair->second.old_rank << ":"
		      << the_pair->second.priority << ")";
	  }
	}
	//@}

	/**
	 * 检测几何体是否会在 new_rank 上存在。
	 */
	template <class GEO> bool 
	  is_on_this_new_rank(const GEO& geo,
			      int new_rank) const {
	  rank_map_t * p_map = this->get_rank_map(geo);
	  if (p_map == NULL) return false;
	  typename rank_map_t::iterator
	    the_pair = p_map->find(new_rank);
	  return (the_pair != p_map->end());
	}

	/**
	 * 判断一个几何体在新分区上的存储是否由本进程完成。
	 * 
	 * @param map 几何体的 rank_map 性质指针
	 * @param new_rank 几何体在新分区上的秩
	 * 
	 * @return 在本进程存储为真，否则为假 
	 */
	template <class GEO> bool
	  is_save_on_this_rank(const GEO& geo,
			       int new_rank) const {
	  rank_map_t * p_map = this->get_rank_map(geo);
	  assert (p_map != NULL);
	  typename rank_map_t::iterator
	    the_pair = p_map->find(new_rank);
	  /// 如果几何体的 rank_map 性质中没有 new_rank 的入口
	  if (the_pair == p_map->end()) return false;

	  /// 如果几何体对 new_rank 分区的存储在本进程上完成，返回真
	  return (the_pair->second.old_rank == _forest->rank()); /// 否则返回假
	}

	/** 
	 * 一个几何体的多个备份归并到一个秩上时，完成其相应指针的归并。存储
	 * 的时候，因为一个几何体的多个备份事实上只存储了一次，其它分区上的
	 * 这个指针在载入的时候将会仅仅获得全局标号。因此，如果该几何体事实
	 * 上存储在此归档中，我们会记录下该指针的位置；如果该几何体不存储于
	 * 此归档中，我们直接取出该指针。
	 * 
	 * @param type 类型(2: 存储在此文档；4: 不存储在此文档中)
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

	      /**
	       * 对不在第一个秩上存储的指针进行合并。
	       */
	      std::map<unsigned long, std::list<void *> >::iterator
		the_entry = _global_pointer_to_merge_later.find(global_idx);
	      if (the_entry != _global_pointer_to_merge_later.end()) {
		std::list<void *>::iterator
		  the_ptr = the_entry->second.begin(),
		  end_ptr = the_entry->second.end();
		for (;the_ptr != end_ptr;++ the_ptr) {
		  GEO ** geo_ptr_ptr = (GEO **)(*the_ptr);
		  (*geo_ptr_ptr) = p_geo;
		}
		_global_pointer_to_merge_later.erase(the_entry);
	      }
	    } else {
	      assert (the_pair->second == p_geo);
	    }
	  } else {
	    assert (type == 4);
	    if (the_pair != _global_pointer_to_merge.end()) {
	      p_geo = (GEO *)(the_pair->second);
	    } else {
	      _global_pointer_to_merge_later[global_idx].push_back((void *)(&p_geo));
	    }
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
	    std::pair<int,void *>(p_geo->dimension, (void *)(p_geo));
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
#if 0
	    std::cerr << "sending " << the_pair->second.size()
		      << " items from " << _forest->rank()
		      << " to " << the_pair->first << std::endl;
#endif
	    typename tuple_map_t::iterator
	      the_tuple = the_pair->second.begin(),
	      end_tuple = the_pair->second.end();
	    for (;the_tuple != end_tuple;++ the_tuple) {
	      const unsigned long& global_idx = the_tuple->first;
	      int& dimension = the_tuple->second.first;
	      void *& p_obj = the_tuple->second.second;
	      os << global_idx << dimension << p_obj;
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
#if 0
	    std::cerr << "received " << n_item
		      << " items from " << *the_rank
		      << " to " << _forest->rank() << std::endl;
#endif
	    for (u_int i = 0;i < n_item;++ i) {
	      unsigned long global_idx;
	      int dimension;
	      void * remote_obj;
	      is >> global_idx >> dimension >> remote_obj;
	      std::pair<int,void *>& the_pair = 
		_global_pointer_to_share[*the_rank][global_idx];
	      assert ((dimension == the_pair.first) && 
		      (dimension >= 0) && 
		      (dimension <= 3));

	      if (dimension == 0) {
#define GDIM 0
#define GEO HGeometry<GDIM,dow>
#define OBJ Shared_object<GEO>
		GEO * p_geo = (GEO *)(the_pair.second);
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
		GEO * p_geo = (GEO *)(the_pair.second);
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
		GEO * p_geo = (GEO *)(the_pair.second);
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
		GEO * p_geo = (GEO *)(the_pair.second);
		OBJ * p_info = _forest->get_shared_info(*p_geo);
		if (p_info == NULL) p_info = _forest->new_shared_info(*p_geo);
		p_info->add_clone(*the_rank, (GEO *)remote_obj);
#undef OBJ
#undef GEO
#undef GDIM
	      }
	    }
	  }

#if 0 
	  /**
	   * 检查是否所有发出去信息的几何体都收到了共享信息。
	   */
	  the_pair = _global_pointer_to_share.begin();
	  for (;the_pair != end_pair;++ the_pair) {
	    typename tuple_map_t::iterator
	      the_tuple = the_pair->second.begin(),
	      end_tuple = the_pair->second.end();
	    for (;the_tuple != end_tuple;++ the_tuple) {
	      int& dimension = the_tuple->second.first;
	      void *& p_obj = the_tuple->second.second;
	      if (dimension == 0) {
#define GDIM 0
#define GEO HGeometry<GDIM,dow>
		GEO * p_geo = (GEO *)p_obj;
		assert (_forest->get_shared_info(*p_geo) != NULL);
#undef GEO
#undef GDIM
	      } else if (dimension == 1) {
#define GDIM 1
#define GEO HGeometry<GDIM,dow>
		GEO * p_geo = (GEO *)p_obj;
		assert (_forest->get_shared_info(*p_geo) != NULL);
#undef GEO
#undef GDIM
	      } else if (dimension == 2) {
#define GDIM 2
#define GEO HGeometry<GDIM,dow>
		GEO * p_geo = (GEO *)p_obj;
		assert (_forest->get_shared_info(*p_geo) != NULL);
#undef GEO
#undef GDIM
	      } else if (dimension == 3) {
#define GDIM 3
#define GEO HGeometry<GDIM,dow>
		GEO * p_geo = (GEO *)p_obj;
		assert (_forest->get_shared_info(*p_geo) != NULL);
#undef GEO
#undef GDIM
	      }
	    }
	  }
#endif
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
		      double (LOADER::*value)(const GeometryBM&)) {
	  if (&(ir_mesh.getForest()) != _forest) {
	    std::cerr << "The mesh should be on the same forest of the HLoadBalance." 
		      << std::endl;
	    abort();
	  }

	  RegularMesh<dim,dow>& mesh = ir_mesh.regularMesh();
	  if (! _is_configed) {
	    new_property_id(_pid_loading); 
	    _is_configed = true;
	    new_property_id(_pid_leaf_loading); 
	  }
	  u_int n_ele = mesh.n_geometry(dim);
	  for (u_int i = 0;i < n_ele;++ i) {
	    HGeometry<dim,dow> * p_geo = mesh.template h_geometry<dim>(i);
	    double * p_loading = p_geo->get_property(_pid_loading);
	    double loading = (loader.*value)(mesh.geometry(dim,i));
	    if (p_loading == NULL) {
	      p_loading = p_geo->new_property(_pid_loading);
	      if (p_geo->get_property(_pid_leaf_loading) == NULL) {
		p_geo->new_property(_pid_leaf_loading);
	      }
	      (*p_loading) = loading;
	    } else {
	      (*p_loading) += loading;
	    }
	  }
	}

  private:
	struct _const_loader_t {
	  double _loading;
	  double value(const GeometryBM&) { return _loading; }
	};

  public:
	/** 
	 * 为 ir_mesh 中的每个单元设置上值为 loading 的负载。
	 * 
	 * @param ir_mesh 网格
	 * @param loading 负载值
	 */
	void config(birdview_t& ir_mesh,
		    double loading = 1.0) {
	  _const_loader_t _loader;
	  _loader._loading = loading;
	  this->template config<_const_loader_t>(ir_mesh, _loader, &_const_loader_t::value);
	}

  private:
	struct _var_loader_t {
	  double (*_loading)(const GeometryBM&);
	  double value(const GeometryBM& geo) { return (*_loading)(geo); }
	};

  public:
	/** 
	 * 使用函数 loading 来为 ir_mesh 中的每个单元设置负载。
	 * 
	 * @param ir_mesh 网格
	 * @param loading 设置负载的函数指针
	 */
	void config(birdview_t& ir_mesh,
		    double (*loading)(const GeometryBM&)) {
	  _var_loader_t _loader;
	  _loader._loading = loading;
	  this->template config<_var_loader_t>(ir_mesh, _loader, &_var_loader_t::value);
	}

  private:
	template <class LOADER>
	  struct _const_obj_loader_t {
	    const LOADER * loader;
	    double (LOADER::*fun_ptr)(const GeometryBM&) const;
	    _const_obj_loader_t(const LOADER& _loader) { loader = &_loader; }
	    double value(const GeometryBM& geo) { return (loader->*fun_ptr)(geo); }
	  };

  public:
	template <class LOADER>
	  void config(birdview_t& ir_mesh,
		      const LOADER& loader,
		      double (LOADER::*value)(const GeometryBM&) const) {
	  typedef _const_obj_loader_t<LOADER> loader_t;
	  loader_t _loader(loader);
	  _loader.fun_ptr = value;
	  this->template config<loader_t>(ir_mesh, _loader, &loader_t::value);
	}
	//@}

	//@{
  public:
	/**
	 * 对区域进行重新拆分，获取使得负载平衡的拆分点位置。这个功能由
	 * partition, lump_loading, get_loading 三个函数完成。
	 *
	 * 此处引进了新的参数 tol，其目标是为了使得树结构的根节点可以进行调
	 * 整。如果一个节点的负载大于这个值，则将细分根节点。
	 *
	 * 参数 is_renumerate_root 决定是否对根单元进行HSFC重新排序。
	 *
	 * 参数 f 是一个函数指针，它将三维的坐标进行一个坐标变换，然后使用变
	 * 换后的坐标来完成排序。
	 */
	void partition(u_int n_new_rank = 0,
		       double tol_percent = 0.2,
		       bool is_renumerate_root = false,
		       void (*f)(const double*, double*) = NULL) {
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

	  MPI_Bcast(&total_loading, 1, MPI_DOUBLE, n_rank - 1, comm);
	  double mean_loading = n_new_rank;
	  mean_loading = total_loading/mean_loading;

	  /**
	   * 如果一个根单元的负载比 mean_loading*tol_percent 大，则将根单元细分。
	   */
	  this->refine_root_element(mean_loading*tol_percent, is_renumerate_root, f);

	  _cut_point.clear();
	  _new_rank.clear();

	  u_int idx = 0;
	  _cut_point.push_back(idx); /// put a "0" at beginning

	  loading  = partial_loading - loading;
	  u_int current_rank = (u_int)(floor(loading/mean_loading) + 0.001);
          current_rank = std::min(current_rank, n_new_rank - 1);
	  _new_rank.push_back(current_rank);

	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    idx += 1;
	    double * p_ele_loading = the_ele->get_property(_pid_loading);
	    assert ((p_ele_loading != NULL));
	    loading += (*p_ele_loading);
	    u_int the_rank = (u_int)(floor(loading/mean_loading) + 0.001);
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
	  _is_configed = false;
	  free_property_id(_pid_leaf_loading); 
	}

  private:
	/**
	 * 根节点计算负载。利用get_loading 函数的对后代递归的性质将正则网格的
	 * 单元上的负载集中到背景网格的单元上。
	 */
	double lump_loading() {
	  double loading = 0.0;
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    loading += this->get_loading(*the_ele);
	  }
	  return loading;
	}

  private:

	/**
	 * 计算几何体 geo 上的负载。此函数将计算所有后代的负载，并将其累加到
	 * 自身的负载上。
	 */
	template <class GEO> double
	  get_loading(const GEO& geo) const {
	  if (&geo == NULL) return 0; /// 对空指针返回 0
	  double * p_loading = geo.get_property(_pid_loading);
	  if (p_loading == NULL) {
	    p_loading = geo.new_property(_pid_loading);
	    (*p_loading) = 0;
	  }

	  if (geo.isRefined()) {
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      (*p_loading) += get_loading(*geo.child[i]);
	    }
	    if (geo.get_property(_pid_leaf_loading) != NULL) {
	      for (u_int i = 0;i < geo.n_child;++ i) {
		geo.child[i]->free_property(_pid_loading);
	      }
	    }
	  }
	  return (*p_loading);
	}
	//@}

	//@{
	/**
	 * 对树的根单元进行调整。在数据存储以前，我们对根单元进行加细，使得新
	 * 分区能够更加合理一些。在数据装载以后，我们对根单元进行粗化，使得分
	 * 区上的自适应能够更加灵活一些，特别是在粗化的时候能够充分放粗。
	 */

	/**
	 * 对于负载大于 tol 的根单元进行递归细分，使得新构造的根单元的负载都
	 * 小于 tol。
	 */
	void refine_root_element(double tol, bool is_renum_root, 
				 void (*f)(const double*, double *)) {
	  _nre.clear();

	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    this->refine_root_element(tol, *the_ele, _nre);
	  }
	  _forest->rootElement().swap(_nre);
	  if (is_renum_root) {
	    _forest->renumerateRootElement(f);
	  }
	}

	void refine_root_element(double tol,
				 HGeometry<dim,dow>& geo,
				 typename forest_t::container_t& nre) const {
	  double * p_loading = geo.get_property(_pid_loading);
	  double loading = 0.0;
	  if (p_loading != NULL) loading = *p_loading;
	  if ((loading > tol) &&  /// 负载足够大
	      (geo.isRefined()) && /// 几何体有被细分
	      (geo.child[0]->get_property(_pid_loading) != NULL)) { /// 孩子上也有负载
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      this->refine_root_element(tol, *geo.child[i], nre);
	    }
	  } else {
	    nre.push_back(&geo);
	  }
	}
            
	/**
	 * 对新载入的树的根单元进行粗化。如果一个根单元的父节点不是共享几何
	 * 体，则说明其所有兄弟均在此分区上，则取其父节点为根单元。
	 */
	void coarse_root_element() {
	  property_id_t<bool> pid_is_found;
	  new_property_id(pid_is_found);

	  _nre.clear();
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    this->coarse_root_element(pid_is_found, *the_ele, _nre);
	  }
	  _forest->rootElement().swap(_nre);
	}

	void coarse_root_element(const property_id_t<bool>& pid,
				 HGeometry<dim,dow>& geo,
				 typename forest_t::container_t& nre) const {
	  if (geo.get_property(pid) != NULL) return;
	  if ((geo.parent != NULL) &&
	      (_forest->get_shared_info(*geo.parent) == NULL)) {
	    this->coarse_root_element(pid, *geo.parent, nre);
	  } else {
	    nre.push_back(&geo);
	    geo.new_property(pid);
	  }
	}
	//@}


	//@{
  private:
	/**
	 * 对整个网格的几何体进行标识，判断几何体在负载平衡完成后每个几何体所
	 * 处的分区。这个功能由 set_new_rank_down, set_new_rank_up,
	 * pack_collect_rank, unpack_collect_rank 四个函数完成。
	 */
	void set_new_rank() {
	  new_property_id(_pid_rank_map);

	  set_new_rank_down();
	  set_new_rank_up();

	  if (_forest->n_rank() <= 1) return;

	  /**
	   * 将各个分区上所加的分区标识进行汇总。
	   */
	  Shared_type_filter::only<0> type_filter;
	  MPI_Comm comm = _forest->communicator();
#define SYNC_DATA(D)							\
	  if (dim >= D) {						\
	    sync_data(comm, _forest->template get_shared_list<D>(), *this, \
		      &this_t::template pack_set_new_rank<D>,		\
		      &this_t::template unpack_set_new_rank<D>,		\
		      type_filter);					\
	  }
	  SYNC_DATA(0);
	  SYNC_DATA(1);
	  SYNC_DATA(2);
	  SYNC_DATA(3);
#undef SYNC_DATA

	  //check_new_rank(); /// 调试时用
	}

  public:
	template <int GDIM> void
	  pack_set_new_rank(HGeometry<GDIM,dow> * geo,
			    int remote_rank,
			    Migration::ostream<>& os) {
	  assert ((_forest->get_shared_info(*geo) != NULL));

	  rank_map_t * p_map = get_rank_map(*geo);
	  if (p_map == NULL) {
	    p_map = new_rank_map(*geo);
	  }

	  std::size_t n = p_map->size();
	  os << n;
	  typename rank_map_t::iterator
	    the_pair = p_map->begin(),
	    end_pair = p_map->end();
	  for (;the_pair != end_pair;++ the_pair) {
	    const int& new_rank = the_pair->first;
	    const int& old_rank = the_pair->second.old_rank;
	    int& priority = the_pair->second.priority;
	    if (is_refined_on_new_rank(*geo, new_rank)) {
	      priority = 2;
	    } else if (_forest->is_dummy(*geo)) {
	      priority = 0;
	    } else {
	      priority = 1;
	    }

	    os << new_rank << old_rank << priority;
	  }
	}

	template <int GDIM> void
	  unpack_set_new_rank(HGeometry<GDIM,dow> * geo,
			      int remote_rank,
			      Migration::istream<>& is) {
	  assert ((_forest->get_shared_info(*geo) != NULL));

	  /**
	   * 在每个新分区上，对于每个几何体我们确保只有一个旧分区进行一次的数
	   * 据输出。对于旧分区的选择，遵循以下的原则：
	   *
	   * 1. 如果一个旧分区会输出该几何体，并且不是哑几何体，则具有优先权；
	   * 2. 如果两个旧分区具有完全相同的条件，则秩比较小的那个分区具有优
	   *    先权；
	   */
	  rank_map_t * p_map = get_rank_map(*geo);
	  assert (p_map != NULL);

	  int new_rank;
	  std::size_t n;
	  is >> n;
	  for (u_int i = 0;i < n;++ i) {
	    rank_map_val remote_rank_map;
	    is >> new_rank 
	       >> remote_rank_map.old_rank
	       >> remote_rank_map.priority;
	    typename rank_map_t::iterator 
	      the_pair = p_map->find(new_rank);
	    if (the_pair == p_map->end()) {
	      (*p_map)[new_rank] = remote_rank_map;
	    } else {
	      if (the_pair->second.priority == remote_rank_map.priority) {
		/// 取较小的输出秩
		if (the_pair->second.old_rank > remote_rank_map.old_rank) {
		  the_pair->second.old_rank = remote_rank_map.old_rank;
		}
	      } else if (remote_rank_map.priority > the_pair->second.priority) {
		the_pair->second = remote_rank_map;
	      }
	    }
	  }
	}

  private:
	/**
	 * 对所有共享几何体的存储秩进行验证，确保每个新秩上的每个几何体只被
	 * 存储一次。
	 */
	void check_new_rank() {
	  Shared_type_filter::only<0> type_filter;
	  MPI_Comm comm = _forest->communicator();

#define SYNC_DATA(D)							\
	  if (dim >= D) {						\
	    sync_data(comm, _forest->template get_shared_list<D>(), *this, \
		      &this_t::template pack_check_new_rank<D>,		\
		      &this_t::template unpack_check_new_rank<D>,	\
		      type_filter);					\
	  }
	  SYNC_DATA(0);
	  SYNC_DATA(1);
	  SYNC_DATA(2);
	  SYNC_DATA(3);
#undef SYNC_DATA
	}

  public:
	template <int GDIM> void
	  pack_check_new_rank(HGeometry<GDIM,dow> * geo,
			      int remote_rank,
			      Migration::ostream<>& os) {
	  rank_map_t * p_map = get_rank_map(*geo);
	  assert (p_map != NULL);

	  os << p_map->size();
	  typename rank_map_t::iterator
	    the_pair = p_map->begin(),
	    end_pair = p_map->end();
	  for (;the_pair != end_pair;++ the_pair) {
	    const int& new_rank = the_pair->first;
	    const int& old_rank = the_pair->second.old_rank;
	    int& priority = the_pair->second.priority;
	    os << new_rank << old_rank << priority;
	  }
	}

	template <int GDIM> void
	  unpack_check_new_rank(HGeometry<GDIM,dow> * geo,
				int remote_rank,
				Migration::istream<>& is) {
	  assert ((_forest->get_shared_info(*geo) != NULL));

	  rank_map_t * p_map = get_rank_map(*geo);
	  assert (p_map != NULL);

	  std::size_t n;
	  is >> n;

	  rank_map_t remote_rank_map;
	  for (u_int i = 0;i < n;++ i) {
	    int new_rank;
	    rank_map_val remote_rank_map_val;
	    is >> new_rank 
	       >> remote_rank_map_val.old_rank
	       >> remote_rank_map_val.priority;
	    remote_rank_map[new_rank] = remote_rank_map_val;
	  }

	  typename rank_map_t::iterator
	    the_remote_rank_map = remote_rank_map.begin(),
	    end_remote_rank_map = remote_rank_map.end();
	  for (;the_remote_rank_map != end_remote_rank_map;++ the_remote_rank_map) {
	    const int& remote_new_rank = the_remote_rank_map->first;
	    const int& remote_old_rank = the_remote_rank_map->second.old_rank;
	    const int& remote_priority = the_remote_rank_map->second.priority;

	    typename rank_map_t::iterator 
	      the_pair = p_map->find(remote_new_rank);
	    int errcode = 1;
	    if (the_pair == p_map->end()) {
	      errcode *= 2;
	    } else {
	      if (the_pair->second.old_rank != remote_old_rank) {
		errcode *= 3;
	      }
	      if (the_pair->second.priority != remote_priority) {
		errcode *= 5;
	      }
	    }
	    if (errcode != 1) {
	      std::cerr << " [errcode:" << errcode << "]";
	      print_geo_rank_map(*geo);
	      std::cerr << "-[" << remote_rank 
			<< ":" << GDIM;
	      print_rank_map(&remote_rank_map);
	      std::cerr << "] ";
	      MPI_Abort(_forest->communicator(), errcode);
	    }
	  }
	}

  private:
	/**
	 * 对几何体标记为将位于新分区的第 new_rank 个进程上，并且由本进程进行
	 * 存储。此函数对其所有的顶点、边界和后代做递归。
	 */
	template <class GEO> void
	  geometry_set_new_rank_down(const GEO& geo, int new_rank) const {
	  rank_map_t * p_map = get_rank_map(geo);
	  if (p_map == NULL) {
	    p_map = new_rank_map(geo);
	  } else if (p_map->find(new_rank) != p_map->end()) {
	    return; /// 已经设置过的情形
	  }

	  (*p_map)[new_rank].old_rank = _forest->rank(); /// 将存储秩设为当前的秩
    
	  /// 然后对顶点、边界和后代做递归    
	  for (u_int i = 0;i < geo.n_vertex;++ i) {
	    geometry_set_new_rank_down(*geo.vertex[i], new_rank);
	  }
	  for (u_int i = 0;i < geo.n_boundary;++ i) {
	    geometry_set_new_rank_down(*geo.boundary[i], new_rank);
	  }
	  if (geo.isRefined()) {
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      geometry_set_new_rank_down(*geo.child[i], new_rank);
	    }
	  }
	}

	/**
	 * 根据根单元所属的分区确定根单元所有后代的新分区。
	 */
	void set_new_rank_down() {
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement();
	  u_int n_part = _new_rank.size();
	  for (u_int i = 0;i < n_part;++ i) {
	    for (u_int j = _cut_point[i];j < _cut_point[i + 1];++ j, ++ the_ele) {
	      geometry_set_new_rank_down(*the_ele, _new_rank[i]);
	    }
	  }
	  assert (the_ele == _forest->endRootElement());
	}

	/**
	 * 对几何体标记为将位于新分区的第 new_rank 个进程上。当一个几何体位于
	 * 某个进程上时，其所有的低维组件和祖先当然也在此几何体上，并且我们使
	 * 得其所有的兄弟也在此几何体上，因此此函数对其所有的顶点、边界、祖先
	 * 和兄弟做递归。
	 */
	template <class GEO> void
	  geometry_set_new_rank_up(const GEO& geo, int new_rank) const {
	  rank_map_t * p_map = get_rank_map(geo);
	  if (p_map == NULL) p_map = new_rank_map(geo);

	  (*p_map)[new_rank].old_rank = _forest->rank(); /// 将存储秩设为当前的秩

	  /// 然后对顶点、边界、祖先和兄弟做递归    
	  for (u_int i = 0;i < geo.n_vertex;++ i) {
	    geometry_set_new_rank_up(*geo.vertex[i], new_rank);
	  }
	  for (u_int i = 0;i < geo.n_boundary;++ i) {
	    geometry_set_new_rank_up(*geo.boundary[i], new_rank);
	  }
	  GEO* const& p_parent = geo.parent;
	  if (p_parent != NULL) {
	    if (! is_on_this_new_rank(*p_parent, new_rank)) {
	      geometry_set_new_rank_up(*p_parent, new_rank);
	    }

	    for (u_int i = 0;i < p_parent->n_child;++ i) {
	      GEO* const& p_sibling = p_parent->child[i];
	      if (p_sibling == &geo) continue; /// 是自己
	      if (! is_on_this_new_rank(*p_sibling, new_rank)) {
		geometry_set_new_rank_up(*p_sibling, new_rank);
	      }
	    }
	  }
	}

	/**
	 * 对于根单元以上层次的几何体设定新分区。
	 */
	void set_new_rank_up() {
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement();
	  u_int n_part = _new_rank.size();
	  for (u_int i = 0;i < n_part;++ i) {
	    for (u_int j = _cut_point[i];j < _cut_point[i + 1];++ j, ++ the_ele) {
	      geometry_set_new_rank_up(*the_ele, _new_rank[i]);
	    }
	  }
	  assert (the_ele == _forest->endRootElement());
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
	    HGeometry<dim,dow> * p_geo = &(*the_ele);
	    do {
	      if (p_geo->parent != NULL) {
		p_geo = p_geo->parent;
	      } else break;
	    } while (true);
	    geometry_global_index(*p_geo, idx);
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
	    HGeometry<dim,dow> * p_geo = &(*the_ele);
	    do {
	      if (p_geo->parent != NULL) {
		p_geo = p_geo->parent;
	      } else break;
	    } while (true);
	    geometry_global_index(*p_geo, idx);
	  }
    
	  if (_forest->n_rank() <= 1) return;
	  /**
	   * 将各个分区上具有多个拷贝的几何体的变号进行同步，我们取那个最小
	   * 的标号。
	   */
	  Shared_type_filter::only<0> type_filter;
#define SYNC_DATA(D)							\
	  if (dim >= D) {						\
	    sync_data(comm, _forest->template get_shared_list<D>(), *this, \
		      &this_t::template pack_global_index<D>,		\
		      &this_t::template unpack_global_index<D>,		\
		      type_filter);					\
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
	  assert ((_forest->get_shared_info(*geo) != NULL));

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

	  rank_map_t * p_map = get_rank_map(geo);
	  if (_forest->get_shared_info(geo) != NULL || /// 原本是共享的几何体
	      (p_map != NULL && p_map->size() > 1)) { /// 或者新秩个数大于 1 的情形都编号
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
	void write_config_file(const std::string& dirname, int tag) {
	  char filename[1024];
	  sprintf(filename, "%s/.config", dirname.c_str());
	  std::ofstream os(filename);
	  os << tag << "\t# data file tag\n"
	     << _forest->n_rank() << "\t# number of old rank" 
	     << std::endl;
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

	std::vector<u_int> _n_orphans;
	std::vector<std::list<void *> > _orphans_ptr;

	template <class GEO> void
	  geometry_collect_orphans(GEO& geo, 
				   const property_id_t<>& pid,
				   int new_rank) {
	  if (geo.get_property(pid) != NULL) return;
	  geo.new_property(pid);

	  if ((geo.get_property(_pid_is_saved) == NULL) &&
	      is_save_on_this_rank(geo, new_rank)) {
	    _n_orphans[GEO::dim] += 1;
	    _orphans_ptr[GEO::dim].push_back((void *)(&geo));
	  }

	  for (u_int i = 0;i < geo.n_vertex;++ i) {
	    geometry_collect_orphans(*geo.vertex[i], pid, new_rank);
	  }
	  for (u_int i = 0;i < geo.n_boundary;++ i) {
	    geometry_collect_orphans(*geo.boundary[i], pid, new_rank);
	  }
	  if (geo.isRefined()) {
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      geometry_collect_orphans(*geo.child[i], pid, new_rank);
	    }
	  }
	}    

	void collect_orphans(int new_rank) {
	  _n_orphans.clear();
	  _n_orphans.resize(dim + 1, 0);

	  _orphans_ptr.clear();
	  _orphans_ptr.resize(dim + 1);

	  property_id_t<> pid_is_collected;
	  new_property_id(pid_is_collected);

	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    HGeometry<dim,dow> * p_geo = &(*the_ele);
	    do {
	      if (p_geo->parent != NULL) {
		p_geo = p_geo->parent;
	      } else break;
	    } while (true);
	    geometry_collect_orphans(*p_geo, pid_is_collected, new_rank);
	  }
	}

	template <int D> void 
	  save_dim_orphans(HGeometry_oarchive<this_t>& oa) {
	  typedef HGeometry<D,dow> GEO;
	  u_int n_orphans = _n_orphans[D];
#if 1
	  std::cerr << D << "D("  << n_orphans << ") ";
#endif
	  assert (n_orphans == _orphans_ptr[D].size());
	  oa & boost::serialization::make_binary_object(&n_orphans, sizeof(u_int));
	  typename std::list<void *>::iterator
	    the_ptr = _orphans_ptr[D].begin(), 
	    end_ptr = _orphans_ptr[D].end();
	  for (;the_ptr != end_ptr;++ the_ptr) {
	    GEO * p_geo = (GEO *)(*the_ptr);
	    oa & p_geo;
	  }
	}

	void save_orphans(HGeometry_oarchive<this_t>& oa) {
#if 1
	  std::cerr << "saving orphans: ";
#endif
	  collect_orphans(oa.new_rank());
	  if (dim >= 0) this->template save_dim_orphans<0>(oa);
	  if (dim >= 1) this->template save_dim_orphans<1>(oa);
	  if (dim >= 2) this->template save_dim_orphans<2>(oa);
	  if (dim >= 3) this->template save_dim_orphans<3>(oa);
#if 1
	  std::cerr << std::endl;
#endif
	}

  public:
	void save_data(const std::string& dirname, 
		       birdview_set_t& bvs) {
	  set_new_rank();
	  global_index();

	  int dummy;
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
	  Migration::sync_data_buffer(*_forest);

	  /// 产生一个随机整数作为数据标识
	  int tag = get_comm_int(_forest->communicator());

	  typedef HGeometry_oarchive<this_t> archive_t;
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement();
	  u_int n_part = _new_rank.size();
	  for (u_int i = 0;i < n_part;++ i) {
	    sprintf(command, "mkdir -p %s && mkdir -p %s/%d", 
		    dirname.c_str(), dirname.c_str(), _new_rank[i]);
	    dummy = system(command);
	    sprintf(filename, "%s/%d/%d.dat", dirname.c_str(), 
		    _new_rank[i], _forest->rank());
	    std::ofstream os(filename, std::ios::binary);
	    os.write((char *)&tag, sizeof(tag)); /// 先写出标识

	    {
	      archive_t oa(*this, _new_rank[i], os);
	      new_property_id(_pid_is_saved);

	      u_int n_ele = _cut_point[i + 1] - _cut_point[i];
	      oa & boost::serialization::make_binary_object(&n_ele, sizeof(u_int));
	      std::cerr << "saving data in " << filename 
			<< ", # macro element = " << n_ele << std::endl;

	      for (u_int j = _cut_point[i];j < _cut_point[i + 1];++ j) {
		HGeometry<dim,dow> * p_ele = *(the_ele ++);
		oa & p_ele;
	      }

	      save_orphans(oa);
	      free_property_id(_pid_is_saved);
	    }

	    os.close();
	  }

	  if (_forest->rank() == 0) write_config_file(dirname, tag);
	  _forest->rootElement().swap(_nre); /// 将树的根单元恢复回来
	  _nre.clear();
	  MPI_Barrier(_forest->communicator()); /// 最后做一次同步
	}
	//@}

	//@{
  public:
	int load_config_file(const std::string& dirname, int& tag) {
	  Migration::load_config(dirname);

	  char filename[1024];
	  sprintf(filename, "%s/.config", dirname.c_str());
	  filtering_istream is;
	  Migration::ensured_open_filtering_stream(filename, is);

	  int old_rank;
	  is >> tag >> old_rank;
	  return old_rank;
	}

  private:
	void reconstruct_birdview(HElement<dim,dow>& ele,
				  Migration::data_id_t data_id,
				  u_int idx) {
	  typedef Migration::buffer_t buffer_t;
	  buffer_t& hbuf = ele.h_element->buffer;
	  typename buffer_t::iterator it = hbuf.find(data_id);
	  bool is_active = true;
	  ele.value = 0;
	  if (it == hbuf.end()) {
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
	    assert (ele.h_element->isRefined());
	    ele.value = 1;
	    ele.refine();
	    for (int i = 0;i < ele.n_child;++ i) {
	      reconstruct_birdview(*ele.child[i], data_id, idx);
	    }
	  }
	}

	void reconstruct_birdview(birdview_t& ir_mesh,
				  Migration::data_id_t data_id,
				  u_int idx) {
	  typename birdview_t::RootIterator
	    the_ele = ir_mesh.beginRootElement(),
	    end_ele = ir_mesh.endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    reconstruct_birdview(*the_ele, data_id, idx);
	  }
	}

	/**
	 * 如果几何体是共享几何体，则将几何体自身设置成哑几何体，并对其边界和
	 * 后代做递归。
	 */
	template <class GEO>
	  void set_dummy_flag(GEO& geo) const {
	  if ((_forest->get_shared_info(geo) != NULL) &&
	      (! _forest->is_dummy(geo))) {
	    _forest->dummy(geo);
	  }

	  for (u_int i = 0;i < geo.n_boundary;++ i) {
	    this->set_dummy_flag(*geo.boundary[i]);
	  }
	  if (geo.isRefined()) {
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      this->set_dummy_flag(*geo.child[i]);
	    }
	  }
	}

	/**
	 * 清除几何体身上的哑旗标，并对边界和后代做递归。
	 */
	template <class GEO>
	  void clear_dummy_flag(GEO& geo) const {
	  if (_forest->is_dummy(geo)) {
	    _forest->undummy(geo);
	  }

	  for (u_int i = 0;i < geo.n_boundary;++ i) {
	    this->clear_dummy_flag(*geo.boundary[i]);
	  }
	  if (geo.isRefined()) {
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      this->clear_dummy_flag(*geo.child[i]);
	    }
	  }
	}

	/**
	 * 通过两个步骤来将部分的共享几何体上设置上哑旗标：
	 *
	 * 1. 对所有共享几何体都设置上哑旗标；
	 *
	 * 2. 对当前根单元及其后代都清除哑旗标；
	 *
	 * 对一个单元几何体进行一个操作的时候，对其所有边界和后代都递归进行。
	 * 
	 */
	void set_dummy_flag() {
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    const HGeometry<dim,dow> * p_geo = &(*the_ele);
	    do {
	      if (p_geo->parent != NULL) {
		p_geo = p_geo->parent;
	      } else break;
	    } while (true);
	    set_dummy_flag(*p_geo);
	  }

	  the_ele = _forest->beginRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    clear_dummy_flag(*the_ele);
	  }
	}

  private:
	template <int D>
	  void load_dim_orphans(HGeometry_iarchive<this_t>& ia) const {
	  typedef HGeometry<D,dow> GEO;
	  u_int n_orphans;
	  ia & boost::serialization::make_binary_object(&n_orphans, sizeof(u_int));
#if 1
	  std::cerr << D << "D("  << n_orphans << "), ";
#endif
	  for (u_int i = 0;i < n_orphans;++ i) {
	    GEO * p_geo;
	    ia & p_geo;
	  }
	}

	void load_orphans(HGeometry_iarchive<this_t>& ia) const {
#if 1
	  std::cerr << "loading orphans: ";
#endif
	  if (dim >= 0) this->template load_dim_orphans<0>(ia);
	  if (dim >= 1) this->template load_dim_orphans<1>(ia);
	  if (dim >= 2) this->template load_dim_orphans<2>(ia);
	  if (dim >= 3) this->template load_dim_orphans<3>(ia);
#if 1
	  std::cerr << std::endl;
#endif
	}

	/// 调试用代码：检查树的完整性
	template <class GEO> void
	  geometry_check_integrity(const GEO& geo) const {
	  /**
	   * 如果父亲不是空指针，检查自己确实是父亲的孩子。
	   */
	  if (geo.parent != NULL) {
	    u_int i = 0;
	    for (;i < geo.parent->n_child;++ i) {
	      if (&geo == geo.parent->child[i]) break;
	    }
	    assert ((i < geo.parent->n_child));
	  }

	  /**
	   * 对边界几何体做递归。
	   */
	  for (u_int i = 0;i < geo.n_boundary;++ i) {
	    assert ((geo.boundary[i] != NULL));
	    geometry_check_integrity(*geo.boundary[i]);
	  }
	  if (geo.isRefined()) {
	    /**
	     * 如果自己加密了，确认所有孩子都不是空指针，且所有孩子的父亲指
	     * 针都是自己。
	     */
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      assert ((geo.child[i] != NULL));
	      assert ((&geo == geo.child[i]->parent));
	      geometry_check_integrity(*geo.child[i]);
	    }
	  } else {
	    /**
	     * 如果自己没有加密，确认所有的孩子指针都是空指针。
	     */
	    for (u_int i = 0;i < geo.n_child;++ i) {
	      assert (geo.child[i] == NULL);
	    }
	  }
	}

  public:
	/// 调试用代码：检查树的完整性
	void check_integrity() const {
	  typename forest_t::RootIterator
	    the_ele = _forest->beginRootElement(),
	    end_ele = _forest->endRootElement();
	  for (;the_ele != end_ele;++ the_ele) {
	    HGeometry<dim,dow> * p_geo = &(*the_ele);
	    do {
	      if (p_geo->parent != NULL) {
		p_geo = p_geo->parent;
	      } else break;
	    } while (true);
	    geometry_check_integrity(*p_geo);
	  }
	}

  public:
	void load_data(const std::string& dirname,
		       birdview_set_t& bvs,
		       bool is_has_orphans = false) {
	  clear();

	  int tag = -1;
	  int n_old_rank = load_config_file(dirname, tag);

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
		if (! _global_pointer_to_merge_later.empty()) {
		  std::cerr << "Rank " << _forest->rank() << ": ";
		  std::map<unsigned long, std::list<void *> >::iterator
		    the_entry = _global_pointer_to_merge_later.begin(),
		    end_entry = _global_pointer_to_merge_later.end();
		  for (;the_entry != end_entry;++ the_entry) {
		    std::cerr << "(" << the_entry->first 
			      << "->" << the_entry->second.size()
			      << ")";
		  }
		  std::cerr << std::endl;
		  assert (false);
		}

		this->share_global_pointer();
		this->coarse_root_element();
		this->set_dummy_flag();
		_forest->set_is_used_up();

		//_forest->verify_shared_object();
		//this->check_integrity(); /// 调试时候使用

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

	      /// check the existence of the file
	      if (access(filename, F_OK) != 0) continue;
	      is.open(filename, std::ios::in | std::ios::binary); /// 则打开文件为流
	      int itag;
	      is.read((char *)&itag, sizeof(itag));
	      if (itag != tag) { /// 检查数据文件的标签是否正确
		std::cerr << "*************************************************\n"
			  << "*   DATA TAG NOT MATCH! DATA MAYBE CORRUPTED!!  *\n"
			  << "*************************************************\n"
			  << "Rank: " << _forest->rank() 
			  << ", File: " << filename
			  << ", Tag read: " << itag
			  << ", Tag config: " << tag
			  << std::endl;
		MPI_Abort(_forest->communicator(), 1);
	      }
	      break; /// 打开正确，则开始读入
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

	      if (is_has_orphans) load_orphans(ia);
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
    hlb.load_data(dirname, bvs, is_has_orphans);
  }

  template <class FOREST>
    void load_mesh(const std::string& dirname,
                   FOREST& forest,
                   BirdView<FOREST>& mesh,
                   bool is_has_orphans = false) {
    HLoadBalance<FOREST> hlb(forest);
    BirdViewSet<FOREST> bvs;
    bvs.add(mesh);
    hlb.load_data(dirname, bvs, is_has_orphans);
  }

  template <class FOREST>
    void load_mesh(const std::string& dirname,
                   bool is_has_orphans,
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

    hlb.load_data(dirname, bvs, is_has_orphans);
  }

  template <class FOREST>
    void load_mesh_set(const std::string& dirname,
                       FOREST& forest,
                       BirdViewSet<FOREST>& bvs, 
                       bool is_has_orphans = false) {
    HLoadBalance<FOREST> hlb(forest);
    hlb.load_data(dirname, bvs, is_has_orphans);
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

#endif // __MPI_ULoadBalance_h__

/**
 * end of file
 * 
 */
