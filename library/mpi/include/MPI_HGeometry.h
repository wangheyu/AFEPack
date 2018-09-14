 /**
  * @file   MPI_HGeometry.h
  * @author Ruo Li <rli@aztec>
  * @date   Fri Sep 11 16:23:39 2009
  * 
  * @brief  
  * 
  * 
  */

#ifndef __MPI_HGeometry_h__
#define __MPI_HGeometry_h__

#include <set>
#include "../../include/HGeometry.h"
#include "MPI.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {

  template <int DIM, int DOW, class MATCHER> class HGeometryForest;
  template <class FOREST> class BirdView;
  template <class FOREST> class BirdViewSet;
  template <class FOREST> class HGeometryMatcher;
  template <class FOREST> class HLoadBalance;

  template <int DOW>
    struct PointDistance {
      int value(const Point<DOW>& p0,
                const Point<DOW>& p1,
                double tol) const {
        return (distance(p0, p1) < tol)?0:-1;
      }
    };

  /**
   * HGeometryForest 作为 HGeometryTree 的派生类，并不是一系列的
   * HGeometryTree 的组合，而只是对 HGeometryTree 作了进一步的包装，使
   * 得其可以完成不同分区之间的数据交换。
   *
   * 其重新实现的 readMesh 函数可以实现对具有分区类型的网格的读入。
   */
  template <int DIM, int DOW=DIM, class MATCHER=PointDistance<DOW> >
    class HGeometryForest : public HGeometryTree<DIM,DOW> {
  private:
  typedef HGeometryForest<DIM,DOW,MATCHER> this_t;
  public:
  enum {dim = DIM, dow = DOW};

  typedef HGeometryTree<DIM,DOW> tree_t;
  typedef MATCHER matcher_t;

  template <int D>
  const Shared_ptr_list<HGeometry<D,DOW> >&
  get_shared_list() const {
    const void * p_list;
    switch (D) {
    case 0: p_list = &_shared_list_0d; break;
    case 1: p_list = &_shared_list_1d; break;
    case 2: p_list = &_shared_list_2d; break;
    case 3: p_list = &_shared_list_3d; break;
    }
    return *((const Shared_ptr_list<HGeometry<D,DOW> > *)p_list);
  }

  /**
   * 一个几何体的多份拷贝中，其中部分拷贝所在的进程的秩是所有拷贝中最
   * 小的，这些拷贝所在的秩我们称为`首秩'(primary rank)。此函数返回一
   * 个几何体是否是在首秩上。GEO 需要为 HGeometry 类型。
   */
  template <int GDIM>
  bool is_on_primary_rank(const HGeometry<GDIM,DOW>& geo) const {
    typedef HGeometry<GDIM,DOW> geo_t;
    typedef Shared_object<geo_t> obj_t;
    obj_t * p_obj = this->get_shared_info(geo);
    if (p_obj == NULL) {
      return true;
    } else {
      return p_obj->is_on_primary_rank(_rank);
    }
  }

  /**
   * 返回一个几何体的首秩。
   */
  template <int GDIM>
  int primary_rank(const HGeometry<GDIM,DOW>& geo) const {
    typedef HGeometry<GDIM,DOW> geo_t;
    typedef Shared_object<geo_t> obj_t;
    obj_t * p_obj = this->get_shared_info(geo);
    if (p_obj == NULL) {
      return _rank;
    } else {
      return p_obj->primary_rank(_rank);
    }
  }

  /**
   * 判断一个几何体是否是主对象。
   */
  template <int GDIM>
  bool is_primary_geometry(const HGeometry<GDIM,DOW>& geo) const {
    typedef HGeometry<GDIM,DOW> geo_t;
    typedef Shared_object<geo_t> obj_t;
    obj_t * p_obj = this->get_shared_info(geo);
    if (p_obj == NULL) {
      return true;
    } else {
      return p_obj->is_primary_object(_rank);
    }
  }

  /**
   * 作为来源于本几何遗传树的正则网格 mesh，返回其 dim 维的第 idx 个几
   * 何体是否是共享的几何体。
   */
  bool is_geometry_shared(const RegularMesh<DIM,DOW>& mesh,
                          int dim, int idx) const {
    switch(dim) {
    case 0: return (get_shared_info(*mesh.template h_geometry<0>(idx)) != NULL);
    case 1: return (get_shared_info(*mesh.template h_geometry<1>(idx)) != NULL);
    case 2: return (get_shared_info(*mesh.template h_geometry<2>(idx)) != NULL);
    case 3: return (get_shared_info(*mesh.template h_geometry<3>(idx)) != NULL);
    }
    return false;
  }

  const matcher_t& matcher() const { return _matcher; }
  matcher_t& matcher() { return _matcher; }

  private:
  MPI_Comm _comm;
  int _rank;
  int _n_rank;
  matcher_t _matcher;

  /**
   * 表示几何体的匹配信息是否已经发送的性质。几何体对所有应该进行匹配
   * 的进程都已经发送了自身的匹配信息之后，此性质将会被设置。函数
   * is_shared_info_sent 和 set_shared_info_sent 对此性质进行设置和查
   * 询。
   */
  property_id_t<> _pid_shared_info_sent; 
  template <class GEO> bool
  is_shared_info_sent(GEO& geo) const {
    return (geo.get_property(_pid_shared_info_sent) != NULL);
  }
  template <class GEO> void
  set_shared_info_sent(GEO& geo) const {
    if (! is_shared_info_sent(geo)) {
      geo.new_property(_pid_shared_info_sent);
    }
  }
  template <class GEO> void
  clear_shared_info_sent(GEO& geo) const {
    geo.free_property(_pid_shared_info_sent);
  }

  /**
   * 在使用不是基于最初背景单元的负载平衡的时候(MPI_ULoadbalance.h)，有
   * 些几何体会在分区上存储为哑几何体。在这样的情况下，每个分区上存储的
   * 几何体中，有部分将永远不可能出现在网格中。所谓哑几何体，我们定义其
   * 为将不可能会在分区的网格中出现的几何体，包括根单元的前辈几何体和根
   * 单元的那些不是根单元兄弟，以及这些哑单元几何体的部分边界和顶点。
   */
  property_id_t<> _pid_dummy;
  template <class GEO> bool
  is_dummy(GEO& geo) const {
    return (geo.get_property(_pid_dummy) != NULL);
  }
  template <class GEO> void
  dummy(GEO& geo) const {
    geo.new_property(_pid_dummy);
  }
  template <class GEO> void
  undummy(GEO& geo) const {
    geo.free_property(_pid_dummy);
  }

  /**
   * 几何体的共享信息性质。几何体的共享信息直接附着在几何体上，而在
   * HGeometryForest 上存储这个性质的指针的列表，便于对所有共享几何体
   * 同时进行数据同步。
   *
   * 根据几何体的维数不同，相应的这一组信息的变量名也不相同，下面按维
   * 数进行定义。
   */
#define GDIM 0
#define GEO HGeometry<GDIM,DOW>
#define OBJ Shared_object<GEO>
  private:
  property_id_t<OBJ> _pid_so_0d;
  Shared_ptr_list<GEO> _shared_list_0d;

  public:
  OBJ * get_shared_info(const GEO& geo) const {
    return geo.get_property(_pid_so_0d);
  }
  OBJ * new_shared_info(GEO& geo) {
    OBJ * p_obj = geo.new_property(_pid_so_0d);
    p_obj->local_pointer() = &geo;
    _shared_list_0d.push_back(p_obj);
    return p_obj;
  }
  void erase_shared_info(GEO& geo) {
    OBJ * p_obj = geo.get_property(_pid_so_0d);
    if (p_obj == NULL) return;
    _shared_list_0d.erase(std::find(_shared_list_0d.begin_ptr(), 
                                    _shared_list_0d.end_ptr(), 
                                    p_obj));
  }
#undef OBJ
#undef GEO
#undef GDIM

#define GDIM 1
#define GEO HGeometry<GDIM,DOW>
#define OBJ Shared_object<GEO>
  private:
  property_id_t<OBJ> _pid_so_1d;
  Shared_ptr_list<GEO> _shared_list_1d;

  public:
  OBJ * get_shared_info(const GEO& geo) const {
    return geo.get_property(_pid_so_1d);
  }
  OBJ * new_shared_info(GEO& geo) {
    OBJ * p_obj = geo.new_property(_pid_so_1d);
    p_obj->local_pointer() = &geo;
    _shared_list_1d.push_back(p_obj);
    return p_obj;
  }
  void erase_shared_info(GEO& geo) {
    OBJ * p_obj = geo.get_property(_pid_so_1d);
    if (p_obj == NULL) return;
    _shared_list_1d.erase(std::find(_shared_list_1d.begin_ptr(), 
                                    _shared_list_1d.end_ptr(), 
                                    p_obj));
  }
#undef OBJ
#undef GEO
#undef GDIM

#define GDIM 2
#define GEO HGeometry<GDIM,DOW>
#define OBJ Shared_object<GEO>
  private:
  property_id_t<OBJ> _pid_so_2d;
  Shared_ptr_list<GEO> _shared_list_2d;

  public:
  OBJ * get_shared_info(const GEO& geo) const {
    return geo.get_property(_pid_so_2d);
  }
  OBJ * new_shared_info(GEO& geo) {
    OBJ * p_obj = geo.new_property(_pid_so_2d);
    p_obj->local_pointer() = &geo;
    _shared_list_2d.push_back(p_obj);
    return p_obj;
  }
  void erase_shared_info(GEO& geo) {
    OBJ * p_obj = geo.get_property(_pid_so_2d);
    if (p_obj == NULL) return;
    _shared_list_2d.erase(std::find(_shared_list_2d.begin_ptr(), 
                                    _shared_list_2d.end_ptr(), 
                                    p_obj));
  }
#undef OBJ
#undef GEO
#undef GDIM

#define GDIM 3
#define GEO HGeometry<GDIM,DOW>
#define OBJ Shared_object<GEO>
  private:
  property_id_t<OBJ> _pid_so_3d;
  Shared_ptr_list<GEO> _shared_list_3d;

  public:
  OBJ * get_shared_info(const GEO& geo) const {
    return geo.get_property(_pid_so_3d);
  }
  OBJ * new_shared_info(GEO& geo) {
    OBJ * p_obj = geo.new_property(_pid_so_3d);
    p_obj->local_pointer() = &geo;
    _shared_list_3d.push_back(p_obj);
    return p_obj;
  }
  void erase_shared_info(GEO& geo) {
    OBJ * p_obj = geo.get_property(_pid_so_3d);
    if (p_obj == NULL) return;
    _shared_list_3d.erase(std::find(_shared_list_3d.begin_ptr(), 
                                    _shared_list_3d.end_ptr(), 
                                    p_obj));
  }
#undef OBJ
#undef GEO
#undef GDIM

  bool lock() { return HGeometryTree<DIM,DOW>::lock(); }
  void unlock() { HGeometryTree<DIM,DOW>::unlock(); }

  public:
  HGeometryForest() { new_property(); }
  HGeometryForest(MPI_Comm comm) {
    set_communicator(comm);
    new_property();
  }
  virtual ~HGeometryForest() { this->clear(); }

  void clear() {
    _shared_list_0d.clear();
    _shared_list_1d.clear();
    _shared_list_2d.clear();
    _shared_list_3d.clear();

    free_property();

    /**
     * 考虑到在使用 ULoadBalance 进行负载平衡以后，根单元可能不是最粗
     * 的背景单元了，直接调用基类的 clear 函数可能造成内存泄漏，我们此
     * 处先将最粗的背景单元列表进行重建，然后再调用基类的 clear 函数对
     * 内存进行清空的操作。
     */
    typename tree_t::container_t nre;
    property_id_t<bool> pid_is_found;
    new_property_id(pid_is_found);
    typename tree_t::RootIterator
    the_ele = this->beginRootElement(),
    end_ele = this->endRootElement();
    for (;the_ele != end_ele;++ the_ele) {
      this->collect_parent_element(pid_is_found, *the_ele, nre);
    }
    free_property_id(pid_is_found);

    /// 使用重建的根单元
    this->rootElement().swap(nre);
    nre.clear();
    tree_t::clear(); /// 然后进行内存清除

    new_property();
  }

  private:
  void collect_parent_element(const property_id_t<bool>& pid,
                              HGeometry<dim,dow>& geo,
                              typename tree_t::container_t& nre) const {
    if (geo.get_property(pid) != NULL) return;
    if (geo.parent != NULL) {
      this->collect_parent_element(pid, *geo.parent, nre);
    } else {
      nre.push_back(&geo);
      geo.new_property(pid);
    }
  }

  public:
  void new_property() {
    new_property_id(_pid_shared_info_sent);
    new_property_id(_pid_dummy);
    if (DIM >= 0) new_property_id(_pid_so_0d);
    if (DIM >= 1) new_property_id(_pid_so_1d);
    if (DIM >= 2) new_property_id(_pid_so_2d);
    if (DIM >= 3) new_property_id(_pid_so_3d);
  }

  void free_property() {
    free_property_id(_pid_shared_info_sent);
    free_property_id(_pid_dummy);
    if (DIM >= 0) free_property_id(_pid_so_0d);
    if (DIM >= 1) free_property_id(_pid_so_1d);
    if (DIM >= 2) free_property_id(_pid_so_2d);
    if (DIM >= 3) free_property_id(_pid_so_3d);
  }

  void set_communicator(MPI_Comm comm) {
    _comm = comm;
    MPI_Comm_rank(_comm, &_rank);
    MPI_Comm_size(_comm, &_n_rank);
  }
  MPI_Comm communicator() const { return _comm; }
  int rank() const { return _rank; }
  int n_rank() const { return _n_rank; }
  void readMesh(const std::string&);

  void eraseRootElement(u_int level = 1);
  void renumerateRootElement(void (*f)(const double *, double *) = NULL);

  /**
   * 从根单元开始，向祖先方向设置几何体被使用的旗标。
   */
  void set_is_used_up();
  template <class HGEO> void
  set_is_used_up(HGEO& geo) const;

  /**
   * 将所有的共享几何体的共享信息发送旗标(_pid_shared_info_sent)设为已
   * 经发送过了。此函数在数据从存储上读入以后使用，使得旧有的共享信息
   * 不必重新在几何体匹配的时候发送。
   */
  void set_shared_info_sent() const {
    if (dim >= 3) set_dim_shared_info_sent<3>();
    if (dim >= 2) set_dim_shared_info_sent<2>();
    if (dim >= 1) set_dim_shared_info_sent<1>();
    if (dim >= 0) set_dim_shared_info_sent<0>();
  }
  template <int D>
  void set_dim_shared_info_sent() const {
    typedef HGeometry<D,dow> geo_t;
    typedef Shared_object<geo_t> obj_t;
    typedef Shared_ptr_list<geo_t> list_t;

    const list_t& geo_ptr_list = this->template get_shared_list<D>();
    typename list_t::const_iterator
    the_geo_ptr = geo_ptr_list.begin(),
    end_geo_ptr = geo_ptr_list.end();
    for (;the_geo_ptr != end_geo_ptr;++ the_geo_ptr) {
      geo_t * p_geo = the_geo_ptr->local_pointer();
      this->set_shared_info_sent(*p_geo);
    }
  }

  /**
   * 调试用代码：对分区之间的共享几何体的信息进行验证，确证一个几何体的
   * 共享列表中的其它拷贝也将自己作为共享几何体。
   */
  public:
  template <int GDIM> void
  pack_verify_shared_object(HGeometry<GDIM,dow> * geo,
                            int remote_rank,
                            Migration::ostream<>& os) {
    typedef HGeometry<GDIM,dow> geo_t;
    os << geo; /// 发送出去本身的指针

    Shared_object<geo_t> * p_info = this->get_shared_info(*geo);
    std::size_t n = p_info->size();
    os << n;
    typename std::multimap<int,Remote_pointer<geo_t> >::iterator
    the_entry = p_info->begin(), end_entry = p_info->end();
    for (;the_entry != end_entry;++ the_entry) {
      os << the_entry->first
         << the_entry->second.type
         << the_entry->second.ptr;
    }
  }

  template <int GDIM> void
  unpack_verify_shared_object(HGeometry<GDIM,dow> * geo,
                              int remote_rank,
                              Migration::istream<>& is) {
    typedef HGeometry<GDIM,dow> geo_t;
    geo_t * p_remote_geo;
    is >> p_remote_geo; /// 读入远程几何体的指针

    Shared_object<geo_t> * p_info = this->get_shared_info(*geo);
    assert (p_info != NULL);

    std::size_t n;
    is >> n;
    Remote_pointer<geo_t> p_remote_ptr;
    for (std::size_t i = 0;i < n;++ i) {
      int orank;
      is >> orank >> p_remote_ptr.type >> p_remote_ptr.ptr;
      if ((orank == this->rank()) &&
          (p_remote_ptr.ptr == geo)) {
        orank = remote_rank;
        p_remote_ptr.ptr = p_remote_geo;
      }
      if (! p_info->is_duplicate_entry(orank, p_remote_ptr)) {
        MPI_Abort(this->communicator(), 0);
      }
    }
  }

  void verify_shared_object() {
    if (this->n_rank() <= 1) return;
    Shared_type_filter::all type_filter;
#define SYNC_DATA(D)                                                    \
    if (dim >= D) {                                                     \
      sync_data(this->communicator(),                                   \
                this->template get_shared_list<D>(), *this,             \
                &this_t::template pack_verify_shared_object<D>,         \
                &this_t::template unpack_verify_shared_object<D>,       \
                type_filter);                                           \
    }
    SYNC_DATA(0);
    SYNC_DATA(1);
    SYNC_DATA(2);
    SYNC_DATA(3);
#undef SYNC_DATA
  }

  private:
  void eraseRootElementOneLevel();
  template <class GEO>
  void nullParent(GEO& geo) {
    geo.parent = NULL;
    for (u_int i = 0;i < GEO::n_boundary;++ i) {
      this->nullParent(*geo.boundary[i]);
    }
  }
  template <class GEO>
  void tryDeleteGeometry(GEO * p_geo, 
                         const property_id_t<bool>& pid) {
    for (u_int i = 0;i < GEO::n_boundary;++ i) {
      if (p_geo->boundary[i]->get_property(pid) != NULL) {
        p_geo->boundary[i] = NULL;
      } else {
        this->tryDeleteGeometry(p_geo->boundary[i], pid);
      }
    }
    bool * p_prp = p_geo->get_property(pid);
    assert (p_prp == NULL);
    p_geo->new_property(pid);
    this->erase_shared_info(*p_geo); /// 顺手将其共享信息列表中入口删除
  }
  template <class GEO>
  void deleteGeometry(GEO * p_geo) {
    for (u_int i = 0;i < GEO::n_boundary;++ i) {
      if (p_geo->boundary[i] != NULL) {
        this->deleteGeometry(p_geo->boundary[i]);
      }
    }
    delete p_geo;
  }

  friend class BirdView<this_t>;
  friend class HGeometryMatcher<this_t>;
  friend class HLoadBalance<this_t>;

  public:
  typedef BirdView<this_t> birdview_t;
  typedef BirdViewSet<this_t> birdview_set_t;
  typedef HLoadBalance<this_t> load_balancer_t;

  }; // HGeometryForest

  /**
   * BirdView 是 IrregularMesh 的派生类，其中 semiregularize 函数被重
   * 新实现以完成数据交换，使得不同分区之间的网格满足半正则化的要求。
   */
  template <class FOREST>
    class BirdView : public IrregularMesh<FOREST::dim,FOREST::dow> {
  public:
    enum {dim = FOREST::dim, dow = FOREST::dow};
    typedef IrregularMesh<dim, dow> base_t;
    typedef FOREST forest_t;
  public:
  BirdView() : base_t() {}
    explicit BirdView(forest_t& forest) : base_t(forest) {}
    virtual ~BirdView() {}

    forest_t& getForest() const {
      return dynamic_cast<forest_t&>(this->geometryTree());
    }
    virtual void semiregularize();
    void eraseRootElement(u_int level = 1);
  private:
    void eraseRootElementOneLevel();
  };

  /**
   * 一些 BirdView 的集合，在进行数据存储和载入中使用。数据迁移将对此
   * 集合中的网格进行操作。
   */
  template <class FOREST>
    class BirdViewSet : public std::list<BirdView<FOREST>*> {
  private:
    typedef BirdView<FOREST> birdview_t;
    typedef std::list<birdview_t*> base_t;
  public:
    typedef _Deref_iterator<typename base_t::iterator, birdview_t> iterator;
    typedef _Deref_iterator<typename base_t::const_iterator, const birdview_t> const_iterator;
  public:
    void add(birdview_t& bv) { base_t::push_back(&bv); }
    void erase(birdview_t& bv) { base_t::remove(&bv); }

    iterator begin() { return base_t::begin(); }
    iterator end() { return base_t::end(); }

    const_iterator begin() const { return base_t::begin(); }
    const_iterator end() const { return base_t::end(); }

    typename base_t::iterator begin_ptr() { return base_t::begin(); }
    typename base_t::iterator end_ptr() { return base_t::end(); }

    typename base_t::const_iterator begin_ptr() const { return base_t::begin(); }
    typename base_t::const_iterator end_ptr() const { return base_t::end(); }
  };

  template <class FOREST>
    BirdViewSet<FOREST> 
    make_set(BirdView<FOREST>& mesh0, 
             BirdView<FOREST>& mesh1) {
    BirdViewSet<FOREST> bvs;
    bvs.add(mesh0);
    bvs.add(mesh1);
    return bvs;
  }
  template <class FOREST>
    BirdViewSet<FOREST> 
    make_set(BirdView<FOREST>& mesh0, 
             BirdView<FOREST>& mesh1,
             BirdView<FOREST>& mesh2) {
    BirdViewSet<FOREST> bvs;
    bvs.add(mesh0);
    bvs.add(mesh1);
    bvs.add(mesh2);
    return bvs;
  }
  template <class FOREST>
    BirdViewSet<FOREST> 
    make_set(BirdView<FOREST>& mesh0, 
             BirdView<FOREST>& mesh1,
             BirdView<FOREST>& mesh2,
             BirdView<FOREST>& mesh3) {
    BirdViewSet<FOREST> bvs;
    bvs.add(mesh0);
    bvs.add(mesh1);
    bvs.add(mesh2);
    bvs.add(mesh3);
    return bvs;
  }

  /**
   * HGeometryMatcher 完成不同分区上最主要的数据匹配和交换的工作，这包
   * 括不同分区上的共享几何体的匹配，以及不同分区上几何体上附着的信息
   * 的交换。完成不同分区上共享几何体的匹配的过程中遵循下面的原则：
   *
   * 1. 每个几何体由谁来匹配：
   *
   *    a. 如果一个几何体有父亲几何体，则由其父亲完成匹配，因为有父亲的
   *       几何体的共享结构必定与其父亲相同；
   *
   *    b. 如果一个几何体没有父亲几何体，则由其高一维的几何体完成匹配，
   *       因为对没有父亲的几何体，其共享结构必定与其高一维几何体相同；
   *
   *    c. 对于点几何体，由其所属的线段几何体对其进行匹配，对于自适应
   *       过程中产生的点几何体，其共享结构必定与其所属的线段几何体相
   *       同；
   *
   * 2. 几何体的匹配信息可能发送多次，但是都是完全一致的信息；
   *
   * 下面实现的另一个交换是对几何体上的是否处于 USED 状态进行特定的同步，
   * 使得其中的一个拷贝的成为 USED 状态时，其他的拷贝都会成为 USED 状态。
   * 
   */
  template <class FOREST>
    class HGeometryMatcher {
  public:
    enum {dim = FOREST::dim, dow = FOREST::dow};
    typedef FOREST forest_t;
  private:
    typedef HGeometryMatcher<forest_t> this_t;
    typedef typename forest_t::matcher_t matcher_t;

    HTools _tools;
    forest_t * _forest;
    bool _is_operated;

    const matcher_t& matcher() const { return _forest->matcher(); }

    template <class GEO> bool 
      is_shared_info_sent(GEO& geo) const {
      return _forest->is_shared_info_sent(geo);
    }
    template <class GEO> void 
      set_shared_info_sent(GEO& geo) const {
      _forest->set_shared_info_sent(geo);
    }
    template <class GEO> Shared_object<GEO> * 
      new_shared_info(GEO& geo) const {
      return _forest->new_shared_info(geo);
    }
    template <class GEO> Shared_object<GEO> *
      get_shared_info(GEO& geo) const {
      return _forest->get_shared_info(geo);
    }

  public:
  HGeometryMatcher(forest_t& forest) : _forest(&forest) {}

    /**
     * 对几何遗传树中的几何体进行匹配，使得几何体的复制品在不同的分区
     * 中具有完全相同的树结构。在树结构被修改了以后需要使用此操作进行
     * 匹配。函数 pack_match_geometry 和 unpack_match_geometry 作为对
     * 此操作中的函数指针的支持。返回是否确实进行了操作。
     */
    bool match_geometry();
    template <int GDIM> bool 
      is_pack_match_geometry(HGeometry<GDIM,dow> * geo) const;
    template <int GDIM> void 
      pack_match_geometry(HGeometry<GDIM,dow> * geo,
                          int remote_rank,
                          Migration::ostream<>& os);
    template <int GDIM> void 
      unpack_match_geometry(HGeometry<GDIM,dow> * geo,
                            int remote_rank,
                            Migration::istream<>& is);

  public:
    /**
     * 设置点匹配器。
     */
    void set_matcher(const matcher_t& _matcher) { matcher = _matcher; }

    /**
     * 对分区间共享的几何体的 index 进行同步化，如果有一个备份的设为了
     * USED 状态，就将所有备份都设置为 USED 状态。函数
     * pack_sync_is_used 和 unpack_sync_is_used 作为对 此操作中的函数指
     * 针的支持。
     */
    bool sync_is_used();
    template <int GDIM> void 
      pack_sync_is_used(HGeometry<GDIM,dow> * geo,
                        int remote_rank,
                        Migration::ostream<>& os);
    template <int GDIM> void 
      unpack_sync_is_used(HGeometry<GDIM,dow> * geo,
                          int remote_rank,
                          Migration::istream<>& is);

  private:
    /**
     * 计算一个几何体做匹配时候的参考坐标点。
     */
    template <class GEO>
      void geometry_reference_point(GEO& geo, 
                                    Point<dow>& pnt) const {
      int n_vtx = geo.n_vertex;
      for (int n = 0;n < dow;++ n) {
        pnt[n] = 0.0;
        for (int i = 0;i < n_vtx;++ i) {
          pnt[n] += (*geo.vertex[i])[n];
        }
        pnt[n] /= n_vtx;
      }
    }
    void geometry_reference_point(HGeometry<0,dow>& geo, 
                                  Point<dow>& pnt) const {
      for (int n = 0;n < dow;++ n) {
        pnt[n] = geo[n];
      }
    }
  };


#include "MPI_HGeometry.templates.h"

} // namespace MPI

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_HGeometry_h__

/**
 * end of file
 * 
 */

