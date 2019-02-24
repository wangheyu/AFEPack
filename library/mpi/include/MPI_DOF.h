/**
 * @file   MPI_DOF.h
 * @author Ruo Li <rli@aztec>
 * @date   Mon Oct 26 16:56:21 2009
 * 
 * @brief  进行全局自由度排序和分区间的自由度匹配。
 * 
 * 
 */

#ifndef __MPI_DOF_h__
#define __MPI_DOF_h__

#include "MPI_HGeometry.h"
#include "MPI_SyncProp.h"

AFEPACK_OPEN_NAMESPACE

namespace MPI {
  namespace DOF {

    template <class FOREST, class FESPACE>
      class GlobalIndex : private std::vector<int> {
    public:
    typedef FOREST forest_t;
    typedef FESPACE fe_space_t;
    private:
    typedef std::vector<int> base_t;
    typedef GlobalIndex<forest_t,fe_space_t> this_t;

    const forest_t * _forest;
    const fe_space_t * _fe_sp;

    u_int _n_global_dof;
    u_int _n_primary_dof;
    int _first_idx;

    property_id_t<int> _pid_min_rank;
    property_id_t<BinaryBuffer<> > _pid_global_idx;

    enum { UNUSED_IDX = -1 };

    public:
    GlobalIndex() : _forest(NULL), _fe_sp(NULL) {}
    explicit GlobalIndex(const forest_t& forest) :
    _forest(&forest), _fe_sp(NULL) {}
    GlobalIndex(const forest_t& forest, const fe_space_t& fe_sp) :
    _forest(&forest), _fe_sp(&fe_sp) {}
    GlobalIndex(const this_t& gi) :
    base_t(gi), _forest(gi._forest), _fe_sp(gi._fe_sp) {}
    
    this_t& operator=(const this_t& gi) {
      base_t::operator=(gi);
      _forest = gi._forest;
      _fe_sp = gi._fe_sp;
    }

    //@{
    /**
     * 取森林。
     */
    const forest_t& forest() const { return *_forest; }
    void set_forest(const forest_t& forest) { _forest = &forest; }
    //@}
    
    //@{
    /**
     * 取有限元空间。
     */
    const fe_space_t& fe_space() const { return *_fe_sp; }
    void set_fe_space(const fe_space_t& fe_sp) { _fe_sp = &fe_sp; }
    //@}

    //@{
    /**
     * 取全局自由度个数，局部的在首秩上的自由度个数，局部的自由度个数，
     * 以及第一个在首秩上的自由度的全局指标。
     * 
     */
    u_int n_global_dof() const { return _n_global_dof; } 
    u_int n_primary_dof() const { return _n_primary_dof; }
    u_int n_local_dof() const { return base_t::size(); }
    int first_index() const { return _first_idx; }
    int last_index() const { return _first_idx + _n_primary_dof - 1; }
    //@}

    /**
     * 建立全局自由度指标。
     */
    void build();
    /**
     * 从本地自由度中抽出所有在首秩上的自由度，放入数组 idx_ptr 中。
     */
    void build_primary_index(int * idx_ptr) const;

    //@{
    /**
     * 建立 Trilinos 的 Epetra 的映射。
     */
    template <class MAP> void build_epetra_map(MAP& map) const;
    template <class MAP, class INVMAP> void 
      build_epetra_map(MAP& map, INVMAP& inv_map) const;
    //@}

    //@{
    /**
     * 将 Hypre 的全局向量分发为局部向量。
     *
     * Hypre 的全局向量仅仅能够读出在本进程上的元素，而分布式的有限元
     * 函数需要部分远程的元素才能够构建完整的信息。本函数将 Hypre 的向
     * 量中的元素进行一些分发，使得分布式的有限元函数能够获得其需要的
     * 远程元素的值，从而其能够被 AFEPack 的内部使用。
     *
     * 模板参数 INVEC 需为 HYPRE_IJVector， OUTVEC 需为 Vector<double>
     * 或者其派生类。
     *
     */
    template <class INVEC, class OUTVEC> void
    scatter_hypre_vector(INVEC& iv, OUTVEC& ov) const;
    //@}

    //@{
    /**
     * 取局部的第 i 个自由度的全局指标。
     */
    const int& operator()(u_int i) const { return base_t::operator[](i); }
    int& operator()(u_int i) { return base_t::operator[](i); }
    //@}

    /**
     * 从本地到全局指标数组的转换。
     */
    template <class LC, class GC> void
      translate(const LC& lc, GC& gc) const;

    //@{
    /**
     * 内部使用。在建立全局指标时打包和解包数据。
     */
    template <int GDIM> void
      pack_global_idx(HGeometry<GDIM,fe_space_t::dow> * geo,
                      int remote_rank,
                      Migration::ostream<>& os);
    template <int GDIM> void
      unpack_global_idx(HGeometry<GDIM,fe_space_t::dow> * geo,
                        int remote_rank,
                        Migration::istream<>& is);
    //@}

    /// 查询一个自由度是否在主几何体上
    bool is_dof_on_minimal_geometry(u_int i) const;

    /// 查询一个自由度的首秩
    int get_dof_minimal_rank(u_int i) const;

    private:
    void set_dof_minimal_rank(); /// 设置所有自由度到最小秩

    void sync_idx(); /// 进程间进行自由度指标同步

    std::vector<bool> _idopg; /// 记录每个自由度是否在主几何体上
    };

#include "MPI_DOF.templates.h"
  }
}

AFEPACK_CLOSE_NAMESPACE

#endif // __MPI_DOF_h__

/**
 * end of file
 * 
 */
