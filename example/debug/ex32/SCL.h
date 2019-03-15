/**
 * @file   SCL.h
 * @author Ruo Li <rli@aztec>
 * @date   Fri Oct 23 23:08:42 2009
 * 
 * @brief  
 * 
 * 
 */

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <cstdlib>

#include <functional>
#include <algorithm>
#include <numeric>

#include <AFEPack/Miscellaneous.h>
#include <AFEPack/Operator.h>
#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/mpi/MPI_LoadBalance.h>
#include <AFEPack/mpi/MPI_FaceData.h>
#include <AFEPack/mpi/MPI_Controller.h>

template <int DIM>
struct GeometryCache 
{
  int n_q_pnt;
  AFEPack::Point<DIM> bc;
  double vol;
  std::vector<double> Jxw;
  std::vector<AFEPack::Point<DIM> > q_pnt;
};

template <typename value_type, int DIM>
  struct EdgeCache : public GeometryCache<DIM>
{
  std::vector<std::vector<value_type> > bas_val[2];
  std::vector<std::vector<double> > un;
  std::vector<value_type> u_h_val[2];
};

template <typename value_type, int DIM>
  struct ElementCache : public GeometryCache<DIM>
{
  std::vector<std::vector<value_type> > bas_val;
  std::vector<value_type> u_h_val;
};

#define DIM 2
#define GET_EDGE_CACHE(dgele)                                           \
  EdgeCache<double,DIM>& dgele_cache = edge_cache[(dgele).index()];     \
  double& vol = dgele_cache.vol;                                        \
  std::vector<double>& Jxw = dgele_cache.Jxw;                           \
  int& n_q_pnt = dgele_cache.n_q_pnt;                                   \
  std::vector<AFEPack::Point<DIM> >& q_pnt = dgele_cache.q_pnt;                  \
  std::vector<std::vector<double> >& bas_val0 = dgele_cache.bas_val[0]; \
  std::vector<std::vector<double> >& bas_val1 = dgele_cache.bas_val[1]; \
  std::vector<std::vector<double> >& un = dgele_cache.un;               \
  std::vector<double>& u_h_val0 = dgele_cache.u_h_val[0];               \
  std::vector<double>& u_h_val1 = dgele_cache.u_h_val[1];
#define GET_ELEMENT_CACHE(ele)                                         \
  ElementCache<double,DIM>& ele_cache = element_cache[(ele).index()];  \
  double& vol = ele_cache.vol;                                         \
  std::vector<double>& Jxw = ele_cache.Jxw;                            \
  int& n_q_pnt = ele_cache.n_q_pnt;                                    \
  std::vector<AFEPack::Point<DIM> >& q_pnt = ele_cache.q_pnt;                   \
  std::vector<std::vector<double> >& bas_val = ele_cache.bas_val;    \
  std::vector<double>& u_h_val = ele_cache.u_h_val;

class TimeFunction : public Function<double>
{
 private:
  double t;
  double (*f)(const double *, const double&);
 public:
 TimeFunction(const double& _t, 
              double (*_f)(const double *, 
                           const double&)) : 
  t(_t), f(_f) {}
  virtual ~TimeFunction() {}
 public:
  virtual double value(const double * p) const 
  { return (*f)(p, t); }
};

class SCL
{
 public:
  typedef AFEPack::MPI::HGeometryForest<DIM> forest_t;
  typedef AFEPack::MPI::BirdView<forest_t> ir_mesh_t;
  typedef RegularMesh<DIM> mesh_t;
  typedef DGFEMSpace<double,DIM> fe_space_t;  
  typedef Element<double,DIM> element_t;
  typedef DGElement<double,DIM> dg_element_t;
  typedef FEMFunction<double,DIM> fe_func_t;
  typedef AFEPack::MPI::FaceData::FESpace<std::vector<double>,fe_space_t> data_packer_t;
  typedef data_packer_t::packer_t packer_t;
  typedef AFEPack::MPI::FaceData::Syncer<forest_t,packer_t> syncer_t;
  typedef Vector<double> vector_t;

  std::string meshfile;
  double t;
  double end_t;
  double (*_u_0_)(const double *, const double&);
  std::vector<double> (*_f_)(const AFEPack::Point<DIM>&, const double&);
  double (*_v_)(const AFEPack::Point<DIM>&, const double&);
  void (*_u_b_)(const AFEPack::Point<DIM>&, double&, const std::vector<double>&, const double&, const int&);


 private:
  TemplateGeometry<DIM>                      triangle_template_geometry;
  CoordTransform<DIM>                        triangle_coord_transform;
  UnitOutNormal<DIM>                         triangle_unit_out_normal;
  TemplateDOF<DIM>                           triangle_template_dof;
  BasisFunctionAdmin<double,DIM>             triangle_basis_function;

  TemplateGeometry<DIM>                      twin_triangle_template_geometry;
  CoordTransform<DIM>                        twin_triangle_coord_transform;
  UnitOutNormal<DIM>                         twin_triangle_unit_out_normal;
  TemplateDOF<DIM>                           twin_triangle_template_dof;
  BasisFunctionAdmin<double,DIM>             twin_triangle_basis_function;

  TemplateGeometry<DIM-1>                    interval_template_geometry;
  CoordTransform<DIM-1,DIM>                  interval_to2d_coord_transform;
  
  std::vector<TemplateElement<double,DIM> >  template_element;
  std::vector<TemplateDGElement<DIM-1,DIM> > dg_template_element;
  
  forest_t                                   htree;
  ir_mesh_t *                                ir_mesh;
  fe_space_t *                               fem_space;
  fe_func_t *                                u_h;
  double                                     dt;

  Indicator<DIM>                             indicator;

  std::vector<EdgeCache<double,DIM> >        edge_cache; // 边上的数据缓冲
  std::vector<ElementCache<double,DIM> >     element_cache; /// 单元上数据缓冲

  /**
   * 需要两个新的类成员变量：
   *   1. 数据交换器，使其内部可以存储交换映射图；
   *   2. 一个旗标表征数据交换映射图是否需要更新；
   *
   * 在每次进行负载平衡或者是网格自适应后，有限元空间被重新建立，此时
   * 就需要更新数据交换映射图。
   */
  syncer_t                                   syncer; /// 数据交换器
  bool                                       is_update_transmit_map;

 public:
  void run();
  void save_data();

 private:
  void initialize();
  void time_step_length();
  void step_forward();
  void build_fe_space();
  void init_edge_cache();
  void init_element_cache();
  void update_edge_cache(fe_func_t&);

  void adapt_mesh();
  void get_indicator();

  void load_balance();

  void buildMeshInfo(int, int);
  void get_rhs(fe_func_t&, vector_t&);
  void bound_value(const AFEPack::Point<DIM>&, 
                   double&, 
                   const std::vector<double>&, 
                   const dg_element_t&); 
  void bound_value(const std::vector<AFEPack::Point<DIM> >&, 
                   std::vector<double>&, 
                   const std::vector<std::vector<double> >&, 
                   const dg_element_t&);
  double flux(const AFEPack::Point<DIM>&, double, double,
              const std::vector<double>& n) const;
};

/**
 * end of file
 * 
 */


