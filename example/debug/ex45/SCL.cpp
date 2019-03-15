/**
 * @file   SCL.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Oct 23 23:12:36 2009
 * 
 * @brief  
 * 
 * 
 */

#include "SCL.h"

void SCL::load_balance()
{
  /// 首先判断是否需要做负载平衡
  RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
  int n_ele = mesh.n_geometry(DIM);
  int n_max_ele;
  MPI_Allreduce(&n_ele, &n_max_ele, 1, MPI_INT, MPI_MAX,
                htree.communicator());
  int is_rank_load_balance = (0.9*n_max_ele > n_ele)?1:0; /// 单元个数差太多
  int is_load_balance;
  MPI_Allreduce(&is_rank_load_balance, &is_load_balance, 1, MPI_INT, MPI_SUM,
                htree.communicator());
  if (is_load_balance == 0) return;

  dump(htree.n_rank(), ".hlb");
  load(".hlb");
}

void SCL::dump(int n_rank, const std::string& dir) {
  /// 然后对数据进行重新分区存储
  AFEPack::MPI::HLoadBalance<SCL::forest_t> hlb(htree);
  hlb.config(*ir_mesh); /// 使用 ir_mesh 进行负载计算
  hlb.partition(n_rank); /// 重新分区

  AFEPack::MPI::BirdViewSet<SCL::forest_t> bvs;
  bvs.add(*ir_mesh); /// 存储时对 ir_mesh 也进行存储

  Migration::clear_data_buffer(htree); /// 清除不用的数据
  Migration::data_id_t id = Migration::name_to_id("u_h");
  if (! Migration::is_valid(id)) { /// 获取数据ID，如果没有则新建一个
    id = Migration::register_data_name("u_h");
  }
  Migration::export_fe_func(*u_h, id); /// 将有限元函数放到几何遗传树上

  if (htree.rank() == 0) {
    std::cout << "Dumping the data ... " << std::endl;
    char cmd[8192];
    sprintf(cmd, "rm -rf %s", dir.c_str());
    system(cmd);
  }
  MPI_Barrier(htree.communicator()); /// 同步以后
  hlb.save_data(dir, bvs); /// 存储数据
}

void SCL::load(const std::string& dir) {
  /// 清除内存
  if (u_h != NULL) delete u_h; /// 有限元函数和有限元空间均会重造
  if (fem_space != NULL) delete fem_space; /// 所以这两个变量都直接删除
  if (ir_mesh != NULL) delete ir_mesh; /// 网格和几何遗传树会重新载入
  htree.clear(); /// 因此清除其中的数据空间

  /// 数据载入
  AFEPack::MPI::HLoadBalance<SCL::forest_t> hlb(htree);
  AFEPack::MPI::BirdViewSet<SCL::forest_t> bvs;
  ir_mesh = new ir_mesh_t(htree);
  bvs.add(*ir_mesh); /// 存储时对 ir_mesh 也进行存储

  if (htree.rank() == 0) {
    std::cout << "Loading the data ... " << std::endl;
  }
  hlb.load_data(dir, bvs, true); /// 载入数据
  ir_mesh->semiregularize(); /// 重造正则网格
  ir_mesh->regularize(false);

  build_fe_space(); /// 重建有限元空间
  u_h = new fe_func_t(*fem_space); /// 和有限元函数
  Migration::data_id_t id = Migration::name_to_id("u_h");
  if (! Migration::is_valid(id)) { /// 获取数据ID，如果没有则新建一个
    id = Migration::register_data_name("u_h");
  }
  Migration::import_fe_func(*u_h, id); /// 将有限元函数从几何遗传树上读出
  Migration::clear_data_buffer(htree); /// 将几何遗传树上的数据清除

  update_edge_cache(*u_h); /// 更新边界上数据缓冲
}

void SCL::initialize()
{
  htree.set_communicator(MPI_COMM_WORLD);

  /**< 对数据迁移全局环境进行初始化 */
  Migration::initialize(htree.communicator());
  AFEPack::MPI::load_forest(meshfile, htree);
  //htree.readMesh(meshfile);
  ir_mesh = new ir_mesh_t(htree);
  ir_mesh->semiregularize();
  ir_mesh->regularize(false);
  
  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_coord_transform.readData("triangle.crd_trs");
  triangle_unit_out_normal.readData("triangle.out_nrm");
  triangle_template_dof.reinit(triangle_template_geometry);
  triangle_template_dof.readData("triangle.0.tmp_dof");
  triangle_basis_function.reinit(triangle_template_dof);
  triangle_basis_function.readData("triangle.0.bas_fun");

  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  twin_triangle_coord_transform.readData("twin_triangle.crd_trs");
  twin_triangle_unit_out_normal.readData("twin_triangle.out_nrm");
  twin_triangle_template_dof.reinit(twin_triangle_template_geometry);
  twin_triangle_template_dof.readData("twin_triangle.0.tmp_dof");
  twin_triangle_basis_function.reinit(twin_triangle_template_dof);
  twin_triangle_basis_function.readData("twin_triangle.0.bas_fun");

  template_element.resize(2);
  template_element[0].reinit(triangle_template_geometry,
                             triangle_template_dof,
                             triangle_coord_transform,
                             triangle_basis_function,
                             triangle_unit_out_normal);
  template_element[1].reinit(twin_triangle_template_geometry,
                             twin_triangle_template_dof,
                             twin_triangle_coord_transform,
                             twin_triangle_basis_function,
                             twin_triangle_unit_out_normal);

  interval_template_geometry.readData("interval.tmp_geo");
  interval_to2d_coord_transform.readData("interval.to2d.crd_trs");

  dg_template_element.resize(1);
  dg_template_element[0].reinit(interval_template_geometry,
				interval_to2d_coord_transform);

  build_fe_space();

  u_h = new fe_func_t(*fem_space);
  TimeFunction u(t, _u_0_);
  Operator::L2Project(u, *u_h, Operator::LOCAL_LEAST_SQUARE, 3);

  update_edge_cache(*u_h);

  for (u_int i = 0;i < 5;++ i) {
    adapt_mesh();
    Operator::L2Project(u, *u_h, Operator::LOCAL_LEAST_SQUARE, 3);
  }
}

void SCL::get_indicator()
{
  RegularMesh<DIM>& mesh = ir_mesh->regularMesh();
  indicator.reinit(mesh);

  fe_space_t::DGElementIterator 
    the_dgele = fem_space->beginDGElement(),
    end_dgele = fem_space->endDGElement();
  for (;the_dgele != end_dgele;++ the_dgele) {
    EdgeCache<double, DIM>& ec = edge_cache[the_dgele->index()];
    element_t * nei_ele0 = the_dgele->p_neighbourElement(0);
    element_t * nei_ele1 = the_dgele->p_neighbourElement(1);
    double a = fabs(ec.u_h_val[0][0] - ec.u_h_val[1][0])*ec.vol;
    indicator[nei_ele0->index()] += a;
    if (nei_ele1 != NULL) {
      indicator[nei_ele1->index()] += a;
    }
  }
}

void SCL::adapt_mesh() {
  get_indicator(); /// 计算自适应指示子

  double convergence_order = 1.0;
  double adapt_tolerence = 5.0e-04;
  double ADAPT_RATIO = 0.03;
  double rtol = (1.33333*pow(2.0, DIM + 1.0))*adapt_tolerence;
  double ctol = adapt_tolerence/(1.33333*pow(2.0, DIM + 1.0));

  size_t n_ele, n_ele_to_refine, n_ele_to_coarse;
  n_ele = fem_space->n_element();
  n_ele_to_refine = std::count_if(indicator.begin(), indicator.end(), 
                                  std::bind2nd(std::greater<double>(), rtol));
  n_ele_to_coarse = std::count_if(indicator.begin(), indicator.end(), 
                                  std::bind2nd(std::less<double>(), ctol));
  /// 如果需要加密的网格数太少，则不做自适应
  int is_rank_adapt = (n_ele_to_refine + n_ele_to_coarse < ADAPT_RATIO*n_ele)?0:1;
  int is_adapt; /// 是否做自适应需要在所有进程间进行同步
  MPI_Allreduce(&is_rank_adapt, &is_adapt, 1, MPI_INT, MPI_SUM, 
                htree.communicator());
  if (is_adapt == 0) return;

  ir_mesh_t * old_ir_mesh = ir_mesh;
  fe_space_t * old_fem_space = fem_space;
  fe_func_t * old_u_h = u_h;

  ir_mesh = new ir_mesh_t(*old_ir_mesh);
  MeshAdaptor<DIM> mesh_adaptor(*ir_mesh);
  mesh_adaptor.convergenceOrder() = convergence_order;
  mesh_adaptor.refineStep() = 1;
  mesh_adaptor.setIndicator(indicator);
  mesh_adaptor.tolerence() = adapt_tolerence;
  mesh_adaptor.adapt();

  /**< 正则化自适应后的网格 */
  ir_mesh->semiregularize();
  ir_mesh->regularize(false);

  /**< 建立新的有限元空间和将解更新到新的空间 */
  build_fe_space();
  u_h = new fe_func_t(*fem_space);
  Operator::L2Project(*old_u_h, *u_h, Operator::LOCAL_LEAST_SQUARE, 3);

  /**< 释放旧的网格、有限元空间以及解的内存 */
  delete old_u_h;
  delete old_fem_space;
  delete old_ir_mesh;

  update_edge_cache(*u_h); /// 更新边界缓冲
}

void SCL::build_fe_space() {
  mesh_t& mesh = ir_mesh->regularMesh();
  fem_space = new fe_space_t(mesh, template_element, dg_template_element);
  u_int n_ele = mesh.n_geometry(DIM);
  fem_space->element().resize(n_ele);
  for (u_int i = 0;i < n_ele;i ++) {
    if (mesh.geometry(DIM,i).n_vertex() == 3) {
      fem_space->element(i).reinit(*fem_space, i, 0);
    } else {
      fem_space->element(i).reinit(*fem_space, i, 1);
    }
  }
  fem_space->buildElement();
  fem_space->buildDof();
  fem_space->buildDofBoundaryMark();

  u_int n_dg_ele = mesh.n_geometry(DIM-1);
  fem_space->dgElement().resize(n_dg_ele);
  for (u_int i = 0;i < n_dg_ele;i ++)
    fem_space->dgElement(i).reinit(*fem_space, i, 0);
  fem_space->buildDGElement();

  /**< 初始化缓冲 */
  init_element_cache();
  init_edge_cache();
}

void SCL::run()
{
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  double start_time = MPI_Wtime();
  double last_time = start_time;
  double this_time;
  initialize();
  do {
    step_forward();
    adapt_mesh();
    load_balance(); /// 完成负载平衡
    AFEPack::MPI::getControl();

    if (htree.rank() == 0) {
      last_time = this_time;
      this_time = MPI_Wtime();
      std::cout << "t = " << t 
                << ", dt = " << dt
                << ", wall time = " << this_time - start_time
                << ", step time = " << this_time - last_time
                << std::endl;
    }
  } while(end_t - t > 2.0e-16);
}

void SCL::init_edge_cache()
{
  edge_cache.resize(fem_space->n_DGElement());

  fe_space_t::DGElementIterator 
    the_dgele = fem_space->beginDGElement(),
    end_dgele = fem_space->endDGElement();
  for (;the_dgele != end_dgele;++ the_dgele) {
    GET_EDGE_CACHE(*the_dgele);

    double tmp_vol = the_dgele->templateElement().volume();
    const QuadratureInfo<DIM-1>& qi = the_dgele->findQuadratureInfo(3);
    n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jacobian = the_dgele->local_to_global_jacobian(qi.quadraturePoint());
    q_pnt = the_dgele->local_to_global(qi.quadraturePoint());
    const element_t& neigh0 = the_dgele->neighbourElement(0);
    bas_val0 = neigh0.basis_function_value(q_pnt);
    un = unitOutNormal(q_pnt, neigh0, *the_dgele);
    if (the_dgele->p_neighbourElement(1) != NULL) {
      const element_t& neigh1 = the_dgele->neighbourElement(1);
      bas_val1 = neigh1.basis_function_value(q_pnt);
    }
    vol = 0;
    Jxw.resize(n_q_pnt);
    for (int l = 0;l < n_q_pnt;l ++) {
      Jxw[l] = tmp_vol*jacobian[l]*qi.weight(l);
      vol += Jxw[l];
    }
  }
}

void SCL::init_element_cache()
{
  element_cache.resize(fem_space->n_element());

  fe_space_t::ElementIterator 
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele;++ the_ele) {
    GET_ELEMENT_CACHE(*the_ele);

    double tmp_vol = the_ele->templateElement().volume();
    const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(1);
    n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jacobian = the_ele->local_to_global_jacobian(qi.quadraturePoint());
    q_pnt = the_ele->local_to_global(qi.quadraturePoint());
    bas_val = the_ele->basis_function_value(q_pnt);

    vol = 0;
    Jxw.resize(n_q_pnt);
    for (int l = 0;l < n_q_pnt;l ++) {
      Jxw[l] = tmp_vol*jacobian[l]*qi.weight(l);
      vol += Jxw[l];
    }
  }
}

void SCL::time_step_length()
{
  double CFL = 0.25;
  mesh_t& mesh = ir_mesh->regularMesh();
  double rank_dt = 1.0;
  fe_space_t::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele;++ the_ele) {
    double h = sqrt(2*element_cache[the_ele->index()].vol);
    const std::vector<int>& ele_dof = the_ele->dof();
    double v = (*_v_)(AFEPack::Point<DIM>(), (*u_h)(ele_dof[0]));
    rank_dt = std::min(rank_dt, h/v);
  }
  MPI_Allreduce(&rank_dt, &dt, 1, MPI_DOUBLE, 
                MPI_MIN, htree.communicator());
  dt *= CFL;
  if (t + dt > end_t) dt = end_t - t;
}

void SCL::step_forward()
{
  time_step_length();

  double old_t = t;
  vector_t rhs;
  fe_func_t u_h1(*u_h);

  // 3 Step TVD Runge-Kutta Scheme
  get_rhs(u_h1, rhs);
  u_h1.add(dt, rhs);
  t = old_t + dt;

  get_rhs(u_h1, rhs);
  u_h1.sadd(0.25, 0.75, *u_h, dt/4, rhs);
  t = old_t + dt/2;

  get_rhs(u_h1, rhs);
  u_h->sadd(1./3, 2./3, u_h1, 2*dt/3, rhs);
  t = old_t + dt;
}

void SCL::update_edge_cache(fe_func_t& v_h) 
{
  /// 为数据交换准备空间
  typedef std::vector<double> vec_t;
  property_id_t<vec_t> pid_buf;
  new_property_id(pid_buf);

  fe_space_t::DGElementIterator 
    the_dgele = fem_space->beginDGElement(),
    end_dgele = fem_space->endDGElement();
  for (;the_dgele != end_dgele;++ the_dgele) {
    GET_EDGE_CACHE(*the_dgele);

    const element_t& neigh0 = the_dgele->neighbourElement(0);
    u_h_val0 = v_h.value(bas_val0, neigh0);
    if (the_dgele->p_neighbourElement(1) != NULL) { /// 内部边
      const element_t& neigh1 = the_dgele->neighbourElement(1);
      u_h_val1 = v_h.value(bas_val1, neigh1);
    } else if (the_dgele->boundaryMark() <= 0) { /// 分区间的内边界
      /**
       * 这里的程序是本例的重点，我们取出数据发送缓冲区，然后将要发送
       * 的数据填入其中。
       */
      vec_t * p_buf = fem_space->new_property<vec_t>(*the_dgele,
                                                     pid_buf);
      (*p_buf) = u_h_val0;
    } else { /// 物理区域的边界
      u_h_val1 = u_h_val0;
      bound_value(q_pnt, u_h_val1, un, *the_dgele);
    }
  }

  AFEPack::MPI::PropSyncer<forest_t,vec_t> syncer(htree, pid_buf);
  syncer.sync<DIM-1>(); /// 完成分区间的数据交换

  the_dgele = fem_space->beginDGElement();
  for (;the_dgele != end_dgele;++ the_dgele) {
    GET_EDGE_CACHE(*the_dgele);
    /**
     * 第一个邻居单元为空，并且自身的边界标识为 0 的情形就是分区间的内
     * 边界。首先拿出数据接收缓冲区，然后将其中的数据拷贝到变量中。
     */
    if (the_dgele->p_neighbourElement(1) == NULL &&
        the_dgele->boundaryMark() <= 0) {
      vec_t * p_buf = fem_space->get_property<vec_t>(*the_dgele,
                                                     pid_buf);
      u_h_val1 = (*p_buf);
    }
  }
}

double SCL::flux(const AFEPack::Point<DIM>& p,
                 double u_l, double u_r,
                 const std::vector<double>& n) const 
{
  std::vector<double> f_l((*_f_)(p, u_l));
  std::vector<double> f_r((*_f_)(p, u_r));
  double alpha0 = (*_v_)(p, u_l);
  double alpha1 = (*_v_)(p, u_r);
  double alpha = std::max(alpha0, alpha1);
  return 0.5*((innerProduct(f_l, n) + innerProduct(f_r, n)) -
              alpha*(u_r - u_l));
}

void SCL::get_rhs(fe_func_t& v_h, vector_t& rhs)
{
  update_edge_cache(v_h);

  double alpha = 1.0;
  rhs.reinit(fem_space->n_dof());
  fe_space_t::DGElementIterator 
    the_dgele = fem_space->beginDGElement(),
    end_dgele = fem_space->endDGElement();
  for (;the_dgele != end_dgele;++ the_dgele) {
    GET_EDGE_CACHE(*the_dgele);

    const element_t& neigh0 = the_dgele->neighbourElement(0);
    const std::vector<int>& ele_dof0 = neigh0.dof();
    if (the_dgele->p_neighbourElement(1) != NULL) { // 内部边
      const element_t& neigh1 = the_dgele->neighbourElement(1);
      const std::vector<int>& ele_dof1 = neigh1.dof();
      for (int l = 0;l < n_q_pnt;l ++) {
        double a = Jxw[l]*flux(q_pnt[l], u_h_val0[l], u_h_val1[l], un[l]);
        rhs(ele_dof0[0]) -= a;
        rhs(ele_dof1[0]) += a;
      }
    } else { // 分区之间的边界和区域的物理边界上的边
      for (int l = 0;l < n_q_pnt;l ++) {
        double a = Jxw[l]*flux(q_pnt[l], u_h_val0[l], u_h_val1[l], un[l]);
        rhs(ele_dof0[0]) -= a;
      }
    }
  }

  fe_space_t::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for (;the_ele != end_ele;++ the_ele) {
    const std::vector<int>& ele_dof = the_ele->dof();
    rhs(ele_dof[0]) /= element_cache[the_ele->index()].vol;
  }
}

void SCL::bound_value(const AFEPack::Point<DIM>& q_point,
                      double& v_h, 
                      const std::vector<double>& n,
                      const dg_element_t& dgele)
{
  (*_u_b_)(q_point, v_h, n, t, dgele.boundaryMark());
}

void SCL::bound_value(const std::vector<AFEPack::Point<DIM> >& q_point, 
                      std::vector<double>& v_h, 
                      const std::vector<std::vector<double> >& n,
                      const dg_element_t& dgele)
{
  for (int i = 0;i < q_point.size();i ++)
    (*_u_b_)(q_point[i], v_h[i], n[i], t, dgele.boundaryMark());
}

void SCL::save_data()
{
  char filename[1024];
  sprintf(filename, "u_h%d.dx", htree.rank());
  u_h->writeOpenDXData(filename);
  if (htree.rank() == 0) {
    std::cout << "Data saved in u_h${rank}.h." << std::endl;
  }
}

/**
 * end of file
 * 
 */


