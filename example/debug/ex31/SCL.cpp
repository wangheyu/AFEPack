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

void SCL::initialize()
{
  htree.set_communicator(MPI_COMM_WORLD);
  htree.readMesh(meshfile);
  ir_mesh.reinit(htree);
  ir_mesh.globalRefine(1);
  ir_mesh.semiregularize();
  ir_mesh.regularize(false);
  
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
  init_element_cache();
  init_edge_cache();

  u_h.reinit(fem_space);
  TimeFunction u(t, _u_0_);
  Operator::L2Project(u, u_h, Operator::LOCAL_LEAST_SQUARE, 3);
}

void SCL::build_fe_space() {
  mesh_t& mesh = ir_mesh.regularMesh();
  fem_space.reinit(mesh, template_element, dg_template_element);
  fem_space.element().resize(mesh.n_geometry(DIM));
  for (u_int i = 0;i < mesh.n_geometry(DIM);i ++) {
    if (mesh.geometry(DIM,i).n_vertex() == 3) {
      fem_space.element(i).reinit(fem_space, i, 0);
    } else {
      fem_space.element(i).reinit(fem_space, i, 1);
    }
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  fem_space.dgElement().resize(mesh.n_geometry(DIM-1));
  for (u_int i = 0;i < mesh.n_geometry(DIM-1);i ++)
    fem_space.dgElement(i).reinit(fem_space, i, 0);
  fem_space.buildDGElement();
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
  edge_cache.resize(fem_space.n_DGElement());

  fe_space_t::DGElementIterator 
    the_dgele = fem_space.beginDGElement(),
    end_dgele = fem_space.endDGElement();
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
  element_cache.resize(fem_space.n_element());

  fe_space_t::ElementIterator 
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
  mesh_t& mesh = ir_mesh.regularMesh();
  double rank_dt = 1.0;
  fe_space_t::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    double h = sqrt(2*element_cache[the_ele->index()].vol);
    const std::vector<int>& ele_dof = the_ele->dof();
    double v = (*_v_)(AFEPack::Point<DIM>(), u_h(ele_dof[0]));
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
  fe_func_t u_h1(u_h);

  // 3 Step TVD Runge-Kutta Scheme
  get_rhs(u_h1, rhs);
  u_h1.add(dt, rhs);
  t = old_t + dt;

  get_rhs(u_h1, rhs);
  u_h1.sadd(0.25, 0.75, u_h, dt/4, rhs);
  t = old_t + dt/2;

  get_rhs(u_h1, rhs);
  u_h.sadd(1./3, 2./3, u_h1, 2*dt/3, rhs);
  t = old_t + dt;
}

void SCL::update_edge_cache(fe_func_t& v_h) 
{
  /// 为数据交换准备空间
  syncer_t syncer(htree, data_packer_t::get_packer(fem_space));

  fe_space_t::DGElementIterator 
    the_dgele = fem_space.beginDGElement(),
    end_dgele = fem_space.endDGElement();
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
      syncer.attach_send_buffer(*the_dgele, &u_h_val0);
      syncer.attach_recv_buffer(*the_dgele, &u_h_val1);
    } else { /// 物理区域的边界
      u_h_val1 = u_h_val0;
      bound_value(q_pnt, u_h_val1, un, *the_dgele);
    }
  }

  syncer.sync(); /// 完成分区间的数据交换
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
  rhs.reinit(fem_space.n_dof());
  fe_space_t::DGElementIterator 
    the_dgele = fem_space.beginDGElement(),
    end_dgele = fem_space.endDGElement();
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
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
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
  u_h.writeOpenDXData(filename);
  if (htree.rank() == 0) {
    std::cout << "Data saved in u_h*.h." << std::endl;
  }
}

/**
 * end of file
 * 
 */


