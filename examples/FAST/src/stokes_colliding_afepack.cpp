/**
 * @file stokes_colliding_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：碰撞流算例。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 碰撞流算例：设置相向流动边界条件，观察压力场和速度场在中心区域的耦合。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define main stokes_poiseuille_afepack_hidden_main
#include "stokes_poiseuille_afepack.cpp"
#undef main

namespace
{
  struct CollidingResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       divergence_l2_norm = 0.0;
  };

  double reference_x(const double *p)
  {
    return 2.0 * p[0] - 1.0;
  }

  double reference_y(const double *p)
  {
    return 2.0 * p[1] - 1.0;
  }

  double colliding_velocity_x_exact(const double *p)
  {
    const double x = reference_x(p);
    const double y = reference_y(p);
    return 20.0 * x * y * y * y;
  }

  double colliding_velocity_y_exact(const double *p)
  {
    const double x = reference_x(p);
    const double y = reference_y(p);
    return 5.0 * x * x * x * x - 5.0 * y * y * y * y;
  }

  std::vector<double> colliding_velocity_x_gradient(const double *p)
  {
    const double x = reference_x(p);
    const double y = reference_y(p);
    return {40.0 * y * y * y, 120.0 * x * y * y};
  }

  std::vector<double> colliding_velocity_y_gradient(const double *p)
  {
    const double x = reference_x(p);
    const double y = reference_y(p);
    return {40.0 * x * x * x, -40.0 * y * y * y};
  }

  double colliding_pressure_exact(const double *p)
  {
    const double x = reference_x(p);
    const double y = reference_y(p);
    return 120.0 * x * x * y - 40.0 * y * y * y;
  }

  std::vector<double> colliding_pressure_gradient(const double *p)
  {
    const double x = reference_x(p);
    const double y = reference_y(p);
    return {480.0 * x * y, 240.0 * (x * x - y * y)};
  }

  double colliding_rhs_x(const double *p)
  {
    (void)p;
    return 0.0;
  }

  double colliding_rhs_y(const double *p)
  {
    (void)p;
    return 0.0;
  }

  double rate(const double previous_error,
              const double current_error,
              const double previous_h,
              const double current_h)
  {
    if (previous_error <= 0.0 || current_error <= 0.0 ||
        previous_h <= 0.0 || current_h <= 0.0 ||
        std::abs(previous_h - current_h) < 1.0e-15)
      return 0.0;
    return std::log(previous_error / current_error) /
           std::log(previous_h / current_h);
  }

  double pressure_l2_error_mod_constant(
    const FEMFunction<double, dim> &pressure)
  {
    const auto &space = pressure.femSpace();
    double area = 0.0;
    double error_integral = 0.0;
    double error_square_integral = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(6);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto values = pressure.value(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double point[dim] = {q_point[q][0], q_point[q][1]};
            const double error =
              values[q] - colliding_pressure_exact(point);
            area += jxw;
            error_integral += jxw * error;
            error_square_integral += jxw * error * error;
          }
      }

    const double best_constant = error_integral / area;
    return std::sqrt(std::max(0.0,
                              error_square_integral -
                                area * best_constant * best_constant));
  }

  CollidingResult solve_colliding_stokes(const std::string &mesh_file,
                                         const std::string &output_prefix,
                                         const double h)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> velocity_geometry;
    CoordTransform<dim, dim> velocity_transform;
    TemplateDOF<dim> velocity_dof;
    BasisFunctionAdmin<double, dim, dim> velocity_basis(velocity_dof);
    auto velocity_template = make_triangle_template(2,
                                                    velocity_geometry,
                                                    velocity_transform,
                                                    velocity_dof,
                                                    velocity_basis);

    TemplateGeometry<dim> pressure_geometry;
    CoordTransform<dim, dim> pressure_transform;
    TemplateDOF<dim> pressure_dof;
    BasisFunctionAdmin<double, dim, dim> pressure_basis(pressure_dof);
    auto pressure_template = make_triangle_template(1,
                                                    pressure_geometry,
                                                    pressure_transform,
                                                    pressure_dof,
                                                    pressure_basis);

    auto velocity_space = build_space(mesh, velocity_template);
    auto pressure_space = build_space(mesh, pressure_template);

    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    std::vector<double> matrix(static_cast<std::size_t>(n_total) * n_total,
                               0.0);
    std::vector<double> rhs(n_total, 0.0);

    for (int e = 0; e < velocity_space.n_element(); ++e)
      {
        const auto &velocity_element = velocity_space.element(e);
        const auto &pressure_element = pressure_space.element(e);
        const double volume = velocity_element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          velocity_element.findQuadratureInfo(6);
        const auto jacobian =
          velocity_element.local_to_global_jacobian(
            quad_info.quadraturePoint());
        const auto q_point =
          velocity_element.local_to_global(quad_info.quadraturePoint());
        const auto velocity_values =
          velocity_element.basis_function_value(q_point);
        const auto velocity_grads =
          velocity_element.basis_function_gradient(q_point);
        const auto pressure_values =
          pressure_element.basis_function_value(q_point);

        const auto &velocity_element_dofs = velocity_element.dof();
        const auto &pressure_element_dofs = pressure_element.dof();
        const int n_local_velocity_dofs = velocity_element_dofs.size();
        const int n_local_pressure_dofs = pressure_element_dofs.size();

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double p[dim] = {q_point[q][0], q_point[q][1]};
            const double f[dim] = {colliding_rhs_x(p), colliding_rhs_y(p)};

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    const double stiffness =
                      viscosity * jxw *
                      innerProduct(velocity_grads[i][q],
                                   velocity_grads[j][q]);
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_total,
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                stiffness);
                  }

                for (int component = 0; component < dim; ++component)
                  rhs[velocity_index(component, vi, n_velocity_dofs)] +=
                    jxw * f[component] * velocity_values[i][q];

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj = pressure_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_total,
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                pressure_index(pj, n_velocity_dofs),
                                -jxw * pressure_values[j][q] *
                                  velocity_grads[i][q][component]);
                  }
              }

            for (int i = 0; i < n_local_pressure_dofs; ++i)
              {
                const int pi_dof = pressure_element_dofs[i];
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_total,
                                pressure_index(pi_dof, n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                -jxw * pressure_values[i][q] *
                                  velocity_grads[j][q][component]);
                  }
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (velocity_space.dofBoundaryMark(i) != 0)
        {
          const auto &point = velocity_space.dofInfo(i).interp_point;
          const double p[dim] = {point[0], point[1]};
          const double values[dim] = {colliding_velocity_x_exact(p),
                                      colliding_velocity_y_exact(p)};
          for (int component = 0; component < dim; ++component)
            apply_dirichlet_constraint(matrix,
                                       rhs,
                                       n_total,
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs),
                                       values[component]);
        }

    const auto &pressure_point = pressure_space.dofInfo(0).interp_point;
    const double pressure_pin_point[dim] = {pressure_point[0],
                                            pressure_point[1]};
    apply_dirichlet_constraint(matrix,
                               rhs,
                               n_total,
                               pressure_index(0, n_velocity_dofs),
                               colliding_pressure_exact(pressure_pin_point));

    solve_dense_system(matrix, rhs, n_total);

    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);
    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        u_x(i) = rhs[velocity_index(0, i, n_velocity_dofs)];
        u_y(i) = rhs[velocity_index(1, i, n_velocity_dofs)];
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      p_h(i) = rhs[pressure_index(i, n_velocity_dofs)];

    u_x.writeOpenDXData(output_prefix + "_ux.dx");
    u_y.writeOpenDXData(output_prefix + "_uy.dx");
    p_h.writeOpenDXData(output_prefix + "_p.dx");

    FunctionFunction<double> exact_u_x(&colliding_velocity_x_exact,
                                       &colliding_velocity_x_gradient);
    FunctionFunction<double> exact_u_y(&colliding_velocity_y_exact,
                                       &colliding_velocity_y_gradient);

    const double ux_l2 = Functional::L2Error(u_x, exact_u_x, 6);
    const double uy_l2 = Functional::L2Error(u_y, exact_u_y, 6);
    const double ux_h1 = Functional::H1SemiError(u_x, exact_u_x, 6);
    const double uy_h1 = Functional::H1SemiError(u_y, exact_u_y, 6);

    return {h,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            pressure_l2_error_mod_constant(p_h),
            divergence_l2_norm(u_x, u_y)};
  }
}

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output-prefix]]\n";
      return 2;
    }

  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes colliding-flow manufactured solution "
                   "(book Example 5.1.4)\n";
      std::cout << "The polynomial solution is mapped from (-1,1)^2 "
                   "to the unit square; pressure is scaled so f=0.\n";

      if (argc > 1)
        {
          const std::string mesh_file = argv[1];
          const std::string output_prefix =
            argc > 2 ? argv[2] :
                       std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                         "/stokes_colliding_afepack";
          const CollidingResult result =
            solve_colliding_stokes(mesh_file, output_prefix, 0.0);

          std::cout << "cells: " << result.cells << '\n';
          std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
          std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
          std::cout << "velocity L2 error: "
                    << result.velocity_l2_error << '\n';
          std::cout << "velocity H1-semi error: "
                    << result.velocity_h1_error << '\n';
          std::cout << "mean-free pressure L2 error: "
                    << result.pressure_l2_error << '\n';
          std::cout << "divergence L2 norm: "
                    << result.divergence_l2_norm << '\n';
          std::cout << "Wrote " << output_prefix << "_ux.dx, "
                    << output_prefix << "_uy.dx, and "
                    << output_prefix << "_p.dx\n";
          return 0;
        }

      const std::vector<std::pair<std::string, double>> meshes = {
        {"unit_square_h0p20", 0.20},
        {"unit_square_h0p10", 0.10},
      };

      std::vector<CollidingResult> results;
      for (const auto &[mesh_name, h] : meshes)
        {
          const std::string mesh_file =
            std::string(AFEPACK_DEFAULT_MESH_DIR) + "/" + mesh_name;
          const std::string output_prefix =
            std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
            "/stokes_colliding_afepack_" + mesh_name;
          results.push_back(solve_colliding_stokes(mesh_file,
                                                   output_prefix,
                                                   h));
        }

      std::cout << "h cells velocity_dofs pressure_dofs "
                   "velocity_L2 velocity_L2_rate velocity_H1 velocity_H1_rate "
                   "pressure_L2 pressure_L2_rate divergence_L2\n";
      for (std::size_t i = 0; i < results.size(); ++i)
        {
          const auto &r = results[i];
          const double l2_rate =
            i == 0 ? 0.0 :
                     rate(results[i - 1].velocity_l2_error,
                          r.velocity_l2_error,
                          results[i - 1].h,
                          r.h);
          const double h1_rate =
            i == 0 ? 0.0 :
                     rate(results[i - 1].velocity_h1_error,
                          r.velocity_h1_error,
                          results[i - 1].h,
                          r.h);
          const double p_rate =
            i == 0 ? 0.0 :
                     rate(results[i - 1].pressure_l2_error,
                          r.pressure_l2_error,
                          results[i - 1].h,
                          r.h);
          std::cout << r.h << ' '
                    << r.cells << ' '
                    << r.velocity_dofs << ' '
                    << r.pressure_dofs << ' '
                    << r.velocity_l2_error << ' '
                    << l2_rate << ' '
                    << r.velocity_h1_error << ' '
                    << h1_rate << ' '
                    << r.pressure_l2_error << ' '
                    << p_rate << ' '
                    << r.divergence_l2_norm << '\n';
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
