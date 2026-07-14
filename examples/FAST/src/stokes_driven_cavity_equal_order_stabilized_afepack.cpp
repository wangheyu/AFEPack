/**
 * @file stokes_driven_cavity_equal_order_stabilized_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：等阶稳定化。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 等阶稳定化：使用速度和压力等阶空间并加入稳定化项以满足离散稳定性要求。
 * - 驱动方腔流：经典顶盖驱动方腔基准问题，用于验证速度-压力耦合和涡结构。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <utility>

#define main stokes_poiseuille_afepack_hidden_main
#include "stokes_poiseuille_afepack.cpp"
#undef main

namespace
{
  constexpr double stabilization_parameter = 1.0;

  struct StabilizedCavityResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    double       velocity_l2_norm = 0.0;
    double       velocity_h1_seminorm = 0.0;
    double       max_speed = 0.0;
    double       pressure_min = 0.0;
    double       pressure_max = 0.0;
    double       pressure_mean = 0.0;
    double       pressure_centerline_max = 0.0;
    double       pressure_antisymmetry_rms = 0.0;
    double       divergence_l2_norm = 0.0;
  };

  double original_x(const double *p)
  {
    return 2.0 * p[0] - 1.0;
  }

  double lid_velocity_x(const double *p)
  {
    const double x = original_x(p);
    return 1.0 - x * x * x * x;
  }

  bool near(const double value, const double target)
  {
    return std::abs(value - target) < boundary_tolerance;
  }

  double cavity_velocity_boundary_value(const int component, const double *p)
  {
    if (component == 0 && near(p[1], 1.0))
      return lid_velocity_x(p);
    return 0.0;
  }

  int nearest_pressure_dof(const FEMSpace<double, dim> &pressure_space,
                           const double x,
                           const double y)
  {
    int    nearest = 0;
    double nearest_distance_square = std::numeric_limits<double>::max();
    for (int i = 0; i < pressure_space.n_dof(); ++i)
      {
        const auto &point = pressure_space.dofInfo(i).interp_point;
        const double dx = point[0] - x;
        const double dy = point[1] - y;
        const double distance_square = dx * dx + dy * dy;
        if (distance_square < nearest_distance_square)
          {
            nearest = i;
            nearest_distance_square = distance_square;
          }
      }
    return nearest;
  }

  double pressure_integral(const FEMFunction<double, dim> &pressure)
  {
    const auto &space = pressure.femSpace();
    double integral = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(4);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto values = pressure.value(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          integral += quad_info.weight(q) * jacobian[q] * volume * values[q];
      }

    return integral;
  }

  double pressure_centerline_max(const FEMFunction<double, dim> &pressure)
  {
    double max_abs = 0.0;
    const auto &space = pressure.femSpace();
    for (int i = 0; i < space.n_dof(); ++i)
      {
        const auto &point = space.dofInfo(i).interp_point;
        if (near(point[0], 0.5))
          max_abs = std::max(max_abs, std::abs(pressure(i)));
      }
    return max_abs;
  }

  double pressure_antisymmetry_rms(const FEMFunction<double, dim> &pressure)
  {
    const auto &space = pressure.femSpace();
    double sum_square = 0.0;
    int    count = 0;
    for (int i = 0; i < space.n_dof(); ++i)
      {
        const auto &point = space.dofInfo(i).interp_point;
        const int mirror =
          nearest_pressure_dof(space, 1.0 - point[0], point[1]);
        const double defect = pressure(i) + pressure(mirror);
        sum_square += defect * defect;
        ++count;
      }
    return std::sqrt(sum_square / std::max(count, 1));
  }

  double max_speed(const FEMFunction<double, dim> &u_x,
                   const FEMFunction<double, dim> &u_y)
  {
    double maximum = 0.0;
    for (std::size_t i = 0; i < u_x.size(); ++i)
      {
        const double speed =
          std::sqrt(u_x(i) * u_x(i) + u_y(i) * u_y(i));
        maximum = std::max(maximum, speed);
      }
    return maximum;
  }

  StabilizedCavityResult solve_stabilized_cavity(
    const std::string &mesh_file,
    const std::string &output_prefix,
    const double h)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> velocity_geometry;
    CoordTransform<dim, dim> velocity_transform;
    TemplateDOF<dim> velocity_dof;
    BasisFunctionAdmin<double, dim, dim> velocity_basis(velocity_dof);
    auto velocity_template = make_triangle_template(1,
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
          velocity_element.findQuadratureInfo(5);
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

        double element_measure = 0.0;
        std::vector<double> pressure_mean(n_local_pressure_dofs, 0.0);
        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            element_measure += jxw;
            for (int i = 0; i < n_local_pressure_dofs; ++i)
              pressure_mean[i] += jxw * pressure_values[i][q];
          }
        for (double &mean : pressure_mean)
          mean /= element_measure;

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;

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

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj = pressure_element_dofs[j];
                    const double stabilized_i =
                      pressure_values[i][q] - pressure_mean[i];
                    const double stabilized_j =
                      pressure_values[j][q] - pressure_mean[j];
                    add_entry(matrix,
                              n_total,
                              pressure_index(pi_dof, n_velocity_dofs),
                              pressure_index(pj, n_velocity_dofs),
                              -stabilization_parameter * jxw *
                                stabilized_i * stabilized_j);
                  }
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (velocity_space.dofBoundaryMark(i) != 0)
        {
          const auto &point = velocity_space.dofInfo(i).interp_point;
          const double p[dim] = {point[0], point[1]};
          for (int component = 0; component < dim; ++component)
            apply_dirichlet_constraint(matrix,
                                       rhs,
                                       n_total,
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs),
                                       cavity_velocity_boundary_value(component,
                                                                      p));
        }

    const int pressure_pin = nearest_pressure_dof(pressure_space, 0.5, 0.5);
    apply_dirichlet_constraint(matrix,
                               rhs,
                               n_total,
                               pressure_index(pressure_pin, n_velocity_dofs),
                               0.0);

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

    const double ux_l2 = Functional::L2Norm(u_x, 5);
    const double uy_l2 = Functional::L2Norm(u_y, 5);
    const double ux_h1 = Functional::H1Seminorm(u_x, 5);
    const double uy_h1 = Functional::H1Seminorm(u_y, 5);

    double pressure_min = std::numeric_limits<double>::max();
    double pressure_max = -std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < p_h.size(); ++i)
      {
        pressure_min = std::min(pressure_min, p_h(i));
        pressure_max = std::max(pressure_max, p_h(i));
      }

    return {h,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            max_speed(u_x, u_y),
            pressure_min,
            pressure_max,
            pressure_integral(p_h),
            pressure_centerline_max(p_h),
            pressure_antisymmetry_rms(p_h),
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
      std::cout << "AFEPack Stokes P1/P1 local-projection stabilized "
                   "regularized driven cavity test\n";
      std::cout << "Book Example 5.1.3 with top lid u_x=1-(2x-1)^4; "
                   "pressure local-projection stabilization parameter beta="
                << stabilization_parameter << ".\n";

      if (argc > 1)
        {
          const std::string mesh_file = argv[1];
          const std::string output_prefix =
            argc > 2 ? argv[2] :
                       std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                         "/stokes_driven_cavity_equal_order_stabilized_afepack";
          const StabilizedCavityResult result =
            solve_stabilized_cavity(mesh_file, output_prefix, 0.0);

          std::cout << "cells: " << result.cells << '\n';
          std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
          std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
          std::cout << "velocity L2 norm: "
                    << result.velocity_l2_norm << '\n';
          std::cout << "velocity H1-semi norm: "
                    << result.velocity_h1_seminorm << '\n';
          std::cout << "max nodal speed: " << result.max_speed << '\n';
          std::cout << "pressure min/max: " << result.pressure_min << " "
                    << result.pressure_max << '\n';
          std::cout << "pressure mean: " << result.pressure_mean << '\n';
          std::cout << "pressure centerline max abs: "
                    << result.pressure_centerline_max << '\n';
          std::cout << "pressure antisymmetry RMS: "
                    << result.pressure_antisymmetry_rms << '\n';
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

      std::vector<StabilizedCavityResult> results;
      for (const auto &[mesh_name, h] : meshes)
        {
          const std::string mesh_file =
            std::string(AFEPACK_DEFAULT_MESH_DIR) + "/" + mesh_name;
          const std::string output_prefix =
            std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
            "/stokes_driven_cavity_equal_order_stabilized_afepack_" +
            mesh_name;
          results.push_back(solve_stabilized_cavity(mesh_file,
                                                    output_prefix,
                                                    h));
        }

      std::cout << "h cells velocity_dofs pressure_dofs "
                   "velocity_L2 velocity_H1 max_speed pressure_min "
                   "pressure_max pressure_centerline_max "
                   "pressure_antisymmetry_RMS divergence_L2\n";
      for (const auto &r : results)
        {
          std::cout << r.h << ' '
                    << r.cells << ' '
                    << r.velocity_dofs << ' '
                    << r.pressure_dofs << ' '
                    << r.velocity_l2_norm << ' '
                    << r.velocity_h1_seminorm << ' '
                    << r.max_speed << ' '
                    << r.pressure_min << ' '
                    << r.pressure_max << ' '
                    << r.pressure_centerline_max << ' '
                    << r.pressure_antisymmetry_rms << ' '
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
