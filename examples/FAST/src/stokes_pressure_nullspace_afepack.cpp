/**
 * @file stokes_pressure_nullspace_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：压力零空间处理。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 压力零空间处理：显式处理压力常数核空间，避免鞍点系统奇异性影响迭代。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define main stokes_taylor_hood_afepack_hidden_main
#include "stokes_taylor_hood_afepack.cpp"
#undef main

namespace
{
  enum class PressureNormalization
  {
    fixed_first_dof,
    zero_mean_constraint
  };

  struct PressureStatistics
  {
    double area = 0.0;
    double numerical_mean = 0.0;
    double exact_mean = 0.0;
    double raw_l2_error = 0.0;
    double mean_free_l2_error = 0.0;
    double best_constant_shift = 0.0;
  };

  struct NullspaceResult
  {
    std::string  normalization;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    unsigned int linear_unknowns = 0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       pressure_mean_free_l2_error = 0.0;
    double       pressure_mean = 0.0;
    double       pressure_constraint_residual = 0.0;
    double       divergence_l2_norm = 0.0;
    double       lagrange_multiplier = 0.0;
  };

  const char *normalization_name(const PressureNormalization normalization)
  {
    switch (normalization)
      {
        case PressureNormalization::fixed_first_dof:
          return "fixed_first_pressure_dof";
        case PressureNormalization::zero_mean_constraint:
          return "zero_mean_pressure";
      }

    return "unknown";
  }

  PressureStatistics pressure_statistics(
    const FEMFunction<double, dim> &pressure)
  {
    const auto &space = pressure.femSpace();
    double area = 0.0;
    double numerical_integral = 0.0;
    double exact_integral = 0.0;
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
            const double exact_value = pressure_exact(point);
            const double error = values[q] - exact_value;

            area += jxw;
            numerical_integral += jxw * values[q];
            exact_integral += jxw * exact_value;
            error_integral += jxw * error;
            error_square_integral += jxw * error * error;
          }
      }

    const double best_constant_shift = error_integral / area;
    return {area,
            numerical_integral / area,
            exact_integral / area,
            std::sqrt(std::max(0.0, error_square_integral)),
            std::sqrt(std::max(0.0,
                               error_square_integral -
                                 area * best_constant_shift *
                                   best_constant_shift)),
            best_constant_shift};
  }

  double weighted_pressure_integral(
    const FEMFunction<double, dim> &pressure,
    const std::vector<double> &pressure_integral_weights)
  {
    double integral = 0.0;
    for (int i = 0; i < pressure.femSpace().n_dof(); ++i)
      integral += pressure_integral_weights[i] * pressure(i);
    return integral;
  }

  NullspaceResult solve_pressure_nullspace_case(
    const std::string &mesh_file,
    const PressureNormalization normalization)
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
    const bool use_zero_mean_constraint =
      normalization == PressureNormalization::zero_mean_constraint;
    const int n_linear = n_total + (use_zero_mean_constraint ? 1 : 0);
    const int pressure_constraint_index = n_total;

    std::vector<double> matrix(static_cast<std::size_t>(n_linear) *
                                 n_linear,
                               0.0);
    std::vector<double> rhs(n_linear, 0.0);
    std::vector<double> pressure_integral_weights(n_pressure_dofs, 0.0);

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

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double p[dim] = {q_point[q][0], q_point[q][1]};
            const double f[dim] = {rhs_x(p), rhs_y(p)};

            for (int i = 0; i < n_local_pressure_dofs; ++i)
              pressure_integral_weights[pressure_element_dofs[i]] +=
                jxw * pressure_values[i][q];

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
                                n_linear,
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
                                n_linear,
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
                                n_linear,
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
        for (int component = 0; component < dim; ++component)
          apply_homogeneous_constraint(matrix,
                                       rhs,
                                       n_linear,
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs));

    if (use_zero_mean_constraint)
      for (int i = 0; i < n_pressure_dofs; ++i)
        {
          const int index = pressure_index(i, n_velocity_dofs);
          add_entry(matrix,
                    n_linear,
                    index,
                    pressure_constraint_index,
                    pressure_integral_weights[i]);
          add_entry(matrix,
                    n_linear,
                    pressure_constraint_index,
                    index,
                    pressure_integral_weights[i]);
        }
    else
      apply_homogeneous_constraint(matrix,
                                   rhs,
                                   n_linear,
                                   pressure_index(0, n_velocity_dofs));

    solve_dense_system(matrix, rhs, n_linear);

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

    FunctionFunction<double> exact_u_x(&velocity_x_exact,
                                       &velocity_x_gradient);
    FunctionFunction<double> exact_u_y(&velocity_y_exact,
                                       &velocity_y_gradient);

    const double ux_l2 = Functional::L2Error(u_x, exact_u_x, 5);
    const double uy_l2 = Functional::L2Error(u_y, exact_u_y, 5);
    const double ux_h1 = Functional::H1SemiError(u_x, exact_u_x, 5);
    const double uy_h1 = Functional::H1SemiError(u_y, exact_u_y, 5);
    const PressureStatistics pressure = pressure_statistics(p_h);

    const double pressure_constraint_residual =
      use_zero_mean_constraint ?
        std::abs(weighted_pressure_integral(p_h,
                                            pressure_integral_weights)) :
        std::abs(p_h(0));
    const double lagrange_multiplier =
      use_zero_mean_constraint ? rhs[pressure_constraint_index] : 0.0;

    return {normalization_name(normalization),
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            static_cast<unsigned int>(n_linear),
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            pressure.raw_l2_error,
            pressure.mean_free_l2_error,
            pressure.numerical_mean,
            pressure_constraint_residual,
            divergence_l2_norm(u_x, u_y),
            lagrange_multiplier};
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p20";

  try
    {
      const std::vector<PressureNormalization> normalizations = {
        PressureNormalization::fixed_first_dof,
        PressureNormalization::zero_mean_constraint,
      };

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes pressure nullspace normalization check\n";
      std::cout << "mesh: " << mesh_file << '\n';
      std::cout << "strategy cells velocity_dofs pressure_dofs "
                   "linear_unknowns velocity_L2 velocity_H1 "
                   "pressure_L2 pressure_L2_mod_constants pressure_mean "
                   "constraint_residual divergence_L2 lagrange_multiplier\n";

      for (const auto normalization : normalizations)
        {
          const NullspaceResult result =
            solve_pressure_nullspace_case(mesh_file, normalization);
          std::cout << result.normalization << ' '
                    << result.cells << ' '
                    << result.velocity_dofs << ' '
                    << result.pressure_dofs << ' '
                    << result.linear_unknowns << ' '
                    << result.velocity_l2_error << ' '
                    << result.velocity_h1_error << ' '
                    << result.pressure_l2_error << ' '
                    << result.pressure_mean_free_l2_error << ' '
                    << result.pressure_mean << ' '
                    << result.pressure_constraint_residual << ' '
                    << result.divergence_l2_norm << ' '
                    << result.lagrange_multiplier << '\n';
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
