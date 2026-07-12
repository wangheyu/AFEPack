/**
 * @file stokes_minres_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：MINRES 求解。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - MINRES 求解：面向对称鞍点系统组织矩阵向量乘和块预条件器，验证 MINRES 收敛性。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <lapacke.h>

#define main stokes_block_iterative_afepack_hidden_main
#include "stokes_block_iterative_afepack.cpp"
#undef main

namespace
{
  struct MinresStokesResult
  {
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    unsigned int system_nonzeros = 0;
    int          minres_iterations = 0;
    double       scaled_relative_residual = 0.0;
    double       original_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       divergence_l2_norm = 0.0;
  };

  std::vector<double> make_inverse_sqrt_block_diagonal(
    const RowSparseMatrix &matrix,
    const std::vector<double> &pressure_mass_diag,
    const int n_velocity_dofs)
  {
    std::vector<double> scale(matrix.size(), 1.0);

    for (int i = 0; i < 2 * n_velocity_dofs; ++i)
      {
        const double diagonal = std::abs(matrix.diagonal(i));
        scale[i] = diagonal > 1.0e-14 ? 1.0 / std::sqrt(diagonal) : 1.0;
      }

    for (std::size_t i = 0; i < pressure_mass_diag.size(); ++i)
      {
        const int row = pressure_index(static_cast<int>(i),
                                       n_velocity_dofs);
        const double mass = std::abs(pressure_mass_diag[i]);
        scale[row] = mass > 1.0e-14 ? 1.0 / std::sqrt(mass) : 1.0;
      }

    return scale;
  }

  std::vector<double> apply_symmetric_scaled_operator(
    const RowSparseMatrix &matrix,
    const std::vector<double> &inverse_sqrt_diagonal,
    const std::vector<double> &src)
  {
    std::vector<double> tmp(src.size(), 0.0);
    for (std::size_t i = 0; i < src.size(); ++i)
      tmp[i] = inverse_sqrt_diagonal[i] * src[i];

    tmp = matrix.vmult(tmp);

    for (std::size_t i = 0; i < tmp.size(); ++i)
      tmp[i] *= inverse_sqrt_diagonal[i];
    return tmp;
  }

  std::vector<double> solve_lanczos_least_squares(
    const std::vector<double> &alpha,
    const std::vector<double> &subdiagonal_beta,
    const double initial_residual_norm)
  {
    const int k = static_cast<int>(alpha.size());
    const int m = k + 1;
    std::vector<double> tbar(static_cast<std::size_t>(m) * k, 0.0);
    std::vector<double> rhs(m, 0.0);
    rhs[0] = initial_residual_norm;

    for (int j = 0; j < k; ++j)
      {
        if (j > 0)
          tbar[static_cast<std::size_t>(j - 1) * k + j] =
            subdiagonal_beta[j - 1];
        tbar[static_cast<std::size_t>(j) * k + j] = alpha[j];
        tbar[static_cast<std::size_t>(j + 1) * k + j] =
          subdiagonal_beta[j];
      }

    const int info = LAPACKE_dgels(LAPACK_ROW_MAJOR,
                                   'N',
                                   m,
                                   k,
                                   1,
                                   tbar.data(),
                                   k,
                                   rhs.data(),
                                   1);
    if (info != 0)
      throw std::runtime_error("LAPACKE_dgels failed in MINRES with info=" +
                               std::to_string(info));

    rhs.resize(k);
    return rhs;
  }

  std::vector<double> combine_basis(
    const std::vector<std::vector<double>> &basis,
    const std::vector<double> &coefficients)
  {
    std::vector<double> result(basis.front().size(), 0.0);
    for (std::size_t j = 0; j < coefficients.size(); ++j)
      add_scaled(result, coefficients[j], basis[j]);
    return result;
  }

  SolveInfo solve_symmetric_scaled_minres(
    const RowSparseMatrix &matrix,
    std::vector<double> &x,
    const std::vector<double> &rhs,
    const std::vector<double> &inverse_sqrt_diagonal,
    const double tolerance,
    const int max_iterations)
  {
    std::vector<double> scaled_rhs(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i)
      scaled_rhs[i] = inverse_sqrt_diagonal[i] * rhs[i];

    const double scaled_rhs_norm =
      std::max(vector_norm(scaled_rhs), 1.0e-30);
    const double rhs_norm = std::max(vector_norm(rhs), 1.0e-30);

    std::vector<std::vector<double>> basis;
    std::vector<double> v_current = scaled_rhs;
    const double initial_residual_norm = vector_norm(v_current);
    if (initial_residual_norm / scaled_rhs_norm <= tolerance)
      return {true, 0, initial_residual_norm / scaled_rhs_norm};

    for (double &entry : v_current)
      entry /= initial_residual_norm;
    basis.push_back(v_current);

    std::vector<double> v_previous(rhs.size(), 0.0);
    std::vector<double> alpha;
    std::vector<double> subdiagonal_beta;

    double previous_beta = 0.0;
    double scaled_relative_residual = 1.0;
    double original_relative_residual = 1.0;

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        std::vector<double> w =
          apply_symmetric_scaled_operator(matrix,
                                          inverse_sqrt_diagonal,
                                          v_current);
        if (iteration > 1)
          add_scaled(w, -previous_beta, v_previous);

        const double current_alpha = vector_dot(v_current, w);
        add_scaled(w, -current_alpha, v_current);

        const double next_beta = vector_norm(w);
        alpha.push_back(current_alpha);
        subdiagonal_beta.push_back(next_beta);

        const auto coefficients =
          solve_lanczos_least_squares(alpha,
                                      subdiagonal_beta,
                                      initial_residual_norm);
        const auto scaled_solution = combine_basis(basis, coefficients);

        std::vector<double> scaled_residual = scaled_rhs;
        auto h_y = apply_symmetric_scaled_operator(matrix,
                                                   inverse_sqrt_diagonal,
                                                   scaled_solution);
        add_scaled(scaled_residual, -1.0, h_y);
        scaled_relative_residual =
          vector_norm(scaled_residual) / scaled_rhs_norm;

        for (std::size_t i = 0; i < x.size(); ++i)
          x[i] = inverse_sqrt_diagonal[i] * scaled_solution[i];

        std::vector<double> original_residual = rhs;
        auto ax = matrix.vmult(x);
        add_scaled(original_residual, -1.0, ax);
        original_relative_residual =
          vector_norm(original_residual) / rhs_norm;

        if (original_relative_residual <= tolerance)
          return {true, iteration, scaled_relative_residual};

        if (next_beta <= 1.0e-14)
          break;

        v_previous = v_current;
        v_current = w;
        for (double &entry : v_current)
          entry /= next_beta;
        basis.push_back(v_current);
        previous_beta = next_beta;
      }

    return {false, max_iterations, scaled_relative_residual};
  }

  MinresStokesResult solve_stokes_minres(const std::string &mesh_file,
                                         const std::string &output_prefix)
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

    RowSparseMatrix matrix(n_total);
    std::vector<double> rhs(n_total, 0.0);
    std::vector<double> pressure_mass_diag(n_pressure_dofs, 0.0);

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
                      matrix.add(velocity_index(component,
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
                      matrix.add(velocity_index(component,
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
                      matrix.add(pressure_index(pi_dof, n_velocity_dofs),
                                 velocity_index(component,
                                                vj,
                                                n_velocity_dofs),
                                 -jxw * pressure_values[i][q] *
                                   velocity_grads[j][q][component]);
                  }

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  if (pi_dof == pressure_element_dofs[j])
                    pressure_mass_diag[pi_dof] +=
                      jxw * pressure_values[i][q] * pressure_values[j][q];
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (velocity_space.dofBoundaryMark(i) != 0)
        for (int component = 0; component < dim; ++component)
          apply_homogeneous_constraint(matrix,
                                       rhs,
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs));

    apply_homogeneous_constraint(matrix,
                                 rhs,
                                 pressure_index(0, n_velocity_dofs));
    pressure_mass_diag[0] = 1.0;

    const auto inverse_sqrt_diagonal =
      make_inverse_sqrt_block_diagonal(matrix,
                                       pressure_mass_diag,
                                       n_velocity_dofs);

    std::vector<double> solution(n_total, 0.0);
    const SolveInfo minres_info =
      solve_symmetric_scaled_minres(matrix,
                                    solution,
                                    rhs,
                                    inverse_sqrt_diagonal,
                                    1.0e-8,
                                    1000);

    std::vector<double> original_residual = rhs;
    auto ax = matrix.vmult(solution);
    add_scaled(original_residual, -1.0, ax);
    const double original_relative_residual =
      vector_norm(original_residual) / std::max(vector_norm(rhs), 1.0e-30);

    if (!minres_info.converged)
      throw std::runtime_error(
        "symmetric scaled MINRES failed after " +
        std::to_string(minres_info.iterations) +
        " iterations; final original relative residual=" +
        std::to_string(original_relative_residual));

    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);
    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        u_x(i) = solution[velocity_index(0, i, n_velocity_dofs)];
        u_y(i) = solution[velocity_index(1, i, n_velocity_dofs)];
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      p_h(i) = solution[pressure_index(i, n_velocity_dofs)];

    u_x.writeOpenDXData(output_prefix + "_ux.dx");
    u_y.writeOpenDXData(output_prefix + "_uy.dx");
    p_h.writeOpenDXData(output_prefix + "_p.dx");

    FunctionFunction<double> exact_u_x(&velocity_x_exact,
                                       &velocity_x_gradient);
    FunctionFunction<double> exact_u_y(&velocity_y_exact,
                                       &velocity_y_gradient);
    FunctionFunction<double> exact_p(&pressure_exact,
                                     &pressure_gradient);

    const double ux_l2 = Functional::L2Error(u_x, exact_u_x, 5);
    const double uy_l2 = Functional::L2Error(u_y, exact_u_y, 5);
    const double ux_h1 = Functional::H1SemiError(u_x, exact_u_x, 5);
    const double uy_h1 = Functional::H1SemiError(u_y, exact_u_y, 5);

    return {static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            matrix.nonzeros(),
            minres_info.iterations,
            minres_info.relative_residual,
            original_relative_residual,
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            Functional::L2Error(p_h, exact_p, 5),
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

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p20";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/stokes_minres_afepack";

  try
    {
      const MinresStokesResult result =
        solve_stokes_minres(mesh_file, output_prefix);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes symmetric scaled MINRES Taylor-Hood "
                   "manufactured solution\n";
      std::cout << "cells: " << result.cells << '\n';
      std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
      std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
      std::cout << "system nonzeros: " << result.system_nonzeros << '\n';
      std::cout << "MINRES iterations: " << result.minres_iterations << '\n';
      std::cout << "scaled relative residual: "
                << result.scaled_relative_residual << '\n';
      std::cout << "original relative residual: "
                << result.original_relative_residual << '\n';
      std::cout << "velocity L2 error: " << result.velocity_l2_error << '\n';
      std::cout << "velocity H1-semi error: "
                << result.velocity_h1_error << '\n';
      std::cout << "pressure L2 error with fixed first pressure dof: "
                << result.pressure_l2_error << '\n';
      std::cout << "divergence L2 norm: "
                << result.divergence_l2_norm << '\n';
      std::cout << "Wrote " << output_prefix << "_ux.dx, "
                << output_prefix << "_uy.dx, and "
                << output_prefix << "_p.dx\n";
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
