/**
 * @file navier_stokes_lsc_viscosity_sweep_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：黏性系数扫描。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 黏性系数扫描：批量改变黏性系数或雷诺数，观察预条件器和非线性迭代的稳健性。
 * - LSC 预条件器：构造最小二乘交换子 Schur 补近似，比较其稳健性和代价。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>

#include <lapacke.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef AFEPACK_DEFAULT_MESH_DIR
#  define AFEPACK_DEFAULT_MESH_DIR "build/meshes/afepack"
#endif

#ifndef AFEPACK_DEFAULT_OUTPUT_DIR
#  define AFEPACK_DEFAULT_OUTPUT_DIR "build"
#endif

namespace
{
  constexpr int    dim = 2;
  constexpr double pi  = 3.141592653589793238462643383279502884;
  constexpr double hp_regularization = 1.0e-8;
  constexpr int    max_picard_iterations = 40;
  constexpr int    max_gmres_iterations = 5000;
  constexpr double picard_tolerance = 1.0e-3;
  constexpr double lsc_linear_tolerance = 1.0e-3;

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct NavierStokesLSCResult
  {
    double       viscosity = 0.0;
    double       reynolds_number = 0.0;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    int          picard_iterations = 0;
    int          total_gmres_iterations = 0;
    double       final_relative_update = 0.0;
    double       final_gmres_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       divergence_l2_norm = 0.0;
  };

  class DenseLU
  {
  public:
    DenseLU() = default;

    DenseLU(const std::vector<double> &matrix, const int size)
    {
      factor(matrix, size);
    }

    void factor(const std::vector<double> &matrix, const int size)
    {
      n = size;
      lu = matrix;
      pivots.assign(n, 0);

      const int info =
        LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, lu.data(), n, pivots.data());
      if (info != 0)
        throw std::runtime_error("LAPACKE_dgetrf failed with info=" +
                                 std::to_string(info));
    }

    std::vector<double> solve(const std::vector<double> &rhs) const
    {
      if (static_cast<int>(rhs.size()) != n)
        throw std::runtime_error("DenseLU solve size mismatch");

      std::vector<double> solution = rhs;
      const int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,
                                      'N',
                                      n,
                                      1,
                                      lu.data(),
                                      n,
                                      pivots.data(),
                                      solution.data(),
                                      1);
      if (info != 0)
        throw std::runtime_error("LAPACKE_dgetrs failed with info=" +
                                 std::to_string(info));

      return solution;
    }

  private:
    int n = 0;
    std::vector<double> lu;
    std::vector<int> pivots;
  };

  double sx(const double *p) { return std::sin(pi * p[0]); }
  double cx(const double *p) { return std::cos(pi * p[0]); }
  double sy(const double *p) { return std::sin(pi * p[1]); }
  double cy(const double *p) { return std::cos(pi * p[1]); }
  double c2x(const double *p) { return std::cos(2.0 * pi * p[0]); }
  double c2y(const double *p) { return std::cos(2.0 * pi * p[1]); }

  double velocity_x_exact(const double *p)
  {
    return 2.0 * pi * sx(p) * sx(p) * sy(p) * cy(p);
  }

  double velocity_y_exact(const double *p)
  {
    return -2.0 * pi * sx(p) * cx(p) * sy(p) * sy(p);
  }

  std::vector<double> velocity_x_gradient(const double *p)
  {
    return {4.0 * pi * pi * sx(p) * cx(p) * sy(p) * cy(p),
            2.0 * pi * pi * sx(p) * sx(p) * c2y(p)};
  }

  std::vector<double> velocity_y_gradient(const double *p)
  {
    return {-2.0 * pi * pi * c2x(p) * sy(p) * sy(p),
            -4.0 * pi * pi * sx(p) * cx(p) * sy(p) * cy(p)};
  }

  double pressure_exact(const double *p)
  {
    return std::cos(pi * p[0]) * std::sin(pi * p[1]);
  }

  std::vector<double> pressure_gradient(const double *p)
  {
    return {-pi * sx(p) * sy(p), pi * cx(p) * cy(p)};
  }

  double laplace_velocity_x(const double *p)
  {
    const double b = sy(p) * cy(p);
    return 4.0 * pi * pi * pi * b * (c2x(p) - 2.0 * sx(p) * sx(p));
  }

  double laplace_velocity_y(const double *p)
  {
    const double a = sx(p) * cx(p);
    return 4.0 * pi * pi * pi * a * (2.0 * sy(p) * sy(p) - c2y(p));
  }

  double oseen_convection_x(const double *p)
  {
    const auto grad_x = velocity_x_gradient(p);
    return velocity_x_exact(p) * grad_x[0] +
           velocity_y_exact(p) * grad_x[1];
  }

  double oseen_convection_y(const double *p)
  {
    const auto grad_y = velocity_y_gradient(p);
    return velocity_x_exact(p) * grad_y[0] +
           velocity_y_exact(p) * grad_y[1];
  }

  double rhs_x(const double *p, const double viscosity)
  {
    return -viscosity * laplace_velocity_x(p) + oseen_convection_x(p) +
           pressure_gradient(p)[0];
  }

  double rhs_y(const double *p, const double viscosity)
  {
    return -viscosity * laplace_velocity_y(p) + oseen_convection_y(p) +
           pressure_gradient(p)[1];
  }

  std::string default_template_path()
  {
    if (const char *env = std::getenv("AFEPACK_PATH"))
      return std::string(env) + "/template/triangle";
    return "/home/steve/local/AFEPack/include/AFEPack/template/triangle";
  }

  void ensure_template_path()
  {
    if (std::getenv("AFEPACK_TEMPLATE_PATH") != nullptr)
      return;

    const std::string path = default_template_path();
    if (::setenv("AFEPACK_TEMPLATE_PATH", path.c_str(), 0) != 0)
      throw std::runtime_error("failed to set AFEPACK_TEMPLATE_PATH");
  }

  std::vector<TemplateElement<double, dim, dim>>
  make_triangle_template(const int order,
                         TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle." + std::to_string(order) + ".tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle." + std::to_string(order) + ".bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(geometry, dof, transform, basis);
    return template_element;
  }

  FEMSpace<double, dim> build_space(
    EasyMesh &mesh,
    std::vector<TemplateElement<double, dim, dim>> &template_element)
  {
    FEMSpace<double, dim> fem_space;
    fem_space.reinit(mesh, template_element);

    const int n_element = mesh.n_geometry(dim);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
      fem_space.element(i).reinit(fem_space, i, 0);

    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    return fem_space;
  }

  int velocity_index(const int component,
                     const int dof,
                     const int n_velocity_dofs)
  {
    return component * n_velocity_dofs + dof;
  }

  int pressure_index(const int dof, const int n_velocity_dofs)
  {
    return 2 * n_velocity_dofs + dof;
  }

  double vector_dot(const std::vector<double> &a,
                    const std::vector<double> &b)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
      result += a[i] * b[i];
    return result;
  }

  double vector_norm(const std::vector<double> &v)
  {
    return std::sqrt(vector_dot(v, v));
  }

  void add_scaled(std::vector<double> &destination,
                  const double coefficient,
                  const std::vector<double> &source)
  {
    for (std::size_t i = 0; i < destination.size(); ++i)
      destination[i] += coefficient * source[i];
  }

  std::vector<double> dense_vmult(const std::vector<double> &matrix,
                                  const std::vector<double> &x,
                                  const int n)
  {
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        y[i] += matrix[static_cast<std::size_t>(i) * n + j] * x[j];
    return y;
  }

  std::vector<double> dense_vmult_rect(const std::vector<double> &matrix,
                                       const std::vector<double> &x,
                                       const int n_rows,
                                       const int n_cols)
  {
    std::vector<double> y(n_rows, 0.0);
    for (int i = 0; i < n_rows; ++i)
      for (int j = 0; j < n_cols; ++j)
        y[i] += matrix[static_cast<std::size_t>(i) * n_cols + j] * x[j];
    return y;
  }

  std::vector<double> dense_tvmult_rect(const std::vector<double> &matrix,
                                        const std::vector<double> &x,
                                        const int n_rows,
                                        const int n_cols)
  {
    std::vector<double> y(n_cols, 0.0);
    for (int i = 0; i < n_rows; ++i)
      for (int j = 0; j < n_cols; ++j)
        y[j] += matrix[static_cast<std::size_t>(i) * n_cols + j] * x[i];
    return y;
  }

  void add_entry(std::vector<double> &matrix,
                 const int n,
                 const int row,
                 const int col,
                 const double value)
  {
    matrix[static_cast<std::size_t>(row) * n + col] += value;
  }

  void apply_homogeneous_constraint(std::vector<double> &matrix,
                                    std::vector<double> &rhs,
                                    const int n,
                                    const int index)
  {
    for (int j = 0; j < n; ++j)
      {
        matrix[static_cast<std::size_t>(index) * n + j] = 0.0;
        matrix[static_cast<std::size_t>(j) * n + index] = 0.0;
      }
    matrix[static_cast<std::size_t>(index) * n + index] = 1.0;
    rhs[index] = 0.0;
  }

  void apply_homogeneous_matrix_constraint(std::vector<double> &matrix,
                                           const int n,
                                           const int index)
  {
    for (int j = 0; j < n; ++j)
      {
        matrix[static_cast<std::size_t>(index) * n + j] = 0.0;
        matrix[static_cast<std::size_t>(j) * n + index] = 0.0;
      }
    matrix[static_cast<std::size_t>(index) * n + index] = 1.0;
  }

  void subtract_arithmetic_mean(std::vector<double> &v)
  {
    double sum = 0.0;
    for (const double value : v)
      sum += value;

    const double mean = sum / static_cast<double>(v.size());
    for (double &value : v)
      value -= mean;
  }

  void apply_inverse_velocity_mass(std::vector<double> &vector,
                                   const std::vector<double> &mass_diagonal)
  {
    for (std::size_t i = 0; i < vector.size(); ++i)
      vector[i] /= mass_diagonal[i];
  }

  std::vector<double> apply_velocity_block(
    const std::vector<double> &velocity_block,
    const std::vector<double> &x,
    const int n_velocity_dofs)
  {
    std::vector<double> y(2 * n_velocity_dofs, 0.0);
    std::vector<double> component_x(n_velocity_dofs);

    for (int component = 0; component < dim; ++component)
      {
        for (int i = 0; i < n_velocity_dofs; ++i)
          component_x[i] = x[velocity_index(component, i, n_velocity_dofs)];

        const auto component_y =
          dense_vmult(velocity_block, component_x, n_velocity_dofs);

        for (int i = 0; i < n_velocity_dofs; ++i)
          y[velocity_index(component, i, n_velocity_dofs)] = component_y[i];
      }

    return y;
  }

  std::vector<double> apply_hp(const std::vector<double> &divergence,
                               const std::vector<double> &velocity_mass_diagonal,
                               const std::vector<double> &src,
                               const int n_pressure_dofs,
                               const int n_velocity_block_dofs)
  {
    auto velocity_tmp = dense_tvmult_rect(divergence,
                                          src,
                                          n_pressure_dofs,
                                          n_velocity_block_dofs);
    apply_inverse_velocity_mass(velocity_tmp, velocity_mass_diagonal);
    return dense_vmult_rect(divergence,
                            velocity_tmp,
                            n_pressure_dofs,
                            n_velocity_block_dofs);
  }

  std::vector<double> apply_gp(const std::vector<double> &velocity_block,
                               const std::vector<double> &divergence,
                               const std::vector<double> &velocity_mass_diagonal,
                               const std::vector<double> &src,
                               const int n_velocity_dofs,
                               const int n_pressure_dofs)
  {
    const int n_velocity_block_dofs = 2 * n_velocity_dofs;
    auto velocity_tmp_1 = dense_tvmult_rect(divergence,
                                            src,
                                            n_pressure_dofs,
                                            n_velocity_block_dofs);
    apply_inverse_velocity_mass(velocity_tmp_1, velocity_mass_diagonal);

    auto velocity_tmp_2 =
      apply_velocity_block(velocity_block, velocity_tmp_1, n_velocity_dofs);
    apply_inverse_velocity_mass(velocity_tmp_2, velocity_mass_diagonal);

    return dense_vmult_rect(divergence,
                            velocity_tmp_2,
                            n_pressure_dofs,
                            n_velocity_block_dofs);
  }

  std::vector<double> build_hp_matrix(
    const std::vector<double> &divergence,
    const std::vector<double> &velocity_mass_diagonal,
    const int n_velocity_dofs,
    const int n_pressure_dofs)
  {
    const int n_velocity_block_dofs = 2 * n_velocity_dofs;
    std::vector<double> matrix(
      static_cast<std::size_t>(n_pressure_dofs) * n_pressure_dofs,
      0.0);

    for (int j = 0; j < n_pressure_dofs; ++j)
      {
        std::vector<double> basis_vector(n_pressure_dofs, 0.0);
        basis_vector[j] = 1.0;

        auto column = apply_hp(divergence,
                               velocity_mass_diagonal,
                               basis_vector,
                               n_pressure_dofs,
                               n_velocity_block_dofs);
        for (int i = 0; i < n_pressure_dofs; ++i)
          matrix[static_cast<std::size_t>(i) * n_pressure_dofs + j] =
            column[i];
      }

    for (int i = 0; i < n_pressure_dofs; ++i)
      matrix[static_cast<std::size_t>(i) * n_pressure_dofs + i] +=
        hp_regularization;

    return matrix;
  }

  std::vector<double> build_gp_matrix(
    const std::vector<double> &velocity_block,
    const std::vector<double> &divergence,
    const std::vector<double> &velocity_mass_diagonal,
    const int n_velocity_dofs,
    const int n_pressure_dofs)
  {
    std::vector<double> matrix(
      static_cast<std::size_t>(n_pressure_dofs) * n_pressure_dofs,
      0.0);

    for (int j = 0; j < n_pressure_dofs; ++j)
      {
        std::vector<double> basis_vector(n_pressure_dofs, 0.0);
        basis_vector[j] = 1.0;

        auto column = apply_gp(velocity_block,
                               divergence,
                               velocity_mass_diagonal,
                               basis_vector,
                               n_velocity_dofs,
                               n_pressure_dofs);
        for (int i = 0; i < n_pressure_dofs; ++i)
          matrix[static_cast<std::size_t>(i) * n_pressure_dofs + j] =
            column[i];
      }

    return matrix;
  }

  class BlockLSCPreconditioner
  {
  public:
    BlockLSCPreconditioner(const int n_velocity_dofs,
                           const int n_pressure_dofs,
                           const std::vector<double> &velocity_block,
                           const std::vector<double> &divergence,
                           const std::vector<double> &velocity_mass_diagonal)
      : n_velocity_dofs(n_velocity_dofs),
        n_pressure_dofs(n_pressure_dofs),
        velocity_solver(velocity_block, n_velocity_dofs),
        hp(build_hp_matrix(divergence,
                           velocity_mass_diagonal,
                           n_velocity_dofs,
                           n_pressure_dofs)),
        gp(build_gp_matrix(velocity_block,
                           divergence,
                           velocity_mass_diagonal,
                           n_velocity_dofs,
                           n_pressure_dofs)),
        hp_solver(hp, n_pressure_dofs)
    {}

    std::vector<double> vmult(const std::vector<double> &src) const
    {
      std::vector<double> dst(src.size(), 0.0);
      std::vector<double> component_rhs(n_velocity_dofs);

      for (int component = 0; component < dim; ++component)
        {
          for (int i = 0; i < n_velocity_dofs; ++i)
            component_rhs[i] =
              src[velocity_index(component, i, n_velocity_dofs)];

          const auto component_solution = velocity_solver.solve(component_rhs);
          for (int i = 0; i < n_velocity_dofs; ++i)
            dst[velocity_index(component, i, n_velocity_dofs)] =
              component_solution[i];
        }

      std::vector<double> pressure_rhs(n_pressure_dofs);
      for (int i = 0; i < n_pressure_dofs; ++i)
        pressure_rhs[i] =
          src[pressure_index(i, n_velocity_dofs)];
      subtract_arithmetic_mean(pressure_rhs);

      auto pressure_tmp = hp_solver.solve(pressure_rhs);
      pressure_tmp = dense_vmult(gp, pressure_tmp, n_pressure_dofs);
      subtract_arithmetic_mean(pressure_tmp);
      auto pressure_solution = hp_solver.solve(pressure_tmp);

      for (int i = 0; i < n_pressure_dofs; ++i)
        dst[pressure_index(i, n_velocity_dofs)] = pressure_solution[i];

      return dst;
    }

  private:
    int n_velocity_dofs = 0;
    int n_pressure_dofs = 0;
    DenseLU velocity_solver;
    std::vector<double> hp;
    std::vector<double> gp;
    DenseLU hp_solver;
  };

  void apply_plane_rotation(double &dx,
                            double &dy,
                            const double cs,
                            const double sn)
  {
    const double temporary = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temporary;
  }

  void generate_plane_rotation(const double dx,
                               const double dy,
                               double &cs,
                               double &sn)
  {
    if (dy == 0.0)
      {
        cs = 1.0;
        sn = 0.0;
      }
    else if (std::abs(dy) > std::abs(dx))
      {
        const double temporary = dx / dy;
        sn = 1.0 / std::sqrt(1.0 + temporary * temporary);
        cs = temporary * sn;
      }
    else
      {
        const double temporary = dy / dx;
        cs = 1.0 / std::sqrt(1.0 + temporary * temporary);
        sn = temporary * cs;
      }
  }

  std::vector<double> solve_upper_hessenberg(
    const std::vector<std::vector<double>> &h,
    const std::vector<double> &g,
    const int size)
  {
    std::vector<double> y(size, 0.0);
    for (int i = size - 1; i >= 0; --i)
      {
        double sum = g[i];
        for (int j = i + 1; j < size; ++j)
          sum -= h[i][j] * y[j];
        y[i] = sum / h[i][i];
      }
    return y;
  }

  template <typename Preconditioner>
  SolveInfo solve_right_preconditioned_gmres(
    const std::vector<double> &matrix,
    std::vector<double> &x,
    const std::vector<double> &rhs,
    const int n,
    const Preconditioner &preconditioner,
    const double tolerance,
    const int restart,
    const int max_iterations)
  {
    const double rhs_norm = std::max(vector_norm(rhs), 1.0e-30);
    auto ax = dense_vmult(matrix, x, n);
    std::vector<double> r = rhs;
    add_scaled(r, -1.0, ax);

    double beta = vector_norm(r);
    double relative_residual = beta / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    int total_iterations = 0;
    while (total_iterations < max_iterations)
      {
        std::vector<std::vector<double>> v(restart + 1,
                                           std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> z(restart,
                                           std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> h(
          restart + 1, std::vector<double>(restart, 0.0));
        std::vector<double> cs(restart, 0.0);
        std::vector<double> sn(restart, 0.0);
        std::vector<double> g(restart + 1, 0.0);

        v[0] = r;
        for (double &entry : v[0])
          entry /= beta;
        g[0] = beta;

        int inner_iterations = 0;
        for (; inner_iterations < restart &&
               total_iterations < max_iterations;
             ++inner_iterations, ++total_iterations)
          {
            z[inner_iterations] =
              preconditioner.vmult(v[inner_iterations]);
            v[inner_iterations + 1] =
              dense_vmult(matrix, z[inner_iterations], n);

            for (int i = 0; i <= inner_iterations; ++i)
              {
                h[i][inner_iterations] =
                  vector_dot(v[inner_iterations + 1], v[i]);
                add_scaled(v[inner_iterations + 1],
                           -h[i][inner_iterations],
                           v[i]);
              }

            h[inner_iterations + 1][inner_iterations] =
              vector_norm(v[inner_iterations + 1]);
            if (h[inner_iterations + 1][inner_iterations] != 0.0)
              for (double &entry : v[inner_iterations + 1])
                entry /= h[inner_iterations + 1][inner_iterations];

            for (int i = 0; i < inner_iterations; ++i)
              apply_plane_rotation(h[i][inner_iterations],
                                   h[i + 1][inner_iterations],
                                   cs[i],
                                   sn[i]);

            generate_plane_rotation(h[inner_iterations][inner_iterations],
                                    h[inner_iterations + 1][inner_iterations],
                                    cs[inner_iterations],
                                    sn[inner_iterations]);
            apply_plane_rotation(h[inner_iterations][inner_iterations],
                                 h[inner_iterations + 1][inner_iterations],
                                 cs[inner_iterations],
                                 sn[inner_iterations]);
            apply_plane_rotation(g[inner_iterations],
                                 g[inner_iterations + 1],
                                 cs[inner_iterations],
                                 sn[inner_iterations]);

            relative_residual =
              std::abs(g[inner_iterations + 1]) / rhs_norm;
            if (relative_residual <= tolerance)
              {
                const int system_size = inner_iterations + 1;
                const auto y =
                  solve_upper_hessenberg(h, g, system_size);
                for (int i = 0; i < system_size; ++i)
                  add_scaled(x, y[i], z[i]);
                return {true, total_iterations + 1, relative_residual};
              }
          }

        const auto y =
          solve_upper_hessenberg(h, g, inner_iterations);
        for (int i = 0; i < inner_iterations; ++i)
          add_scaled(x, y[i], z[i]);

        ax = dense_vmult(matrix, x, n);
        r = rhs;
        add_scaled(r, -1.0, ax);
        beta = vector_norm(r);
        relative_residual = beta / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }

  struct AssembledOseenSystem
  {
    std::vector<double> matrix;
    std::vector<double> rhs;
    std::vector<double> velocity_block;
    std::vector<double> divergence;
    std::vector<double> velocity_mass_diagonal;
  };

  AssembledOseenSystem assemble_oseen_system(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    const FEMFunction<double, dim> &wind_x,
    const FEMFunction<double, dim> &wind_y,
    const double viscosity)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    AssembledOseenSystem system;
    system.matrix.assign(static_cast<std::size_t>(n_total) * n_total, 0.0);
    system.rhs.assign(n_total, 0.0);
    system.velocity_block.assign(
      static_cast<std::size_t>(n_velocity_dofs) * n_velocity_dofs, 0.0);
    system.divergence.assign(
      static_cast<std::size_t>(n_pressure_dofs) * (2 * n_velocity_dofs),
      0.0);
    system.velocity_mass_diagonal.assign(2 * n_velocity_dofs, 0.0);

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
        const auto wind_x_values =
          wind_x.value(q_point, velocity_element);
        const auto wind_y_values =
          wind_y.value(q_point, velocity_element);

        const auto &velocity_element_dofs = velocity_element.dof();
        const auto &pressure_element_dofs = pressure_element.dof();
        const int n_local_velocity_dofs = velocity_element_dofs.size();
        const int n_local_pressure_dofs = pressure_element_dofs.size();

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double p[dim] = {q_point[q][0], q_point[q][1]};
            const double wind[dim] = {wind_x_values[q],
                                      wind_y_values[q]};
            const double f[dim] = {rhs_x(p, viscosity),
                                   rhs_y(p, viscosity)};

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
                double lumped_mass = 0.0;
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    const double diffusion =
                      viscosity * innerProduct(velocity_grads[i][q],
                                               velocity_grads[j][q]);
                    const double convection =
                      (wind[0] * velocity_grads[j][q][0] +
                       wind[1] * velocity_grads[j][q][1]) *
                      velocity_values[i][q];
                    const double velocity_entry =
                      jxw * (diffusion + convection);
                    lumped_mass +=
                      jxw * velocity_values[i][q] * velocity_values[j][q];

                    add_entry(system.velocity_block,
                              n_velocity_dofs,
                              vi,
                              vj,
                              velocity_entry);
                    for (int component = 0; component < dim; ++component)
                      add_entry(system.matrix,
                                n_total,
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                velocity_entry);
                  }

                for (int component = 0; component < dim; ++component)
                  system.velocity_mass_diagonal[velocity_index(
                    component, vi, n_velocity_dofs)] += lumped_mass;

                for (int component = 0; component < dim; ++component)
                  system.rhs[velocity_index(component,
                                            vi,
                                            n_velocity_dofs)] +=
                    jxw * f[component] * velocity_values[i][q];

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj = pressure_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(system.matrix,
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
                      {
                        const double divergence_entry =
                          -jxw * pressure_values[i][q] *
                          velocity_grads[j][q][component];
                        add_entry(system.matrix,
                                  n_total,
                                  pressure_index(pi_dof, n_velocity_dofs),
                                  velocity_index(component,
                                                 vj,
                                                 n_velocity_dofs),
                                  divergence_entry);
                        add_entry(system.divergence,
                                  2 * n_velocity_dofs,
                                  pi_dof,
                                  velocity_index(component,
                                                 vj,
                                                 n_velocity_dofs),
                                  divergence_entry);
                      }
                  }
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (velocity_space.dofBoundaryMark(i) != 0)
        {
          apply_homogeneous_matrix_constraint(system.velocity_block,
                                             n_velocity_dofs,
                                             i);
          for (int component = 0; component < dim; ++component)
            apply_homogeneous_constraint(system.matrix,
                                         system.rhs,
                                         n_total,
                                         velocity_index(component,
                                                        i,
                                                        n_velocity_dofs));
        }

    apply_homogeneous_constraint(system.matrix,
                                 system.rhs,
                                 n_total,
                                 pressure_index(0, n_velocity_dofs));

    for (double &entry : system.velocity_mass_diagonal)
      if (std::abs(entry) < 1.0e-14)
        entry = 1.0;

    return system;
  }

  void assign_solution(const std::vector<double> &solution,
                       FEMFunction<double, dim> &u_x,
                       FEMFunction<double, dim> &u_y,
                       FEMFunction<double, dim> &p_h)
  {
    const int n_velocity_dofs = u_x.femSpace().n_dof();
    const int n_pressure_dofs = p_h.femSpace().n_dof();

    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        u_x(i) = solution[velocity_index(0, i, n_velocity_dofs)];
        u_y(i) = solution[velocity_index(1, i, n_velocity_dofs)];
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      p_h(i) = solution[pressure_index(i, n_velocity_dofs)];
  }

  void set_zero(FEMFunction<double, dim> &function)
  {
    for (std::size_t i = 0; i < function.size(); ++i)
      function(i) = 0.0;
  }

  void copy_velocity(FEMFunction<double, dim> &destination_x,
                     FEMFunction<double, dim> &destination_y,
                     const FEMFunction<double, dim> &source_x,
                     const FEMFunction<double, dim> &source_y)
  {
    for (std::size_t i = 0; i < destination_x.size(); ++i)
      {
        destination_x(i) = source_x(i);
        destination_y(i) = source_y(i);
      }
  }

  double relative_velocity_update(const FEMFunction<double, dim> &old_x,
                                  const FEMFunction<double, dim> &old_y,
                                  const FEMFunction<double, dim> &new_x,
                                  const FEMFunction<double, dim> &new_y)
  {
    double update_square = 0.0;
    double solution_square = 0.0;

    for (std::size_t i = 0; i < new_x.size(); ++i)
      {
        const double dx = new_x(i) - old_x(i);
        const double dy = new_y(i) - old_y(i);
        update_square += dx * dx + dy * dy;
        solution_square += new_x(i) * new_x(i) + new_y(i) * new_y(i);
      }

    return std::sqrt(update_square) /
           std::max(std::sqrt(solution_square), 1.0e-30);
  }

  double divergence_l2_norm(const FEMFunction<double, dim> &u_x,
                            const FEMFunction<double, dim> &u_y)
  {
    const auto &space = u_x.femSpace();
    double norm_square = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(5);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto grad_x = u_x.gradient(q_point, element);
        const auto grad_y = u_y.gradient(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double div_u = grad_x[q][0] + grad_y[q][1];
            norm_square += jxw * div_u * div_u;
          }
      }

    return std::sqrt(norm_square);
  }

  std::vector<NavierStokesLSCResult> solve_navier_stokes_lsc_sweep(
    const std::string &mesh_file,
    const std::string &output_prefix,
    const std::vector<double> &viscosities)
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

    FEMFunction<double, dim> old_x(velocity_space);
    FEMFunction<double, dim> old_y(velocity_space);
    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);

    set_zero(old_x);
    set_zero(old_y);
    set_zero(u_x);
    set_zero(u_y);
    set_zero(p_h);

    FunctionFunction<double> exact_u_x(&velocity_x_exact,
                                       &velocity_x_gradient);
    FunctionFunction<double> exact_u_y(&velocity_y_exact,
                                       &velocity_y_gradient);
    FunctionFunction<double> exact_p(&pressure_exact,
                                     &pressure_gradient);

    std::vector<NavierStokesLSCResult> results;
    results.reserve(viscosities.size());

    for (const double nu : viscosities)
      {
        double final_relative_update = 1.0;
        double final_gmres_relative_residual = 1.0;
        int picard_iteration = 0;
        int total_gmres_iterations = 0;

        for (; picard_iteration < max_picard_iterations; ++picard_iteration)
          {
            auto system = assemble_oseen_system(velocity_space,
                                                pressure_space,
                                                old_x,
                                                old_y,
                                                nu);
            BlockLSCPreconditioner preconditioner(
              n_velocity_dofs,
              n_pressure_dofs,
              system.velocity_block,
              system.divergence,
              system.velocity_mass_diagonal);

            std::vector<double> solution(n_total, 0.0);
            const SolveInfo gmres_info =
              solve_right_preconditioned_gmres(system.matrix,
                                               solution,
                                               system.rhs,
                                               n_total,
                                               preconditioner,
                                               lsc_linear_tolerance,
                                               80,
                                               max_gmres_iterations);
            if (!gmres_info.converged)
              throw std::runtime_error(
                "right-preconditioned LSC GMRES failed for viscosity " +
                std::to_string(nu) + " in Picard step " +
                std::to_string(picard_iteration + 1) + " after " +
                std::to_string(gmres_info.iterations) +
                " iterations; final relative residual=" +
                std::to_string(gmres_info.relative_residual));

            assign_solution(solution, u_x, u_y, p_h);
            final_relative_update =
              relative_velocity_update(old_x, old_y, u_x, u_y);
            final_gmres_relative_residual = gmres_info.relative_residual;
            total_gmres_iterations += gmres_info.iterations;
            ++picard_iteration;

            if (final_relative_update < picard_tolerance)
              break;

            copy_velocity(old_x, old_y, u_x, u_y);
          }

        if (final_relative_update >= picard_tolerance)
          throw std::runtime_error(
            "Picard iteration did not converge for viscosity " +
            std::to_string(nu) + " within " +
            std::to_string(max_picard_iterations) +
            " iterations; final relative update=" +
            std::to_string(final_relative_update));

        const double ux_l2 = Functional::L2Error(u_x, exact_u_x, 5);
        const double uy_l2 = Functional::L2Error(u_y, exact_u_y, 5);
        const double ux_h1 = Functional::H1SemiError(u_x, exact_u_x, 5);
        const double uy_h1 = Functional::H1SemiError(u_y, exact_u_y, 5);

        results.push_back({nu,
                           1.0 / nu,
                           static_cast<int>(mesh.n_geometry(dim)),
                           velocity_space.n_dof(),
                           pressure_space.n_dof(),
                           picard_iteration,
                           total_gmres_iterations,
                           final_relative_update,
                           final_gmres_relative_residual,
                           std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
                           std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
                           Functional::L2Error(p_h, exact_p, 5),
                           divergence_l2_norm(u_x, u_y)});

        copy_velocity(old_x, old_y, u_x, u_y);
      }

    u_x.writeOpenDXData(output_prefix + "_ux.dx");
    u_y.writeOpenDXData(output_prefix + "_uy.dx");
    p_h.writeOpenDXData(output_prefix + "_p.dx");

    return results;
  }

  std::vector<double> parse_viscosities(const std::string &argument)
  {
    if (argument.empty())
      return {1.0, 0.5, 0.25, 0.125};

    std::vector<double> viscosities;
    std::stringstream   stream(argument);
    std::string         token;
    while (std::getline(stream, token, ','))
      {
        const double viscosity = std::stod(token);
        if (!(viscosity > 0.0))
          throw std::runtime_error("viscosity values must be positive");
        viscosities.push_back(viscosity);
      }
    if (viscosities.empty())
      throw std::runtime_error("viscosity list is empty");
    return viscosities;
  }
}

int main(int argc, char **argv)
{
  if (argc > 4)
    {
      std::cerr << "usage: " << argv[0]
                << " [mesh-file [output-prefix [viscosities]]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
                 std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p20";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/navier_stokes_lsc_viscosity_sweep_afepack";
  const auto viscosities =
    parse_viscosities(argc > 3 ? argv[3] : std::string());

  try
    {
      const auto results =
        solve_navier_stokes_lsc_sweep(mesh_file, output_prefix, viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes Picard LSC viscosity continuation "
                   "Taylor-Hood manufactured solution\n";
      if (!results.empty())
        {
          std::cout << "cells: " << results.front().cells << '\n';
          std::cout << "velocity dofs: "
                    << results.front().velocity_dofs << '\n';
          std::cout << "pressure dofs: "
                    << results.front().pressure_dofs << '\n';
        }
      std::cout << "nu Re Picard GMRES GMRES/Picard update linear_residual "
                   "velocity_L2 velocity_H1 pressure_L2 divergence_L2\n";
      for (const auto &result : results)
        std::cout << result.viscosity << ' '
                  << result.reynolds_number << ' '
                  << result.picard_iterations << ' '
                  << result.total_gmres_iterations << ' '
                  << static_cast<double>(result.total_gmres_iterations) /
                       std::max(result.picard_iterations, 1)
                  << ' '
                  << result.final_relative_update << ' '
                  << result.final_gmres_relative_residual << ' '
                  << result.velocity_l2_error << ' '
                  << result.velocity_h1_error << ' '
                  << result.pressure_l2_error << ' '
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
