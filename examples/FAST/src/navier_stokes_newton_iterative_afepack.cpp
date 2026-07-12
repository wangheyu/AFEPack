/**
 * @file navier_stokes_newton_iterative_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：迭代线性求解。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 迭代线性求解：使用 Krylov 迭代和块预条件器替代直接求解，以观察可扩展性。
 * - Newton 线性化：使用 Newton 迭代求解非线性 Navier-Stokes 系统并比较收敛特性。
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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
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
  constexpr double viscosity = 1.0;
  constexpr int    max_newton_iterations = 12;
  constexpr double newton_tolerance = 1.0e-8;

  struct NavierStokesResult
  {
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    unsigned int system_nonzeros = 0;
    int          newton_iterations = 0;
    int          total_gmres_iterations = 0;
    double       final_relative_correction = 0.0;
    double       final_gmres_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       divergence_l2_norm = 0.0;
  };

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
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

  double exact_convection_x(const double *p)
  {
    const auto grad_x = velocity_x_gradient(p);
    return velocity_x_exact(p) * grad_x[0] +
           velocity_y_exact(p) * grad_x[1];
  }

  double exact_convection_y(const double *p)
  {
    const auto grad_y = velocity_y_gradient(p);
    return velocity_x_exact(p) * grad_y[0] +
           velocity_y_exact(p) * grad_y[1];
  }

  double rhs_x(const double *p)
  {
    return -viscosity * laplace_velocity_x(p) + exact_convection_x(p) +
           pressure_gradient(p)[0];
  }

  double rhs_y(const double *p)
  {
    return -viscosity * laplace_velocity_y(p) + exact_convection_y(p) +
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

  class RowSparseMatrix
  {
  public:
    explicit RowSparseMatrix(const int n = 0)
      : n_rows(n), rows(n)
    {}

    void reinit(const int n)
    {
      n_rows = n;
      rows.assign(n, {});
    }

    int size() const
    {
      return n_rows;
    }

    void add(const int row, const int col, const double value)
    {
      if (std::abs(value) == 0.0)
        return;
      rows[row][col] += value;
    }

    void clear_row_and_column(const int index)
    {
      for (auto &row : rows)
        row.erase(index);
      rows[index].clear();
      rows[index][index] = 1.0;
    }

    std::vector<double> vmult(const std::vector<double> &x) const
    {
      std::vector<double> y(n_rows, 0.0);
      for (int i = 0; i < n_rows; ++i)
        for (const auto &[j, value] : rows[i])
          y[i] += value * x[j];
      return y;
    }

    double diagonal(const int index) const
    {
      const auto it = rows[index].find(index);
      return it == rows[index].end() ? 0.0 : it->second;
    }

    unsigned int nonzeros() const
    {
      unsigned int nnz = 0;
      for (const auto &row : rows)
        nnz += static_cast<unsigned int>(row.size());
      return nnz;
    }

  private:
    int n_rows = 0;
    std::vector<std::map<int, double>> rows;
  };

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

  void apply_homogeneous_constraint(RowSparseMatrix &matrix,
                                    std::vector<double> &rhs,
                                    const int index)
  {
    matrix.clear_row_and_column(index);
    rhs[index] = 0.0;
  }

  class BlockDiagonalPreconditioner
  {
  public:
    BlockDiagonalPreconditioner(const RowSparseMatrix &matrix,
                                const std::vector<double> &pressure_mass_diag,
                                const int n_velocity_dofs)
      : inverse_diagonal(matrix.size(), 1.0)
    {
      for (int i = 0; i < 2 * n_velocity_dofs; ++i)
        {
          const double diagonal = matrix.diagonal(i);
          inverse_diagonal[i] =
            std::abs(diagonal) > 1.0e-14 ? 1.0 / diagonal : 1.0;
        }

      for (std::size_t i = 0; i < pressure_mass_diag.size(); ++i)
        {
          const int row = pressure_index(static_cast<int>(i),
                                         n_velocity_dofs);
          const double mass = pressure_mass_diag[i];
          inverse_diagonal[row] =
            std::abs(mass) > 1.0e-14 ? 1.0 / mass : 1.0;
        }
    }

    std::vector<double> vmult(const std::vector<double> &src) const
    {
      std::vector<double> dst(src.size(), 0.0);
      for (std::size_t i = 0; i < src.size(); ++i)
        dst[i] = inverse_diagonal[i] * src[i];
      return dst;
    }

  private:
    std::vector<double> inverse_diagonal;
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
    const RowSparseMatrix &matrix,
    std::vector<double> &x,
    const std::vector<double> &rhs,
    const Preconditioner &preconditioner,
    const double tolerance,
    const int restart,
    const int max_iterations)
  {
    const int n = matrix.size();
    const double rhs_norm = std::max(vector_norm(rhs), 1.0e-30);
    auto ax = matrix.vmult(x);
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
              matrix.vmult(z[inner_iterations]);

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

        ax = matrix.vmult(x);
        r = rhs;
        add_scaled(r, -1.0, ax);
        beta = vector_norm(r);
        relative_residual = beta / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }

  void assemble_newton_correction_system(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    const FEMFunction<double, dim> &current_x,
    const FEMFunction<double, dim> &current_y,
    const FEMFunction<double, dim> &current_p,
    RowSparseMatrix &matrix,
    std::vector<double> &rhs,
    std::vector<double> &pressure_mass_diag)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    matrix.reinit(n_total);
    rhs.assign(n_total, 0.0);
    pressure_mass_diag.assign(n_pressure_dofs, 0.0);

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
        const auto current_x_values =
          current_x.value(q_point, velocity_element);
        const auto current_y_values =
          current_y.value(q_point, velocity_element);
        const auto current_x_grads =
          current_x.gradient(q_point, velocity_element);
        const auto current_y_grads =
          current_y.gradient(q_point, velocity_element);
        const auto current_pressure_values =
          current_p.value(q_point, pressure_element);

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
            const double current_velocity[dim] = {current_x_values[q],
                                                  current_y_values[q]};
            const double current_grad[dim][dim] = {
              {current_x_grads[q][0], current_x_grads[q][1]},
              {current_y_grads[q][0], current_y_grads[q][1]}};
            const double current_pressure = current_pressure_values[q];
            const double current_divergence =
              current_grad[0][0] + current_grad[1][1];
            const double current_convection[dim] = {
              current_velocity[0] * current_grad[0][0] +
                current_velocity[1] * current_grad[0][1],
              current_velocity[0] * current_grad[1][0] +
                current_velocity[1] * current_grad[1][1]};

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
                for (int component = 0; component < dim; ++component)
                  {
                    const double current_diffusion =
                      viscosity *
                      (current_grad[component][0] *
                         velocity_grads[i][q][0] +
                       current_grad[component][1] *
                         velocity_grads[i][q][1]);
                    const double current_pressure_term =
                      -current_pressure *
                      velocity_grads[i][q][component];
                    const double residual =
                      current_diffusion +
                      current_convection[component] *
                        velocity_values[i][q] +
                      current_pressure_term -
                      f[component] * velocity_values[i][q];

                    rhs[velocity_index(component, vi, n_velocity_dofs)] -=
                      jxw * residual;
                  }

                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];

                    for (int test_component = 0; test_component < dim;
                         ++test_component)
                      for (int trial_component = 0; trial_component < dim;
                           ++trial_component)
                        {
                          double entry = 0.0;

                          if (test_component == trial_component)
                            {
                              const double diffusion =
                                viscosity *
                                innerProduct(velocity_grads[i][q],
                                             velocity_grads[j][q]);
                              const double advected_correction =
                                (current_velocity[0] *
                                   velocity_grads[j][q][0] +
                                 current_velocity[1] *
                                   velocity_grads[j][q][1]) *
                                velocity_values[i][q];
                              entry += diffusion + advected_correction;
                            }

                          const double newton_convection =
                            velocity_values[j][q] *
                            current_grad[test_component][trial_component] *
                            velocity_values[i][q];
                          entry += newton_convection;

                          matrix.add(velocity_index(test_component,
                                                    vi,
                                                    n_velocity_dofs),
                                     velocity_index(trial_component,
                                                    vj,
                                                    n_velocity_dofs),
                                     jxw * entry);
                        }
                  }

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
                rhs[pressure_index(pi_dof, n_velocity_dofs)] +=
                  jxw * pressure_values[i][q] * current_divergence;

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
  }

  double add_newton_correction(const std::vector<double> &correction,
                               FEMFunction<double, dim> &u_x,
                               FEMFunction<double, dim> &u_y,
                               FEMFunction<double, dim> &p_h)
  {
    const int n_velocity_dofs = u_x.femSpace().n_dof();
    const int n_pressure_dofs = p_h.femSpace().n_dof();

    double correction_square = 0.0;
    double velocity_square = 0.0;

    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        const double dx = correction[velocity_index(0, i, n_velocity_dofs)];
        const double dy = correction[velocity_index(1, i, n_velocity_dofs)];
        u_x(i) += dx;
        u_y(i) += dy;

        correction_square += dx * dx + dy * dy;
        velocity_square += u_x(i) * u_x(i) + u_y(i) * u_y(i);
      }

    for (int i = 0; i < n_pressure_dofs; ++i)
      p_h(i) += correction[pressure_index(i, n_velocity_dofs)];

    return std::sqrt(correction_square) /
           std::max(std::sqrt(velocity_square), 1.0e-30);
  }

  void set_zero(FEMFunction<double, dim> &function)
  {
    for (std::size_t i = 0; i < function.size(); ++i)
      function(i) = 0.0;
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

  NavierStokesResult solve_navier_stokes(const std::string &mesh_file,
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

    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);

    set_zero(u_x);
    set_zero(u_y);
    set_zero(p_h);

    RowSparseMatrix matrix(n_total);
    std::vector<double> rhs(n_total, 0.0);
    std::vector<double> correction(n_total, 0.0);
    std::vector<double> pressure_mass_diag(n_pressure_dofs, 0.0);
    unsigned int last_system_nonzeros = 0;
    int total_gmres_iterations = 0;
    double final_gmres_relative_residual = 1.0;
    double final_relative_correction = 1.0;
    int newton_iteration = 0;

    for (; newton_iteration < max_newton_iterations; ++newton_iteration)
      {
        assemble_newton_correction_system(velocity_space,
                                          pressure_space,
                                          u_x,
                                          u_y,
                                          p_h,
                                          matrix,
                                          rhs,
                                          pressure_mass_diag);

        last_system_nonzeros = matrix.nonzeros();
        correction.assign(n_total, 0.0);
        BlockDiagonalPreconditioner preconditioner(matrix,
                                                   pressure_mass_diag,
                                                   n_velocity_dofs);
        const SolveInfo gmres_info =
          solve_right_preconditioned_gmres(matrix,
                                           correction,
                                           rhs,
                                           preconditioner,
                                           1.0e-10,
                                           120,
                                           5000);
        if (!gmres_info.converged)
          throw std::runtime_error(
            "right-preconditioned GMRES failed in Newton step " +
            std::to_string(newton_iteration + 1) + " after " +
            std::to_string(gmres_info.iterations) +
            " iterations; final relative residual=" +
            std::to_string(gmres_info.relative_residual));

        total_gmres_iterations += gmres_info.iterations;
        final_gmres_relative_residual = gmres_info.relative_residual;
        final_relative_correction =
          add_newton_correction(correction, u_x, u_y, p_h);

        ++newton_iteration;

        if (final_relative_correction < newton_tolerance)
          break;
      }

    if (final_relative_correction >= newton_tolerance)
      throw std::runtime_error("Newton iteration did not converge within " +
                               std::to_string(max_newton_iterations) +
                               " iterations; final relative correction=" +
                               std::to_string(final_relative_correction));

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
            last_system_nonzeros,
            newton_iteration,
            total_gmres_iterations,
            final_relative_correction,
            final_gmres_relative_residual,
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
                 "/navier_stokes_newton_iterative_afepack";

  try
    {
      const NavierStokesResult result =
        solve_navier_stokes(mesh_file, output_prefix);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes sparse block GMRES Newton "
                   "Taylor-Hood manufactured solution\n";
      std::cout << "cells: " << result.cells << '\n';
      std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
      std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
      std::cout << "system nonzeros: " << result.system_nonzeros << '\n';
      std::cout << "Newton iterations: "
                << result.newton_iterations << '\n';
      std::cout << "total right-preconditioned GMRES iterations: "
                << result.total_gmres_iterations << '\n';
      std::cout << "final relative velocity correction: "
                << result.final_relative_correction << '\n';
      std::cout << "final GMRES relative residual: "
                << result.final_gmres_relative_residual << '\n';
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
