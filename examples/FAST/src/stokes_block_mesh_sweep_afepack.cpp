/**
 * @file stokes_block_mesh_sweep_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：块方法网格扫描。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 块方法网格扫描：在多组网格尺度上评估块预条件器随自由度增长的表现。
 * - 网格尺度扫描：在多组网格输入上重复求解，比较自由度增长对误差和迭代次数的影响。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/AMGSolver.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
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
  constexpr double viscosity = 1.0;

  struct StokesResult
  {
    std::string  preconditioner;
    double       mesh_h = 0.0;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    unsigned int system_nonzeros = 0;
    int          gmres_iterations = 0;
    double       gmres_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       divergence_l2_norm = 0.0;
    bool         has_rates = false;
    double       velocity_l2_rate = 0.0;
    double       velocity_h1_rate = 0.0;
    double       pressure_l2_rate = 0.0;
    double       divergence_rate = 0.0;
  };

  enum class PreconditionerKind
  {
    block_diagonal,
    schur_mass_pcg,
    schur_mass_amg1,
    schur_mass_amg2,
    schur_mass_amg3
  };

  struct MeshCase
  {
    double      h = 0.0;
    std::string mesh_file;
  };

  std::vector<PreconditionerKind>
  preconditioners_for_mesh(const MeshCase &mesh_case)
  {
    std::vector<PreconditionerKind> preconditioners = {
      PreconditionerKind::schur_mass_amg1,
      PreconditionerKind::schur_mass_amg2,
      PreconditionerKind::schur_mass_amg3};

    if (mesh_case.h >= 0.05)
      {
        preconditioners.insert(preconditioners.begin(),
                               PreconditionerKind::schur_mass_pcg);
        preconditioners.insert(preconditioners.begin(),
                               PreconditionerKind::block_diagonal);
      }

    return preconditioners;
  }

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

  double rhs_x(const double *p)
  {
    return -viscosity * laplace_velocity_x(p) + pressure_gradient(p)[0];
  }

  double rhs_y(const double *p)
  {
    return -viscosity * laplace_velocity_y(p) + pressure_gradient(p)[1];
  }

  double zero_boundary_value(const double *)
  {
    return 0.0;
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

    const std::map<int, double> &row_entries(const int index) const
    {
      return rows[index];
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

  std::string preconditioner_name(const PreconditionerKind kind)
  {
    switch (kind)
      {
        case PreconditionerKind::block_diagonal:
          return "diag(A)+diag(Mp)";
        case PreconditionerKind::schur_mass_pcg:
          return "tri(Apcg,Mp-pcg)";
        case PreconditionerKind::schur_mass_amg1:
          return "tri(Aamg1,Mp-pcg)";
        case PreconditionerKind::schur_mass_amg2:
          return "tri(Aamg2,Mp-pcg)";
        case PreconditionerKind::schur_mass_amg3:
          return "tri(Aamg3,Mp-pcg)";
      }
    return "unknown";
  }

  int velocity_amg_smooth_steps(const PreconditionerKind kind)
  {
    switch (kind)
      {
        case PreconditionerKind::schur_mass_amg1:
          return 1;
        case PreconditionerKind::schur_mass_amg2:
          return 2;
        case PreconditionerKind::schur_mass_amg3:
          return 3;
        default:
          return 0;
      }
  }

  double convergence_rate(const double coarse_h,
                          const double fine_h,
                          const double coarse_error,
                          const double fine_error)
  {
    if (coarse_h <= 0.0 || fine_h <= 0.0 ||
        coarse_error <= 0.0 || fine_error <= 0.0 ||
        coarse_h == fine_h)
      return 0.0;

    return std::log(coarse_error / fine_error) /
           std::log(coarse_h / fine_h);
  }

  void add_convergence_rates(std::vector<StokesResult> &results)
  {
    std::map<std::string, std::size_t> previous_result;

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        const auto previous = previous_result.find(results[i].preconditioner);
        if (previous != previous_result.end())
          {
            const StokesResult &coarse = results[previous->second];
            results[i].has_rates = true;
            results[i].velocity_l2_rate =
              convergence_rate(coarse.mesh_h,
                               results[i].mesh_h,
                               coarse.velocity_l2_error,
                               results[i].velocity_l2_error);
            results[i].velocity_h1_rate =
              convergence_rate(coarse.mesh_h,
                               results[i].mesh_h,
                               coarse.velocity_h1_error,
                               results[i].velocity_h1_error);
            results[i].pressure_l2_rate =
              convergence_rate(coarse.mesh_h,
                               results[i].mesh_h,
                               coarse.pressure_l2_error,
                               results[i].pressure_l2_error);
            results[i].divergence_rate =
              convergence_rate(coarse.mesh_h,
                               results[i].mesh_h,
                               coarse.divergence_l2_norm,
                               results[i].divergence_l2_norm);
          }

        previous_result[results[i].preconditioner] = i;
      }
  }

  std::string format_rate(const bool has_rate, const double rate)
  {
    if (!has_rate)
      return "-";

    std::ostringstream output;
    output << std::fixed << std::setprecision(3) << rate;
    return output.str();
  }

  void print_gmres_summary(const std::vector<StokesResult> &results)
  {
    struct Summary
    {
      double min_h = 0.0;
      double max_h = 0.0;
      int min_iterations = 0;
      int max_iterations = 0;
      bool has_rate = false;
      double last_velocity_h1_rate = 0.0;
      double last_pressure_l2_rate = 0.0;
    };

    std::map<std::string, Summary> summaries;
    for (const auto &result : results)
      {
        auto &summary = summaries[result.preconditioner];
        if (summary.max_h == 0.0)
          {
            summary.min_h = result.mesh_h;
            summary.max_h = result.mesh_h;
            summary.min_iterations = result.gmres_iterations;
            summary.max_iterations = result.gmres_iterations;
          }
        else
          {
            summary.min_h = std::min(summary.min_h, result.mesh_h);
            summary.max_h = std::max(summary.max_h, result.mesh_h);
            summary.min_iterations =
              std::min(summary.min_iterations, result.gmres_iterations);
            summary.max_iterations =
              std::max(summary.max_iterations, result.gmres_iterations);
          }

        if (result.has_rates)
          {
            summary.has_rate = true;
            summary.last_velocity_h1_rate = result.velocity_h1_rate;
            summary.last_pressure_l2_rate = result.pressure_l2_rate;
          }
      }

    std::cout << "\nGMRES mesh-independence summary\n";
    std::cout << "preconditioner h_min h_max GMRES_min GMRES_max "
                 "last_velocity_H1_rate last_pressure_L2_rate\n";
    for (const auto &[name, summary] : summaries)
      std::cout << std::setw(17) << name << ' '
                << summary.min_h << ' '
                << summary.max_h << ' '
                << summary.min_iterations << ' '
                << summary.max_iterations << ' '
                << std::setw(9)
                << format_rate(summary.has_rate,
                               summary.last_velocity_h1_rate) << ' '
                << std::setw(9)
                << format_rate(summary.has_rate,
                               summary.last_pressure_l2_rate) << '\n';
  }

  void apply_homogeneous_constraint(RowSparseMatrix &matrix,
                                    std::vector<double> &rhs,
                                    const int index)
  {
    matrix.clear_row_and_column(index);
    rhs[index] = 0.0;
  }

  void apply_homogeneous_constraint(RowSparseMatrix &matrix,
                                    const int index)
  {
    matrix.clear_row_and_column(index);
  }

  void apply_homogeneous_dirichlet_boundary(
    StiffMatrix<dim, double> &matrix,
    FEMFunction<double, dim> &solution,
    Vector<double> &rhs,
    const FEMSpace<double, dim> &space)
  {
    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &zero_boundary_value);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &zero_boundary_value);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &zero_boundary_value);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &zero_boundary_value);

    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.apply(matrix, solution, rhs);
  }

  std::vector<double> inverse_diagonal(const RowSparseMatrix &matrix)
  {
    std::vector<double> inverse(matrix.size(), 1.0);
    for (int i = 0; i < matrix.size(); ++i)
      {
        const double diagonal = matrix.diagonal(i);
        inverse[i] = std::abs(diagonal) > 1.0e-14 ? 1.0 / diagonal : 1.0;
      }
    return inverse;
  }

  SolveInfo solve_pcg(const RowSparseMatrix &matrix,
                      std::vector<double> &x,
                      const std::vector<double> &rhs,
                      const std::vector<double> &inverse_diagonal,
                      const double tolerance,
                      const int max_iterations)
  {
    const double rhs_norm = std::max(vector_norm(rhs), 1.0e-30);
    auto ax = matrix.vmult(x);
    std::vector<double> r = rhs;
    add_scaled(r, -1.0, ax);

    double relative_residual = vector_norm(r) / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    std::vector<double> z(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i)
      z[i] = inverse_diagonal[i] * r[i];
    std::vector<double> p = z;
    double rho = vector_dot(r, z);

    for (int iteration = 0; iteration < max_iterations; ++iteration)
      {
        const auto q = matrix.vmult(p);
        const double pq = vector_dot(p, q);
        if (std::abs(pq) <= 1.0e-30)
          break;

        const double alpha = rho / pq;
        add_scaled(x, alpha, p);
        add_scaled(r, -alpha, q);

        relative_residual = vector_norm(r) / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, iteration + 1, relative_residual};

        for (std::size_t i = 0; i < rhs.size(); ++i)
          z[i] = inverse_diagonal[i] * r[i];

        const double rho_new = vector_dot(r, z);
        if (std::abs(rho) <= 1.0e-30)
          break;

        const double beta = rho_new / rho;
        for (std::size_t i = 0; i < p.size(); ++i)
          p[i] = z[i] + beta * p[i];
        rho = rho_new;
      }

    return {false, max_iterations, relative_residual};
  }

  class BlockDiagonalPreconditioner
  {
  public:
    BlockDiagonalPreconditioner(const RowSparseMatrix &matrix,
                                const std::vector<double> &pressure_mass_diag,
                                const int n_velocity_dofs)
      : inverse_diagonal(matrix.size(), 1.0),
        n_velocity_dofs(n_velocity_dofs)
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
    int n_velocity_dofs = 0;
  };

  class SchurMassPreconditioner
  {
  public:
    SchurMassPreconditioner(const RowSparseMatrix &system_matrix,
                            const RowSparseMatrix &velocity_matrix,
                            const RowSparseMatrix &pressure_mass_matrix,
                            AMGPreconditioner *velocity_amg,
                            const int n_velocity_dofs,
                            const int n_pressure_dofs)
      : system_matrix(system_matrix),
        velocity_matrix(velocity_matrix),
        pressure_mass_matrix(pressure_mass_matrix),
        velocity_amg(velocity_amg),
        velocity_inverse_diagonal(inverse_diagonal(velocity_matrix)),
        pressure_inverse_diagonal(inverse_diagonal(pressure_mass_matrix)),
        n_velocity_dofs(n_velocity_dofs),
        n_pressure_dofs(n_pressure_dofs)
    {}

    std::vector<double> vmult(const std::vector<double> &src) const
    {
      std::vector<double> dst(src.size(), 0.0);

      std::vector<double> velocity_solution(2 * n_velocity_dofs, 0.0);
      for (int component = 0; component < dim; ++component)
        {
          std::vector<double> rhs(n_velocity_dofs, 0.0);
          for (int i = 0; i < n_velocity_dofs; ++i)
            rhs[i] = src[velocity_index(component, i, n_velocity_dofs)];

          std::vector<double> component_solution(n_velocity_dofs, 0.0);
          if (velocity_amg != nullptr)
            {
              Vector<double> amg_src(n_velocity_dofs);
              Vector<double> amg_dst(n_velocity_dofs);
              for (int i = 0; i < n_velocity_dofs; ++i)
                amg_src(i) = rhs[i];
              amg_dst = 0.0;
              velocity_amg->vmult(amg_dst, amg_src);
              for (int i = 0; i < n_velocity_dofs; ++i)
                component_solution[i] = amg_dst(i);
            }
          else
            {
              const SolveInfo info =
                solve_pcg(velocity_matrix,
                          component_solution,
                          rhs,
                          velocity_inverse_diagonal,
                          1.0e-2,
                          80);
              if (!info.converged)
                throw std::runtime_error(
                  "velocity block PCG failed inside Stokes preconditioner");
            }

          for (int i = 0; i < n_velocity_dofs; ++i)
            {
              velocity_solution[velocity_index(component,
                                               i,
                                               n_velocity_dofs)] =
                component_solution[i];
              dst[velocity_index(component, i, n_velocity_dofs)] =
                component_solution[i];
            }
        }

      std::vector<double> pressure_rhs(n_pressure_dofs, 0.0);
      for (int i = 0; i < n_pressure_dofs; ++i)
        {
          double by = 0.0;
          const int pressure_row = pressure_index(i, n_velocity_dofs);
          for (const auto &[column, value] :
               system_matrix.row_entries(pressure_row))
            if (column < 2 * n_velocity_dofs)
              by += value * velocity_solution[column];

          pressure_rhs[i] =
            src[pressure_index(i, n_velocity_dofs)] - by;
        }

      std::vector<double> pressure_solution(n_pressure_dofs, 0.0);
      const SolveInfo pressure_info =
        solve_pcg(pressure_mass_matrix,
                  pressure_solution,
                  pressure_rhs,
                  pressure_inverse_diagonal,
                  1.0e-10,
                  120);
      if (!pressure_info.converged)
        throw std::runtime_error(
          "pressure mass PCG failed inside Stokes preconditioner");

      for (int i = 0; i < n_pressure_dofs; ++i)
        dst[pressure_index(i, n_velocity_dofs)] =
          i == 0 ? pressure_solution[i] : -pressure_solution[i];

      return dst;
    }

  private:
    const RowSparseMatrix &system_matrix;
    const RowSparseMatrix &velocity_matrix;
    const RowSparseMatrix &pressure_mass_matrix;
    AMGPreconditioner *velocity_amg = nullptr;
    std::vector<double> velocity_inverse_diagonal;
    std::vector<double> pressure_inverse_diagonal;
    int n_velocity_dofs = 0;
    int n_pressure_dofs = 0;
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

  StokesResult solve_stokes(const std::string &mesh_file,
                            const double mesh_h,
                            const PreconditionerKind preconditioner_kind,
                            const double gmres_tolerance,
                            const int max_gmres_iterations)
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
    RowSparseMatrix velocity_matrix(n_velocity_dofs);
    RowSparseMatrix pressure_mass_matrix(n_pressure_dofs);
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
                    velocity_matrix.add(vi, vj, stiffness);
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
                  {
                    const int pj_dof = pressure_element_dofs[j];
                    const double mass_entry =
                      jxw * pressure_values[i][q] * pressure_values[j][q];
                    pressure_mass_matrix.add(pi_dof, pj_dof, mass_entry);
                    if (pi_dof == pj_dof)
                      pressure_mass_diag[pi_dof] += mass_entry;
                  }
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (velocity_space.dofBoundaryMark(i) != 0)
        {
          for (int component = 0; component < dim; ++component)
            apply_homogeneous_constraint(matrix,
                                         rhs,
                                         velocity_index(component,
                                                        i,
                                                        n_velocity_dofs));
          apply_homogeneous_constraint(velocity_matrix, i);
        }

    apply_homogeneous_constraint(matrix,
                                 rhs,
                                 pressure_index(0, n_velocity_dofs));
    apply_homogeneous_constraint(pressure_mass_matrix, 0);
    pressure_mass_diag[0] = 1.0;

    std::vector<double> solution(n_total, 0.0);
    SolveInfo gmres_info;
    switch (preconditioner_kind)
      {
        case PreconditionerKind::block_diagonal:
          {
            BlockDiagonalPreconditioner preconditioner(matrix,
                                                       pressure_mass_diag,
                                                       n_velocity_dofs);
            gmres_info =
              solve_right_preconditioned_gmres(matrix,
                                               solution,
                                               rhs,
                                               preconditioner,
                                               gmres_tolerance,
                                               100,
                                               max_gmres_iterations);
            break;
          }
        case PreconditionerKind::schur_mass_pcg:
          {
            SchurMassPreconditioner preconditioner(matrix,
                                                   velocity_matrix,
                                                   pressure_mass_matrix,
                                                   nullptr,
                                                   n_velocity_dofs,
                                                   n_pressure_dofs);
            gmres_info =
              solve_right_preconditioned_gmres(matrix,
                                               solution,
                                               rhs,
                                               preconditioner,
                                               gmres_tolerance,
                                               100,
                                               max_gmres_iterations);
            break;
          }
        case PreconditionerKind::schur_mass_amg1:
        case PreconditionerKind::schur_mass_amg2:
        case PreconditionerKind::schur_mass_amg3:
          {
            StiffMatrix<dim, double> velocity_stiffness(velocity_space);
            velocity_stiffness.algebricAccuracy() = 4;
            velocity_stiffness.build();
            FEMFunction<double, dim> velocity_boundary_solution(
              velocity_space);
            Vector<double> velocity_zero_rhs(n_velocity_dofs);
            velocity_zero_rhs = 0.0;
            apply_homogeneous_dirichlet_boundary(
              velocity_stiffness,
              velocity_boundary_solution,
              velocity_zero_rhs,
              velocity_space);

            AMGPreconditioner velocity_amg(velocity_stiffness, 2, 2);
            velocity_amg.smoothStep() =
              velocity_amg_smooth_steps(preconditioner_kind);
            SchurMassPreconditioner preconditioner(matrix,
                                                   velocity_matrix,
                                                   pressure_mass_matrix,
                                                   &velocity_amg,
                                                   n_velocity_dofs,
                                                   n_pressure_dofs);
            gmres_info =
              solve_right_preconditioned_gmres(matrix,
                                               solution,
                                               rhs,
                                               preconditioner,
                                               gmres_tolerance,
                                               100,
                                               max_gmres_iterations);
            break;
          }
      }
    if (!gmres_info.converged)
      throw std::runtime_error(
        "right-preconditioned GMRES failed after " +
        std::to_string(gmres_info.iterations) +
        " iterations; final relative residual=" +
        std::to_string(gmres_info.relative_residual));

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

    return {preconditioner_name(preconditioner_kind),
            mesh_h,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            matrix.nonzeros(),
            gmres_info.iterations,
            gmres_info.relative_residual,
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            Functional::L2Error(p_h, exact_p, 5),
            divergence_l2_norm(u_x, u_y)};
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-directory]\n";
      return 2;
    }

  const std::string mesh_dir =
    argc > 1 ? argv[1] : std::string(AFEPACK_DEFAULT_MESH_DIR);
  const std::vector<MeshCase> mesh_cases = {
    {0.20, mesh_dir + "/unit_square_h0p20"},
    {0.10, mesh_dir + "/unit_square_h0p10"},
    {0.05, mesh_dir + "/unit_square_h0p05"},
    {0.025, mesh_dir + "/unit_square_h0p025"}};

  const double gmres_tolerance = 1.0e-6;
  const int max_gmres_iterations = 5000;

  try
    {
      std::vector<StokesResult> results;
      for (const auto &mesh_case : mesh_cases)
        for (const auto preconditioner_kind :
             preconditioners_for_mesh(mesh_case))
          results.push_back(solve_stokes(mesh_case.mesh_file,
                                         mesh_case.h,
                                         preconditioner_kind,
                                         gmres_tolerance,
                                         max_gmres_iterations));
      add_convergence_rates(results);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes sparse block GMRES mesh sweep "
                   "Taylor-Hood manufactured solution\n";
      std::cout << "GMRES tolerance: " << gmres_tolerance << '\n';
      std::cout << "Finest h=0.025 mesh runs AMG velocity-block Schur "
                   "variants only; block diagonal and Jacobi-PCG velocity "
                   "blocks are omitted there.\n";
      std::cout << "preconditioner h cells velocity_dofs pressure_dofs "
                   "nonzeros GMRES linear_residual velocity_L2 velocity_H1 "
                   "pressure_L2 divergence_L2 velocity_L2_rate "
                   "velocity_H1_rate pressure_L2_rate divergence_rate\n";
      for (const auto &result : results)
        std::cout << std::setw(17) << result.preconditioner << ' '
                  << result.mesh_h << ' '
                  << result.cells << ' '
                  << result.velocity_dofs << ' '
                  << result.pressure_dofs << ' '
                  << result.system_nonzeros << ' '
                  << result.gmres_iterations << ' '
                  << result.gmres_relative_residual << ' '
                  << result.velocity_l2_error << ' '
                  << result.velocity_h1_error << ' '
                  << result.pressure_l2_error << ' '
                  << result.divergence_l2_norm << ' '
                  << std::setw(9)
                  << format_rate(result.has_rates,
                                 result.velocity_l2_rate) << ' '
                  << std::setw(9)
                  << format_rate(result.has_rates,
                                 result.velocity_h1_rate) << ' '
                  << std::setw(9)
                  << format_rate(result.has_rates,
                                 result.pressure_l2_rate) << ' '
                  << std::setw(9)
                  << format_rate(result.has_rates,
                                 result.divergence_rate) << '\n';
      print_gmres_summary(results);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
