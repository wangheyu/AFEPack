/**
 * @file oseen_block_pcd_iterative_afepack.cpp
 * @brief AFEPack Oseen 线性化方程算例：迭代线性求解。
 *
 * @details
 * 本文件是 FAST 目录中的 Oseen 线性化方程 迁移算例，关注Navier-Stokes 线性化后得到的鞍点系统及其 Schur 补预条件器。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 迭代线性求解：使用 Krylov 迭代和块预条件器替代直接求解，以观察可扩展性。
 * - PCD 预条件器：构造压力对流扩散型 Schur 补近似，评估 PCD 预条件的迭代效率。
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
  double           viscosity = 1.0;
  constexpr double pressure_laplace_regularization = 1.0e-6;

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct OseenResult
  {
    std::string  preconditioner;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    int          gmres_iterations = 0;
    bool         gmres_converged = false;
    double       gmres_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    double       pressure_l2_error = 0.0;
    double       divergence_l2_norm = 0.0;
  };

  enum class PreconditionerKind
  {
    pressure_mass,
    pressure_mass_amg,
    pcd,
    pcd_amg,
    pcd_oseen_amg
  };

  class JacobiBiCGSTABSolver
  {
  public:
    JacobiBiCGSTABSolver() = default;

    JacobiBiCGSTABSolver(const std::vector<double> &matrix, const int size)
    {
      reinit(matrix, size);
    }

    void reinit(const std::vector<double> &matrix, const int size)
    {
      n = size;
      entries = matrix;
      inverse_diagonal.assign(n, 1.0);
      for (int i = 0; i < n; ++i)
        {
          const double diagonal = entries[static_cast<std::size_t>(i) * n + i];
          inverse_diagonal[i] =
            std::abs(diagonal) > 1.0e-14 ? 1.0 / diagonal : 1.0;
        }
    }

    std::vector<double> solve(const std::vector<double> &rhs) const
    {
      if (static_cast<int>(rhs.size()) != n)
        throw std::runtime_error("JacobiBiCGSTABSolver solve size mismatch");

      std::vector<double> x(n, 0.0);
      std::vector<double> r = rhs;
      const std::vector<double> r_hat = r;
      const double rhs_norm = std::max(vector_norm_local(rhs), 1.0e-30);
      if (vector_norm_local(r) / rhs_norm <= tolerance)
        return x;

      std::vector<double> p(n, 0.0);
      std::vector<double> v(n, 0.0);
      double rho_old = 1.0;
      double alpha = 1.0;
      double omega = 1.0;

      for (int iteration = 0; iteration < max_iterations; ++iteration)
        {
          const double rho_new = dot_local(r_hat, r);
          if (std::abs(rho_new) < 1.0e-300)
            break;

          const double beta = (rho_new / rho_old) * (alpha / omega);
          for (int i = 0; i < n; ++i)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);

          const auto p_hat = apply_jacobi(p);
          v = vmult(p_hat);
          const double denominator = dot_local(r_hat, v);
          if (std::abs(denominator) < 1.0e-300)
            break;
          alpha = rho_new / denominator;

          std::vector<double> s = r;
          for (int i = 0; i < n; ++i)
            s[i] -= alpha * v[i];

          if (vector_norm_local(s) / rhs_norm <= tolerance)
            {
              for (int i = 0; i < n; ++i)
                x[i] += alpha * p_hat[i];
              return x;
            }

          const auto s_hat = apply_jacobi(s);
          const auto t = vmult(s_hat);
          const double tt = dot_local(t, t);
          if (std::abs(tt) < 1.0e-300)
            break;
          omega = dot_local(t, s) / tt;

          for (int i = 0; i < n; ++i)
            {
              x[i] += alpha * p_hat[i] + omega * s_hat[i];
              r[i] = s[i] - omega * t[i];
            }

          if (vector_norm_local(r) / rhs_norm <= tolerance)
            return x;

          if (std::abs(omega) < 1.0e-300)
            break;
          rho_old = rho_new;
        }

      return x;
    }

  private:
    static constexpr int    max_iterations = 300;
    static constexpr double tolerance = 1.0e-7;

    std::vector<double> vmult(const std::vector<double> &x) const
    {
      std::vector<double> y(n, 0.0);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          y[i] += entries[static_cast<std::size_t>(i) * n + j] * x[j];
      return y;
    }

    std::vector<double> apply_jacobi(const std::vector<double> &x) const
    {
      std::vector<double> y = x;
      for (int i = 0; i < n; ++i)
        y[i] *= inverse_diagonal[i];
      return y;
    }

    static double dot_local(const std::vector<double> &a,
                            const std::vector<double> &b)
    {
      double result = 0.0;
      for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
      return result;
    }

    static double vector_norm_local(const std::vector<double> &v)
    {
      return std::sqrt(dot_local(v, v));
    }

    int n = 0;
    std::vector<double> entries;
    std::vector<double> inverse_diagonal;
  };

  class JacobiPCGSolver
  {
  public:
    JacobiPCGSolver() = default;

    JacobiPCGSolver(const std::vector<double> &matrix, const int size)
    {
      reinit(matrix, size);
    }

    void reinit(const std::vector<double> &matrix, const int size)
    {
      n = size;
      entries = matrix;
      inverse_diagonal.assign(n, 1.0);
      for (int i = 0; i < n; ++i)
        {
          const double diagonal = entries[static_cast<std::size_t>(i) * n + i];
          inverse_diagonal[i] =
            std::abs(diagonal) > 1.0e-14 ? 1.0 / diagonal : 1.0;
        }
    }

    std::vector<double> solve(const std::vector<double> &rhs) const
    {
      if (static_cast<int>(rhs.size()) != n)
        throw std::runtime_error("JacobiPCGSolver solve size mismatch");

      std::vector<double> x(n, 0.0);
      std::vector<double> r = rhs;
      const double rhs_norm = std::max(vector_norm_local(rhs), 1.0e-30);
      if (vector_norm_local(r) / rhs_norm <= tolerance)
        return x;

      std::vector<double> z = apply_jacobi(r);
      std::vector<double> p = z;
      double rho = dot_local(r, z);

      for (int iteration = 0; iteration < max_iterations; ++iteration)
        {
          const std::vector<double> ap = vmult(p);
          const double denominator = dot_local(p, ap);
          if (std::abs(denominator) < 1.0e-30)
            break;

          const double alpha = rho / denominator;
          for (int i = 0; i < n; ++i)
            {
              x[i] += alpha * p[i];
              r[i] -= alpha * ap[i];
            }

          if (vector_norm_local(r) / rhs_norm <= tolerance)
            return x;

          z = apply_jacobi(r);
          const double rho_new = dot_local(r, z);
          const double beta = rho_new / std::max(rho, 1.0e-300);
          for (int i = 0; i < n; ++i)
            p[i] = z[i] + beta * p[i];
          rho = rho_new;
        }

      return x;
    }

  private:
    static constexpr int    max_iterations = 200;
    static constexpr double tolerance = 1.0e-8;

    std::vector<double> vmult(const std::vector<double> &x) const
    {
      std::vector<double> y(n, 0.0);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          y[i] += entries[static_cast<std::size_t>(i) * n + j] * x[j];
      return y;
    }

    std::vector<double> apply_jacobi(const std::vector<double> &x) const
    {
      std::vector<double> y = x;
      for (int i = 0; i < n; ++i)
        y[i] *= inverse_diagonal[i];
      return y;
    }

    static double dot_local(const std::vector<double> &a,
                            const std::vector<double> &b)
    {
      double result = 0.0;
      for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
      return result;
    }

    static double vector_norm_local(const std::vector<double> &v)
    {
      return std::sqrt(dot_local(v, v));
    }

    int n = 0;
    std::vector<double> entries;
    std::vector<double> inverse_diagonal;
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

  double rhs_x(const double *p)
  {
    return -viscosity * laplace_velocity_x(p) + oseen_convection_x(p) +
           pressure_gradient(p)[0];
  }

  double rhs_y(const double *p)
  {
    return -viscosity * laplace_velocity_y(p) + oseen_convection_y(p) +
           pressure_gradient(p)[1];
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

  void subtract_arithmetic_mean(std::vector<double> &v)
  {
    double sum = 0.0;
    for (const double value : v)
      sum += value;

    const double mean = sum / static_cast<double>(v.size());
    for (double &value : v)
      value -= mean;
  }

  std::vector<double> apply_amg_preconditioner(AMGPreconditioner &amg,
                                               const std::vector<double> &rhs)
  {
    Vector<double> amg_rhs(rhs.size());
    Vector<double> amg_solution(rhs.size());
    for (std::size_t i = 0; i < rhs.size(); ++i)
      amg_rhs(i) = rhs[i];

    amg_solution = 0.0;
    amg.vmult(amg_solution, amg_rhs);

    std::vector<double> solution(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i)
      solution[i] = amg_solution(i);
    return solution;
  }

  std::vector<double> solve_amg_preconditioned_bicgstab(
    const std::vector<double> &matrix,
    const int n,
    const std::vector<double> &rhs,
    AMGPreconditioner &amg)
  {
    constexpr int    max_iterations = 80;
    constexpr double tolerance = 1.0e-7;

    std::vector<double> x(n, 0.0);
    std::vector<double> r = rhs;
    const std::vector<double> r_hat = r;
    const double rhs_norm = std::max(vector_norm(rhs), 1.0e-30);
    if (vector_norm(r) / rhs_norm <= tolerance)
      return x;

    std::vector<double> p(n, 0.0);
    std::vector<double> v(n, 0.0);
    double rho_old = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    for (int iteration = 0; iteration < max_iterations; ++iteration)
      {
        const double rho_new = vector_dot(r_hat, r);
        if (std::abs(rho_new) < 1.0e-300)
          break;

        const double beta = (rho_new / rho_old) * (alpha / omega);
        for (int i = 0; i < n; ++i)
          p[i] = r[i] + beta * (p[i] - omega * v[i]);

        const auto p_hat = apply_amg_preconditioner(amg, p);
        v = dense_vmult(matrix, p_hat, n);
        const double denominator = vector_dot(r_hat, v);
        if (std::abs(denominator) < 1.0e-300)
          break;
        alpha = rho_new / denominator;

        std::vector<double> s = r;
        for (int i = 0; i < n; ++i)
          s[i] -= alpha * v[i];

        if (vector_norm(s) / rhs_norm <= tolerance)
          {
            for (int i = 0; i < n; ++i)
              x[i] += alpha * p_hat[i];
            return x;
          }

        const auto s_hat = apply_amg_preconditioner(amg, s);
        const auto t = dense_vmult(matrix, s_hat, n);
        const double tt = vector_dot(t, t);
        if (std::abs(tt) < 1.0e-300)
          break;
        omega = vector_dot(t, s) / tt;

        for (int i = 0; i < n; ++i)
          {
            x[i] += alpha * p_hat[i] + omega * s_hat[i];
            r[i] = s[i] - omega * t[i];
          }

        if (vector_norm(r) / rhs_norm <= tolerance)
          return x;

        if (std::abs(omega) < 1.0e-300)
          break;
        rho_old = rho_new;
      }

    return x;
  }

  std::string preconditioner_name(const PreconditionerKind kind)
  {
    switch (kind)
      {
        case PreconditionerKind::pressure_mass:
          return "diag(Fbicg,Mp-pcg)";
        case PreconditionerKind::pressure_mass_amg:
          return "diag(Famg,Mp-pcg)";
        case PreconditionerKind::pcd:
          return "diag(Fbicg,PCD)";
        case PreconditionerKind::pcd_amg:
          return "diag(Famg,PCD)";
        case PreconditionerKind::pcd_oseen_amg:
          return "diag(Famg-BiCG,PCD)";
      }
    return "unknown";
  }

  std::string preconditioner_output_suffix(const PreconditionerKind kind)
  {
    switch (kind)
      {
        case PreconditionerKind::pressure_mass:
          return "mass";
        case PreconditionerKind::pressure_mass_amg:
          return "mass_amg";
        case PreconditionerKind::pcd:
          return "pcd";
        case PreconditionerKind::pcd_amg:
          return "pcd_amg";
        case PreconditionerKind::pcd_oseen_amg:
          return "pcd_oseen_amg";
      }
    return "unknown";
  }

  int gmres_max_iterations(const PreconditionerKind kind)
  {
    switch (kind)
      {
        case PreconditionerKind::pressure_mass:
        case PreconditionerKind::pressure_mass_amg:
          return 1000;
        case PreconditionerKind::pcd:
        case PreconditionerKind::pcd_amg:
        case PreconditionerKind::pcd_oseen_amg:
          return 5000;
      }
    return 1000;
  }

  class BlockMassPreconditioner
  {
  public:
    BlockMassPreconditioner(const int n_velocity_dofs,
                            const int n_pressure_dofs,
                            const std::vector<double> &velocity_block,
                            const std::vector<double> &pressure_mass,
                            AMGPreconditioner *velocity_amg = nullptr)
      : n_velocity_dofs(n_velocity_dofs),
        n_pressure_dofs(n_pressure_dofs),
        velocity_solver(velocity_block, n_velocity_dofs),
        mass_solver(pressure_mass, n_pressure_dofs),
        velocity_amg(velocity_amg)
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

          std::vector<double> component_solution(n_velocity_dofs, 0.0);
          if (velocity_amg != nullptr)
            {
              Vector<double> amg_src(n_velocity_dofs);
              Vector<double> amg_dst(n_velocity_dofs);
              for (int i = 0; i < n_velocity_dofs; ++i)
                amg_src(i) = component_rhs[i];
              amg_dst = 0.0;
              velocity_amg->vmult(amg_dst, amg_src);
              for (int i = 0; i < n_velocity_dofs; ++i)
                component_solution[i] = amg_dst(i);
            }
          else
            {
              component_solution = velocity_solver.solve(component_rhs);
            }
          for (int i = 0; i < n_velocity_dofs; ++i)
            dst[velocity_index(component, i, n_velocity_dofs)] =
              component_solution[i];
        }

      std::vector<double> pressure_rhs(n_pressure_dofs);
      for (int i = 0; i < n_pressure_dofs; ++i)
        pressure_rhs[i] =
          src[pressure_index(i, n_velocity_dofs)];
      subtract_arithmetic_mean(pressure_rhs);

      const auto pressure_solution = mass_solver.solve(pressure_rhs);
      for (int i = 0; i < n_pressure_dofs; ++i)
        dst[pressure_index(i, n_velocity_dofs)] = pressure_solution[i];

      return dst;
    }

  private:
    int n_velocity_dofs = 0;
    int n_pressure_dofs = 0;
    JacobiBiCGSTABSolver velocity_solver;
    JacobiPCGSolver mass_solver;
    AMGPreconditioner *velocity_amg = nullptr;
  };

  class BlockPCDPreconditioner
  {
  public:
    BlockPCDPreconditioner(const int n_velocity_dofs,
                           const int n_pressure_dofs,
                           const std::vector<double> &velocity_block,
                           const std::vector<double> &pressure_mass,
                           const std::vector<double> &pressure_laplace,
                           const std::vector<double> &pressure_convdiff,
                           AMGPreconditioner *velocity_amg = nullptr,
                           const bool use_full_velocity_solve = false)
      : n_velocity_dofs(n_velocity_dofs),
        n_pressure_dofs(n_pressure_dofs),
        velocity_block(velocity_block),
        velocity_solver(velocity_block, n_velocity_dofs),
        mass_solver(pressure_mass, n_pressure_dofs),
        laplace_solver(pressure_laplace, n_pressure_dofs),
        pressure_convdiff(pressure_convdiff),
        velocity_amg(velocity_amg),
        use_full_velocity_solve(use_full_velocity_solve)
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

          std::vector<double> component_solution(n_velocity_dofs, 0.0);
          if (velocity_amg != nullptr && use_full_velocity_solve)
            {
              component_solution =
                solve_amg_preconditioned_bicgstab(velocity_block,
                                                  n_velocity_dofs,
                                                  component_rhs,
                                                  *velocity_amg);
            }
          else if (velocity_amg != nullptr)
            {
              Vector<double> amg_src(n_velocity_dofs);
              Vector<double> amg_dst(n_velocity_dofs);
              for (int i = 0; i < n_velocity_dofs; ++i)
                amg_src(i) = component_rhs[i];
              amg_dst = 0.0;
              velocity_amg->vmult(amg_dst, amg_src);
              for (int i = 0; i < n_velocity_dofs; ++i)
                component_solution[i] = amg_dst(i);
            }
          else
            {
              component_solution = velocity_solver.solve(component_rhs);
            }
          for (int i = 0; i < n_velocity_dofs; ++i)
            dst[velocity_index(component, i, n_velocity_dofs)] =
              component_solution[i];
        }

      std::vector<double> pressure_rhs(n_pressure_dofs);
      for (int i = 0; i < n_pressure_dofs; ++i)
        pressure_rhs[i] =
          src[pressure_index(i, n_velocity_dofs)];
      subtract_arithmetic_mean(pressure_rhs);

      auto pressure_tmp = mass_solver.solve(pressure_rhs);
      pressure_tmp =
        dense_vmult(pressure_convdiff, pressure_tmp, n_pressure_dofs);
      subtract_arithmetic_mean(pressure_tmp);
      auto pressure_solution = laplace_solver.solve(pressure_tmp);

      for (int i = 0; i < n_pressure_dofs; ++i)
        dst[pressure_index(i, n_velocity_dofs)] = pressure_solution[i];

      return dst;
    }

  private:
    int n_velocity_dofs = 0;
    int n_pressure_dofs = 0;
    std::vector<double> velocity_block;
    JacobiBiCGSTABSolver velocity_solver;
    JacobiPCGSolver mass_solver;
    JacobiPCGSolver laplace_solver;
    std::vector<double> pressure_convdiff;
    AMGPreconditioner *velocity_amg = nullptr;
    bool use_full_velocity_solve = false;
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
    std::vector<double> pressure_mass;
    std::vector<double> pressure_laplace;
    std::vector<double> pressure_convdiff;
  };

  AssembledOseenSystem assemble_oseen_system(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    AssembledOseenSystem system;
    system.matrix.assign(static_cast<std::size_t>(n_total) * n_total, 0.0);
    system.rhs.assign(n_total, 0.0);
    system.velocity_block.assign(
      static_cast<std::size_t>(n_velocity_dofs) * n_velocity_dofs, 0.0);
    system.pressure_mass.assign(
      static_cast<std::size_t>(n_pressure_dofs) * n_pressure_dofs, 0.0);
    system.pressure_laplace.assign(
      static_cast<std::size_t>(n_pressure_dofs) * n_pressure_dofs, 0.0);
    system.pressure_convdiff.assign(
      static_cast<std::size_t>(n_pressure_dofs) * n_pressure_dofs, 0.0);

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
        const auto pressure_grads =
          pressure_element.basis_function_gradient(q_point);

        const auto &velocity_element_dofs = velocity_element.dof();
        const auto &pressure_element_dofs = pressure_element.dof();
        const int n_local_velocity_dofs = velocity_element_dofs.size();
        const int n_local_pressure_dofs = pressure_element_dofs.size();

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double p[dim] = {q_point[q][0], q_point[q][1]};
            const double wind[dim] = {velocity_x_exact(p),
                                      velocity_y_exact(p)};
            const double f[dim] = {rhs_x(p), rhs_y(p)};

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
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
                      add_entry(system.matrix,
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
                    const double mass =
                      pressure_values[j][q] * pressure_values[i][q];
                    const double laplace =
                      innerProduct(pressure_grads[i][q],
                                   pressure_grads[j][q]);
                    const double convection =
                      (wind[0] * pressure_grads[j][q][0] +
                       wind[1] * pressure_grads[j][q][1]) *
                      pressure_values[i][q];

                    add_entry(system.pressure_mass,
                              n_pressure_dofs,
                              pi_dof,
                              pj,
                              jxw * mass);
                    add_entry(system.pressure_laplace,
                              n_pressure_dofs,
                              pi_dof,
                              pj,
                              jxw * (laplace +
                                     pressure_laplace_regularization * mass));
                    add_entry(system.pressure_convdiff,
                              n_pressure_dofs,
                              pi_dof,
                              pj,
                              jxw * (viscosity * laplace + convection));
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

  template <typename Preconditioner>
  OseenResult solve_oseen_with_preconditioner(
    const EasyMesh &mesh,
    FEMSpace<double, dim> &velocity_space,
    FEMSpace<double, dim> &pressure_space,
    const AssembledOseenSystem &system,
    const Preconditioner &preconditioner,
    const std::string &preconditioner_label,
    const int max_gmres_iterations,
    const std::string &output_prefix)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    std::vector<double> solution(n_total, 0.0);
    const SolveInfo gmres_info =
      solve_right_preconditioned_gmres(system.matrix,
                                       solution,
                                       system.rhs,
                                       n_total,
                                       preconditioner,
                                       1.0e-8,
                                       80,
                                       max_gmres_iterations);

    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);
    assign_solution(solution, u_x, u_y, p_h);

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

    return {preconditioner_label,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            gmres_info.iterations,
            gmres_info.converged,
            gmres_info.relative_residual,
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            Functional::L2Error(p_h, exact_p, 5),
            divergence_l2_norm(u_x, u_y)};
  }

  std::vector<OseenResult> solve_oseen_block_sweep(
    const std::string &mesh_file,
    const std::string &output_prefix,
    const std::vector<PreconditionerKind> &preconditioner_kinds,
    const bool append_preconditioner_suffix = true)
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

    auto system = assemble_oseen_system(velocity_space, pressure_space);

    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    std::vector<OseenResult> results;
    results.reserve(preconditioner_kinds.size());

    for (const PreconditionerKind kind : preconditioner_kinds)
      {
        const std::string label = preconditioner_name(kind);
        const std::string case_output_prefix =
          append_preconditioner_suffix ?
            output_prefix + "_" + preconditioner_output_suffix(kind) :
            output_prefix;
        const int max_gmres_iterations = gmres_max_iterations(kind);
        switch (kind)
          {
            case PreconditionerKind::pressure_mass:
              {
                BlockMassPreconditioner preconditioner(
                  n_velocity_dofs,
                  n_pressure_dofs,
                  system.velocity_block,
                  system.pressure_mass);
                results.push_back(
                  solve_oseen_with_preconditioner(mesh,
                                                  velocity_space,
                                                  pressure_space,
                                                  system,
                                                  preconditioner,
                                                  label,
                                                  max_gmres_iterations,
                                                  case_output_prefix));
                break;
              }
            case PreconditionerKind::pressure_mass_amg:
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
                velocity_amg.smoothStep() = 2;
                BlockMassPreconditioner preconditioner(
                  n_velocity_dofs,
                  n_pressure_dofs,
                  system.velocity_block,
                  system.pressure_mass,
                  &velocity_amg);
                results.push_back(
                  solve_oseen_with_preconditioner(mesh,
                                                  velocity_space,
                                                  pressure_space,
                                                  system,
                                                  preconditioner,
                                                  label,
                                                  max_gmres_iterations,
                                                  case_output_prefix));
                break;
              }
            case PreconditionerKind::pcd:
              {
                BlockPCDPreconditioner preconditioner(
                  n_velocity_dofs,
                  n_pressure_dofs,
                  system.velocity_block,
                  system.pressure_mass,
                  system.pressure_laplace,
                  system.pressure_convdiff);
                results.push_back(
                  solve_oseen_with_preconditioner(mesh,
                                                  velocity_space,
                                                  pressure_space,
                                                  system,
                                                  preconditioner,
                                                  label,
                                                  max_gmres_iterations,
                                                  case_output_prefix));
                break;
              }
            case PreconditionerKind::pcd_amg:
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
                velocity_amg.smoothStep() = 2;
                BlockPCDPreconditioner preconditioner(
                  n_velocity_dofs,
                  n_pressure_dofs,
                  system.velocity_block,
                  system.pressure_mass,
                  system.pressure_laplace,
                  system.pressure_convdiff,
                  &velocity_amg);
                results.push_back(
                  solve_oseen_with_preconditioner(mesh,
                                                  velocity_space,
                                                  pressure_space,
                                                  system,
                                                  preconditioner,
                                                  label,
                                                  max_gmres_iterations,
                                                  case_output_prefix));
                break;
              }
            case PreconditionerKind::pcd_oseen_amg:
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
                velocity_amg.smoothStep() = 2;
                BlockPCDPreconditioner preconditioner(
                  n_velocity_dofs,
                  n_pressure_dofs,
                  system.velocity_block,
                  system.pressure_mass,
                  system.pressure_laplace,
                  system.pressure_convdiff,
                  &velocity_amg,
                  true);
                results.push_back(
                  solve_oseen_with_preconditioner(mesh,
                                                  velocity_space,
                                                  pressure_space,
                                                  system,
                                                  preconditioner,
                                                  label,
                                                  max_gmres_iterations,
                                                  case_output_prefix));
                break;
              }
          }
      }

    return results;
  }

  OseenResult solve_oseen_block_pcd(const std::string &mesh_file,
                                    const std::string &output_prefix)
  {
    auto results = solve_oseen_block_sweep(mesh_file,
                                           output_prefix,
                                           {PreconditionerKind::pcd},
                                           false);
    return results.front();
  }

  struct MeshSweepCase
  {
    const char *label;
    const char *mesh_name;
    const char *output_suffix;
  };

  void print_single_result(const OseenResult &result,
                           const std::string &output_prefix)
  {
    std::cout << "cells: " << result.cells << '\n';
    std::cout << "preconditioner: " << result.preconditioner << '\n';
    std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
    std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
    std::cout << "right-preconditioned GMRES iterations: "
              << result.gmres_iterations << '\n';
    std::cout << "GMRES status: "
              << (result.gmres_converged ? "ok" : "limit") << '\n';
    std::cout << "GMRES relative residual: "
              << result.gmres_relative_residual << '\n';
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

  void print_sweep_row(const MeshSweepCase &mesh_case,
                       const OseenResult &result)
  {
    std::cout << result.preconditioner << ' '
              << mesh_case.label << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.gmres_iterations << ' '
              << (result.gmres_converged ? "ok" : "limit") << ' '
              << result.gmres_relative_residual << ' '
              << result.velocity_l2_error << ' '
              << result.velocity_h1_error << ' '
              << result.pressure_l2_error << ' '
              << result.divergence_l2_norm << '\n';
  }

  void run_default_mesh_sweep()
  {
    const std::vector<MeshSweepCase> mesh_cases = {
      {"2.000000e-01", "unit_square_h0p20", "h0p20"},
      {"1.000000e-01", "unit_square_h0p10", "h0p10"}};

    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    const std::string output_stem =
      output_dir + "/oseen_block_pcd_iterative_afepack";

    std::cout << "AFEPack Oseen block PCD with iterative inner solvers "
                 "Taylor-Hood manufactured solution mesh sweep\n";
    std::cout << "preconditioner h cells velocity_dofs pressure_dofs GMRES "
                 "status relative_residual velocity_L2 velocity_H1 "
                 "pressure_L2 divergence_L2\n";

    for (const MeshSweepCase &mesh_case : mesh_cases)
      {
        const std::string mesh_file =
          mesh_dir + "/" + std::string(mesh_case.mesh_name);
        const std::string output_prefix =
          output_stem + "_" + mesh_case.output_suffix;
        const auto results =
          solve_oseen_block_sweep(mesh_file,
                                  output_prefix,
                                  {PreconditionerKind::pressure_mass,
                                   PreconditionerKind::pressure_mass_amg,
                                   PreconditionerKind::pcd,
                                   PreconditionerKind::pcd_amg,
                                   PreconditionerKind::pcd_oseen_amg});
        for (const OseenResult &result : results)
          print_sweep_row(mesh_case, result);
      }
  }
}

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output-prefix]]\n"
                << "With no arguments, runs the built-in h=0.20,0.10 "
                   "mesh sweep.\n";
      return 2;
    }

  try
    {
      std::cout << std::scientific << std::setprecision(6);
      if (argc == 1)
        run_default_mesh_sweep();
      else
        {
          const std::string mesh_file = argv[1];
          const std::string output_prefix =
            argc > 2 ? argv[2] :
                       std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                         "/oseen_block_pcd_iterative_afepack";
          const OseenResult result =
            solve_oseen_block_pcd(mesh_file, output_prefix);

          std::cout << "AFEPack Oseen block PCD with iterative inner solvers "
                       "Taylor-Hood manufactured solution\n";
          print_single_result(result, output_prefix);
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
