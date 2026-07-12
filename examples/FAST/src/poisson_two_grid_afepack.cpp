/**
 * @file poisson_two_grid_afepack.cpp
 * @brief AFEPack Poisson 方程算例：两重网格算法。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 两重网格算法：展示粗细网格组合求解策略在模型问题中的效果。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
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
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef AFEPACK_DEFAULT_MESH_DIR
#  define AFEPACK_DEFAULT_MESH_DIR "build/meshes/afepack"
#endif

namespace
{
  constexpr int    dim = 2;
  constexpr double pi  = 3.141592653589793238462643383279502884;

  struct Level
  {
    double mesh_h = 0.0;
    EasyMesh mesh;
    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis;
    std::vector<TemplateElement<double, dim, dim>> template_element;
    FEMSpace<double, dim> space;
    std::unique_ptr<StiffMatrix<dim, double>> matrix;
    Vector<double> rhs;
  };

  struct TransferEntry
  {
    int    coarse_index = 0;
    double weight = 0.0;
  };

  struct Transfer
  {
    std::vector<std::vector<TransferEntry>> rows;
    int coarse_size = 0;
  };

  struct SolveInfo
  {
    bool   converged = false;
    int    iterations = 0;
    double relative_residual = 0.0;
  };

  struct TwoGridResult
  {
    double       coarse_h = 0.0;
    double       fine_h = 0.0;
    unsigned int coarse_dofs = 0;
    unsigned int fine_dofs = 0;
    int          cycles = 0;
    int          smoothing_sweeps = 0;
    double       initial_residual = 0.0;
    double       jacobi_residual = 0.0;
    double       two_grid_residual = 0.0;
    double       two_grid_factor = 0.0;
    double       l2_error = 0.0;
  };

  struct MultigridResult
  {
    int          levels = 0;
    double       finest_h = 0.0;
    unsigned int finest_dofs = 0;
    int          cycles = 0;
    int          smoothing_sweeps = 0;
    double       initial_residual = 0.0;
    double       jacobi_residual = 0.0;
    double       multigrid_residual = 0.0;
    double       multigrid_factor = 0.0;
    double       l2_error = 0.0;
  };

  struct KrylovResult
  {
    int          levels = 0;
    double       finest_h = 0.0;
    unsigned int finest_dofs = 0;
    double       tolerance = 0.0;
    SolveInfo    cg;
    SolveInfo    jacobi_pcg;
    SolveInfo    multigrid_pcg;
    double       l2_error = 0.0;
  };

  double exact_solution(const double *p)
  {
    return std::sin(pi * p[0]) * std::sin(pi * p[1]);
  }

  double right_hand_side(const double *p)
  {
    return 2.0 * pi * pi * exact_solution(p);
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
  make_triangle_template(TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle.1.tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle.1.bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(geometry, dof, transform, basis);
    return template_element;
  }

  void build_space(EasyMesh &mesh,
                   std::vector<TemplateElement<double, dim, dim>>
                     &template_element,
                   FEMSpace<double, dim> &fem_space)
  {
    fem_space.reinit(mesh, template_element);

    const int n_element = mesh.n_geometry(dim);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
      fem_space.element(i).reinit(fem_space, i, 0);

    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();
  }

  void apply_dirichlet_boundary(StiffMatrix<dim, double> &stiff_matrix,
                                FEMFunction<double, dim> &solution,
                                Vector<double> &rhs,
                                const FEMSpace<double, dim> &space)
  {
    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &exact_solution);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.apply(stiff_matrix, solution, rhs);
  }

  std::unique_ptr<Level> build_level(const double h,
                                     const std::string &mesh_file)
  {
    auto level = std::make_unique<Level>();
    level->mesh_h = h;
    level->mesh.readData(mesh_file);

    level->template_element =
      make_triangle_template(level->geometry,
                             level->transform,
                             level->dof,
                             level->basis);
    build_space(level->mesh, level->template_element, level->space);

    level->matrix = std::make_unique<StiffMatrix<dim, double>>(level->space);
    level->matrix->algebricAccuracy() = 4;
    level->matrix->build();

    FEMFunction<double, dim> boundary_solution(level->space);
    Operator::L2Discretize(&right_hand_side, level->space, level->rhs, 4);
    apply_dirichlet_boundary(*level->matrix,
                             boundary_solution,
                             level->rhs,
                             level->space);
    return level;
  }

  std::vector<double> to_std_vector(const Vector<double> &vector)
  {
    std::vector<double> result(vector.size(), 0.0);
    for (std::size_t i = 0; i < vector.size(); ++i)
      result[i] = vector(i);
    return result;
  }

  double dot(const std::vector<double> &a, const std::vector<double> &b)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
      result += a[i] * b[i];
    return result;
  }

  double l2_norm(const std::vector<double> &v)
  {
    return std::sqrt(dot(v, v));
  }

  std::vector<double> vmult(const SparseMatrix<double> &matrix,
                            const std::vector<double> &x)
  {
    const auto &pattern = matrix.get_sparsity_pattern();
    const size_t *row_start = pattern.get_rowstart_indices();
    const size_t *columns = pattern.get_column_numbers();

    std::vector<double> y(matrix.m(), 0.0);
    for (std::size_t row = 0; row < matrix.m(); ++row)
      for (std::size_t index = row_start[row]; index < row_start[row + 1];
           ++index)
        y[row] += matrix.global_entry(index) * x[columns[index]];
    return y;
  }

  std::vector<double> residual(const SparseMatrix<double> &matrix,
                               const std::vector<double> &x,
                               const std::vector<double> &rhs)
  {
    std::vector<double> result = rhs;
    const auto ax = vmult(matrix, x);
    for (std::size_t i = 0; i < result.size(); ++i)
      result[i] -= ax[i];
    return result;
  }

  double relative_residual(const SparseMatrix<double> &matrix,
                           const std::vector<double> &x,
                           const std::vector<double> &rhs)
  {
    return l2_norm(residual(matrix, x, rhs)) /
           std::max(l2_norm(rhs), 1.0e-30);
  }

  std::vector<double> apply_jacobi(const SparseMatrix<double> &matrix,
                                   const std::vector<double> &src)
  {
    std::vector<double> dst(src.size(), 0.0);
    for (std::size_t i = 0; i < src.size(); ++i)
      {
        const double diagonal = matrix.diag_element(i);
        dst[i] = std::abs(diagonal) > 1.0e-14 ? src[i] / diagonal : src[i];
      }
    return dst;
  }

  SolveInfo solve_cg(const SparseMatrix<double> &matrix,
                     std::vector<double> &x,
                     const std::vector<double> &rhs,
                     const double tolerance,
                     const int max_iterations)
  {
    const double rhs_norm = std::max(l2_norm(rhs), 1.0e-30);
    std::vector<double> r = residual(matrix, x, rhs);
    double relative = l2_norm(r) / rhs_norm;
    if (relative <= tolerance)
      return {true, 0, relative};

    std::vector<double> z = apply_jacobi(matrix, r);
    std::vector<double> p = z;
    double rz_old = dot(r, z);
    if (std::abs(rz_old) < 1.0e-300)
      return {false, 0, relative};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        const auto ap = vmult(matrix, p);
        const double denominator = dot(p, ap);
        if (std::abs(denominator) < 1.0e-300)
          return {false, iteration - 1, relative};

        const double alpha = rz_old / denominator;
        for (std::size_t i = 0; i < x.size(); ++i)
          {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
          }

        relative = l2_norm(r) / rhs_norm;
        if (relative <= tolerance)
          return {true, iteration, relative};

        z = apply_jacobi(matrix, r);
        const double rz_new = dot(r, z);
        if (std::abs(rz_old) < 1.0e-300)
          return {false, iteration, relative};

        const double beta = rz_new / rz_old;
        for (std::size_t i = 0; i < p.size(); ++i)
          p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
      }

    return {false, max_iterations, relative};
  }

  template <typename Preconditioner>
  SolveInfo solve_pcg(const SparseMatrix<double> &matrix,
                      std::vector<double> &x,
                      const std::vector<double> &rhs,
                      const double tolerance,
                      const int max_iterations,
                      const Preconditioner &preconditioner)
  {
    const double rhs_norm = std::max(l2_norm(rhs), 1.0e-30);
    std::vector<double> r = residual(matrix, x, rhs);
    double relative = l2_norm(r) / rhs_norm;
    if (relative <= tolerance)
      return {true, 0, relative};

    std::vector<double> z = preconditioner(r);
    std::vector<double> p = z;
    double rz_old = dot(r, z);
    if (std::abs(rz_old) < 1.0e-300)
      return {false, 0, relative};

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
      {
        const auto ap = vmult(matrix, p);
        const double denominator = dot(p, ap);
        if (std::abs(denominator) < 1.0e-300)
          return {false, iteration - 1, relative};

        const double alpha = rz_old / denominator;
        for (std::size_t i = 0; i < x.size(); ++i)
          {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
          }

        relative = l2_norm(r) / rhs_norm;
        if (relative <= tolerance)
          return {true, iteration, relative};

        z = preconditioner(r);
        const double rz_new = dot(r, z);
        if (std::abs(rz_old) < 1.0e-300)
          return {false, iteration, relative};

        const double beta = rz_new / rz_old;
        for (std::size_t i = 0; i < p.size(); ++i)
          p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
      }

    return {false, max_iterations, relative};
  }

  void weighted_jacobi_smooth(const SparseMatrix<double> &matrix,
                              std::vector<double> &x,
                              const std::vector<double> &rhs,
                              const int sweeps,
                              const double omega)
  {
    for (int sweep = 0; sweep < sweeps; ++sweep)
      {
        const auto r = residual(matrix, x, rhs);
        for (std::size_t i = 0; i < x.size(); ++i)
          {
            const double diagonal = matrix.diag_element(i);
            if (std::abs(diagonal) > 1.0e-14)
              x[i] += omega * r[i] / diagonal;
          }
      }
  }

  bool barycentric_coordinates(const Point<dim> &p,
                               const Point<dim> &a,
                               const Point<dim> &b,
                               const Point<dim> &c,
                               double &lambda_a,
                               double &lambda_b,
                               double &lambda_c)
  {
    const double determinant =
      (b[1] - c[1]) * (a[0] - c[0]) +
      (c[0] - b[0]) * (a[1] - c[1]);
    if (std::abs(determinant) < 1.0e-30)
      return false;

    lambda_a = ((b[1] - c[1]) * (p[0] - c[0]) +
                (c[0] - b[0]) * (p[1] - c[1])) /
               determinant;
    lambda_b = ((c[1] - a[1]) * (p[0] - c[0]) +
                (a[0] - c[0]) * (p[1] - c[1])) /
               determinant;
    lambda_c = 1.0 - lambda_a - lambda_b;
    return true;
  }

  std::vector<TransferEntry>
  interpolation_row(const FEMSpace<double, dim> &coarse_space,
                    const Point<dim> &point)
  {
    constexpr double tolerance = 1.0e-9;
    const EasyMesh &mesh =
      static_cast<const EasyMesh &>(coarse_space.mesh());

    for (unsigned int cell = 0; cell < mesh.n_geometry(dim); ++cell)
      {
        const GeometryBM &triangle = mesh.geometry(dim, cell);
        if (triangle.n_vertex() != 3)
          continue;

        const int v0 = triangle.vertex(0);
        const int v1 = triangle.vertex(1);
        const int v2 = triangle.vertex(2);

        double l0 = 0.0;
        double l1 = 0.0;
        double l2 = 0.0;
        if (!barycentric_coordinates(point,
                                     mesh.point(v0),
                                     mesh.point(v1),
                                     mesh.point(v2),
                                     l0,
                                     l1,
                                     l2))
          continue;

        if (l0 >= -tolerance && l1 >= -tolerance && l2 >= -tolerance &&
            l0 <= 1.0 + tolerance && l1 <= 1.0 + tolerance &&
            l2 <= 1.0 + tolerance)
          {
            std::vector<double> weights = {std::max(0.0, l0),
                                           std::max(0.0, l1),
                                           std::max(0.0, l2)};
            const double sum = weights[0] + weights[1] + weights[2];
            if (sum <= 0.0)
              throw std::runtime_error("invalid interpolation weights");
            for (double &weight : weights)
              weight /= sum;

            const int vertices[3] = {v0, v1, v2};
            std::vector<TransferEntry> row;
            for (int i = 0; i < 3; ++i)
              {
                if (weights[i] <= 1.0e-14)
                  continue;
                const auto &dofs = coarse_space.geometryDof(0, vertices[i]);
                if (dofs.size() != 1)
                  throw std::runtime_error(
                    "expected one P1 dof on each coarse vertex");
                row.push_back({dofs[0], weights[i]});
              }
            return row;
          }
      }

    throw std::runtime_error("fine-grid dof is outside the coarse mesh");
  }

  Transfer build_transfer(const FEMSpace<double, dim> &fine_space,
                          const FEMSpace<double, dim> &coarse_space)
  {
    Transfer transfer;
    transfer.coarse_size = coarse_space.n_dof();
    transfer.rows.resize(fine_space.n_dof());

    for (int i = 0; i < fine_space.n_dof(); ++i)
      transfer.rows[i] =
        interpolation_row(coarse_space, fine_space.dofInfo(i).interp_point);

    return transfer;
  }

  std::vector<double> restrict_residual(const Transfer &transfer,
                                        const std::vector<double> &fine_vector)
  {
    std::vector<double> coarse_vector(transfer.coarse_size, 0.0);
    for (std::size_t i = 0; i < transfer.rows.size(); ++i)
      for (const auto &entry : transfer.rows[i])
        coarse_vector[entry.coarse_index] += entry.weight * fine_vector[i];
    return coarse_vector;
  }

  std::vector<double> prolongate(const Transfer &transfer,
                                 const std::vector<double> &coarse_vector)
  {
    std::vector<double> fine_vector(transfer.rows.size(), 0.0);
    for (std::size_t i = 0; i < transfer.rows.size(); ++i)
      for (const auto &entry : transfer.rows[i])
        fine_vector[i] += entry.weight * coarse_vector[entry.coarse_index];
    return fine_vector;
  }

  void two_grid_cycle(const Level &fine,
                      const Level &coarse,
                      const Transfer &transfer,
                      std::vector<double> &x,
                      const int pre_smoothing_steps,
                      const int post_smoothing_steps,
                      const double omega)
  {
    weighted_jacobi_smooth(*fine.matrix,
                           x,
                           to_std_vector(fine.rhs),
                           pre_smoothing_steps,
                           omega);

    const auto fine_residual =
      residual(*fine.matrix, x, to_std_vector(fine.rhs));
    const auto coarse_rhs = restrict_residual(transfer, fine_residual);

    std::vector<double> coarse_error(coarse.space.n_dof(), 0.0);
    const SolveInfo coarse_solve =
      solve_cg(*coarse.matrix, coarse_error, coarse_rhs, 1.0e-12, 5000);
    if (!coarse_solve.converged)
      throw std::runtime_error("coarse correction CG did not converge");

    const auto fine_error = prolongate(transfer, coarse_error);
    for (std::size_t i = 0; i < x.size(); ++i)
      x[i] += fine_error[i];

    weighted_jacobi_smooth(*fine.matrix,
                           x,
                           to_std_vector(fine.rhs),
                           post_smoothing_steps,
                           omega);
  }

  void multigrid_v_cycle(
    const std::vector<std::unique_ptr<Level>> &levels,
    const std::vector<Transfer> &transfers,
    const int level_index,
    std::vector<double> &x,
    const std::vector<double> &rhs,
    const int pre_smoothing_steps,
    const int post_smoothing_steps,
    const double omega)
  {
    const Level &level = *levels[level_index];

    if (level_index == 0)
      {
        const SolveInfo coarse_solve =
          solve_cg(*level.matrix, x, rhs, 1.0e-12, 5000);
        if (!coarse_solve.converged)
          throw std::runtime_error("coarsest multigrid CG solve failed");
        return;
      }

    weighted_jacobi_smooth(*level.matrix,
                           x,
                           rhs,
                           pre_smoothing_steps,
                           omega);

    const auto fine_residual = residual(*level.matrix, x, rhs);
    const auto coarse_rhs =
      restrict_residual(transfers[level_index - 1], fine_residual);

    std::vector<double> coarse_error(levels[level_index - 1]->space.n_dof(),
                                     0.0);
    multigrid_v_cycle(levels,
                      transfers,
                      level_index - 1,
                      coarse_error,
                      coarse_rhs,
                      pre_smoothing_steps,
                      post_smoothing_steps,
                      omega);

    const auto fine_error =
      prolongate(transfers[level_index - 1], coarse_error);
    for (std::size_t i = 0; i < x.size(); ++i)
      x[i] += fine_error[i];

    weighted_jacobi_smooth(*level.matrix,
                           x,
                           rhs,
                           post_smoothing_steps,
                           omega);
  }

  void copy_to_function(const std::vector<double> &values,
                        FEMFunction<double, dim> &solution)
  {
    if (values.size() != solution.size())
      throw std::runtime_error("solution size mismatch");
    for (std::size_t i = 0; i < values.size(); ++i)
      solution(i) = values[i];
  }

  TwoGridResult run_pair(const double coarse_h,
                         const std::string &coarse_mesh_file,
                         const double fine_h,
                         const std::string &fine_mesh_file)
  {
    auto coarse = build_level(coarse_h, coarse_mesh_file);
    auto fine = build_level(fine_h, fine_mesh_file);
    const Transfer transfer = build_transfer(fine->space, coarse->space);

    constexpr int    cycles = 8;
    constexpr int    pre_smoothing_steps = 3;
    constexpr int    post_smoothing_steps = 3;
    constexpr double omega = 2.0 / 3.0;
    const int total_smoothing_sweeps =
      cycles * (pre_smoothing_steps + post_smoothing_steps);

    std::vector<double> two_grid_solution(fine->space.n_dof(), 0.0);
    std::vector<double> jacobi_solution(fine->space.n_dof(), 0.0);
    const std::vector<double> fine_rhs = to_std_vector(fine->rhs);
    const double initial_residual =
      relative_residual(*fine->matrix, two_grid_solution, fine_rhs);

    weighted_jacobi_smooth(*fine->matrix,
                           jacobi_solution,
                           fine_rhs,
                           total_smoothing_sweeps,
                           omega);

    for (int cycle = 0; cycle < cycles; ++cycle)
      two_grid_cycle(*fine,
                     *coarse,
                     transfer,
                     two_grid_solution,
                     pre_smoothing_steps,
                     post_smoothing_steps,
                     omega);

    const double jacobi_residual =
      relative_residual(*fine->matrix, jacobi_solution, fine_rhs);
    const double two_grid_residual =
      relative_residual(*fine->matrix, two_grid_solution, fine_rhs);
    if (!std::isfinite(two_grid_residual) ||
        two_grid_residual >= initial_residual)
      throw std::runtime_error("two-grid V-cycle did not reduce the residual");
    if (two_grid_residual >= jacobi_residual)
      throw std::runtime_error(
        "two-grid V-cycle should beat equal-sweep Jacobi smoothing");

    FEMFunction<double, dim> solution(fine->space);
    copy_to_function(two_grid_solution, solution);
    const double l2_error =
      Functional::L2Error(solution, FunctionFunction<double>(&exact_solution), 4);

    return {coarse_h,
            fine_h,
            coarse->space.n_dof(),
            fine->space.n_dof(),
            cycles,
            total_smoothing_sweeps,
            initial_residual,
            jacobi_residual,
            two_grid_residual,
            std::pow(two_grid_residual / initial_residual,
                     1.0 / static_cast<double>(cycles)),
            l2_error};
  }

  std::vector<std::pair<double, std::string>> default_meshes()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.200, mesh_dir + "/unit_square_h0p20"},
            {0.100, mesh_dir + "/unit_square_h0p10"},
            {0.050, mesh_dir + "/unit_square_h0p05"},
            {0.025, mesh_dir + "/unit_square_h0p025"}};
  }

  MultigridResult run_multigrid_case(
    std::vector<std::unique_ptr<Level>> &levels,
    const std::vector<Transfer> &transfers,
    const int finest_level)
  {
    constexpr int    cycles = 8;
    constexpr int    pre_smoothing_steps = 3;
    constexpr int    post_smoothing_steps = 3;
    constexpr double omega = 2.0 / 3.0;
    const int total_smoothing_sweeps =
      cycles * (pre_smoothing_steps + post_smoothing_steps);

    Level &fine = *levels[finest_level];
    const std::vector<double> fine_rhs = to_std_vector(fine.rhs);
    std::vector<double> multigrid_solution(fine.space.n_dof(), 0.0);
    std::vector<double> jacobi_solution(fine.space.n_dof(), 0.0);

    const double initial_residual =
      relative_residual(*fine.matrix, multigrid_solution, fine_rhs);

    weighted_jacobi_smooth(*fine.matrix,
                           jacobi_solution,
                           fine_rhs,
                           total_smoothing_sweeps,
                           omega);

    for (int cycle = 0; cycle < cycles; ++cycle)
      multigrid_v_cycle(levels,
                        transfers,
                        finest_level,
                        multigrid_solution,
                        fine_rhs,
                        pre_smoothing_steps,
                        post_smoothing_steps,
                        omega);

    const double jacobi_residual =
      relative_residual(*fine.matrix, jacobi_solution, fine_rhs);
    const double multigrid_residual =
      relative_residual(*fine.matrix, multigrid_solution, fine_rhs);

    if (!std::isfinite(multigrid_residual) ||
        multigrid_residual >= initial_residual)
      throw std::runtime_error(
        "recursive multigrid V-cycle did not reduce the residual");
    if (multigrid_residual >= jacobi_residual)
      throw std::runtime_error(
        "recursive multigrid should beat equal-sweep Jacobi smoothing");

    FEMFunction<double, dim> solution(fine.space);
    copy_to_function(multigrid_solution, solution);
    const double l2_error =
      Functional::L2Error(solution, FunctionFunction<double>(&exact_solution), 4);

    return {finest_level + 1,
            fine.mesh_h,
            fine.space.n_dof(),
            cycles,
            total_smoothing_sweeps,
            initial_residual,
            jacobi_residual,
            multigrid_residual,
            std::pow(multigrid_residual / initial_residual,
                     1.0 / static_cast<double>(cycles)),
            l2_error};
  }

  std::vector<MultigridResult> run_multigrid_study()
  {
    std::vector<std::unique_ptr<Level>> levels;
    for (const auto &mesh : default_meshes())
      levels.push_back(build_level(mesh.first, mesh.second));

    std::vector<Transfer> transfers;
    for (std::size_t i = 1; i < levels.size(); ++i)
      transfers.push_back(build_transfer(levels[i]->space,
                                         levels[i - 1]->space));

    std::vector<MultigridResult> results;
    for (int finest_level = 1;
         finest_level < static_cast<int>(levels.size());
         ++finest_level)
      results.push_back(run_multigrid_case(levels, transfers, finest_level));
    return results;
  }

  KrylovResult run_krylov_case(
    std::vector<std::unique_ptr<Level>> &levels,
    const std::vector<Transfer> &transfers,
    const int finest_level)
  {
    constexpr int    max_iterations = 500;
    constexpr int    pre_smoothing_steps = 3;
    constexpr int    post_smoothing_steps = 3;
    constexpr double omega = 2.0 / 3.0;
    constexpr double tolerance = 1.0e-8;

    Level &fine = *levels[finest_level];
    const std::vector<double> fine_rhs = to_std_vector(fine.rhs);

    std::vector<double> cg_solution(fine.space.n_dof(), 0.0);
    std::vector<double> jacobi_solution(fine.space.n_dof(), 0.0);
    std::vector<double> multigrid_solution(fine.space.n_dof(), 0.0);

    const auto identity_preconditioner =
      [](const std::vector<double> &src) { return src; };
    const auto jacobi_preconditioner =
      [&fine](const std::vector<double> &src) {
        return apply_jacobi(*fine.matrix, src);
      };
    const auto multigrid_preconditioner =
      [&levels,
       &transfers,
       finest_level,
       pre_smoothing_steps,
       post_smoothing_steps,
       omega](const std::vector<double> &src) {
        std::vector<double> dst(src.size(), 0.0);
        multigrid_v_cycle(levels,
                          transfers,
                          finest_level,
                          dst,
                          src,
                          pre_smoothing_steps,
                          post_smoothing_steps,
                          omega);
        return dst;
      };

    const SolveInfo cg =
      solve_pcg(*fine.matrix,
                cg_solution,
                fine_rhs,
                tolerance,
                max_iterations,
                identity_preconditioner);
    const SolveInfo jacobi_pcg =
      solve_pcg(*fine.matrix,
                jacobi_solution,
                fine_rhs,
                tolerance,
                max_iterations,
                jacobi_preconditioner);
    const SolveInfo multigrid_pcg =
      solve_pcg(*fine.matrix,
                multigrid_solution,
                fine_rhs,
                tolerance,
                max_iterations,
                multigrid_preconditioner);

    if (!cg.converged || !jacobi_pcg.converged || !multigrid_pcg.converged)
      throw std::runtime_error("one of the Poisson PCG solves did not converge");
    if (multigrid_pcg.iterations >= jacobi_pcg.iterations)
      throw std::runtime_error(
        "multigrid-preconditioned PCG should beat Jacobi-PCG");

    FEMFunction<double, dim> solution(fine.space);
    copy_to_function(multigrid_solution, solution);
    const double l2_error =
      Functional::L2Error(solution, FunctionFunction<double>(&exact_solution), 4);

    return {finest_level + 1,
            fine.mesh_h,
            fine.space.n_dof(),
            tolerance,
            cg,
            jacobi_pcg,
            multigrid_pcg,
            l2_error};
  }

  std::vector<KrylovResult> run_krylov_study()
  {
    std::vector<std::unique_ptr<Level>> levels;
    for (const auto &mesh : default_meshes())
      levels.push_back(build_level(mesh.first, mesh.second));

    std::vector<Transfer> transfers;
    for (std::size_t i = 1; i < levels.size(); ++i)
      transfers.push_back(build_transfer(levels[i]->space,
                                         levels[i - 1]->space));

    std::vector<KrylovResult> results;
    for (int finest_level = 1;
         finest_level < static_cast<int>(levels.size());
         ++finest_level)
      results.push_back(run_krylov_case(levels, transfers, finest_level));
    return results;
  }

  void run_default_study()
  {
    const auto meshes = default_meshes();
    std::vector<TwoGridResult> results;
    for (std::size_t i = 1; i < meshes.size(); ++i)
      results.push_back(run_pair(meshes[i - 1].first,
                                 meshes[i - 1].second,
                                 meshes[i].first,
                                 meshes[i].second));

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson geometric two-grid V-cycle diagnostic\n";
    std::cout << "P1 triangles, homogeneous Dirichlet manufactured solution.\n";
    std::cout << std::setw(14) << "H" << std::setw(14) << "h"
              << std::setw(10) << "Nc" << std::setw(10) << "Nf"
              << std::setw(10) << "cycles" << std::setw(10)
              << "Jacobi" << std::setw(14) << "r0" << std::setw(14)
              << "rJacobi" << std::setw(14) << "rTG"
              << std::setw(14) << "TG factor" << std::setw(14)
              << "L2 error" << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(14) << result.coarse_h << std::setw(14)
                  << result.fine_h << std::setw(10) << result.coarse_dofs
                  << std::setw(10) << result.fine_dofs << std::setw(10)
                  << result.cycles << std::setw(10)
                  << result.smoothing_sweeps << std::setw(14)
                  << result.initial_residual << std::setw(14)
                  << result.jacobi_residual << std::setw(14)
                  << result.two_grid_residual << std::setw(14)
                  << result.two_grid_factor << std::setw(14)
                  << result.l2_error << '\n';
      }

    const auto multigrid_results = run_multigrid_study();
    std::cout << "AFEPack Poisson recursive geometric multigrid V-cycle diagnostic\n";
    std::cout << std::setw(10) << "levels" << std::setw(14)
              << "h" << std::setw(10) << "dofs" << std::setw(10)
              << "cycles" << std::setw(10) << "Jacobi"
              << std::setw(14) << "r0" << std::setw(14)
              << "rJacobi" << std::setw(14) << "rMG"
              << std::setw(14) << "MG factor" << std::setw(14)
              << "L2 error" << '\n';
    for (const auto &result : multigrid_results)
      {
        std::cout << std::setw(10) << result.levels << std::setw(14)
                  << result.finest_h << std::setw(10)
                  << result.finest_dofs << std::setw(10)
                  << result.cycles << std::setw(10)
                  << result.smoothing_sweeps << std::setw(14)
                  << result.initial_residual << std::setw(14)
                  << result.jacobi_residual << std::setw(14)
                  << result.multigrid_residual << std::setw(14)
                  << result.multigrid_factor << std::setw(14)
                  << result.l2_error << '\n';
      }

    const auto krylov_results = run_krylov_study();
    std::cout << "AFEPack Poisson recursive multigrid preconditioned CG diagnostic\n";
    std::cout << std::setw(10) << "levels" << std::setw(14)
              << "h" << std::setw(10) << "dofs" << std::setw(14)
              << "tol" << std::setw(10) << "CG" << std::setw(14)
              << "rCG" << std::setw(12) << "Jac-PCG" << std::setw(14)
              << "rJac-PCG" << std::setw(10) << "MG-PCG"
              << std::setw(14) << "rMG-PCG" << std::setw(14)
              << "L2 error" << '\n';
    for (const auto &result : krylov_results)
      {
        std::cout << std::setw(10) << result.levels << std::setw(14)
                  << result.finest_h << std::setw(10)
                  << result.finest_dofs << std::setw(14)
                  << result.tolerance << std::setw(10)
                  << result.cg.iterations << std::setw(14)
                  << result.cg.relative_residual << std::setw(12)
                  << result.jacobi_pcg.iterations << std::setw(14)
                  << result.jacobi_pcg.relative_residual << std::setw(10)
                  << result.multigrid_pcg.iterations << std::setw(14)
                  << result.multigrid_pcg.relative_residual << std::setw(14)
                  << result.l2_error << '\n';
      }
  }
}

int main()
{
  try
    {
      run_default_study();
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
