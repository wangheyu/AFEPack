/**
 * @file convection_diffusion_recirculating_sweep_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：基础离散流程。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 基础离散流程：给出模型问题从网格读取、有限元空间构造、线性系统装配到结果输出的完整流程。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <functional>
#include <numeric>
#include <utility>

#define main convection_diffusion_recirculating_hidden_main
#include "convection_diffusion_recirculating_layer_afepack.cpp"
#undef main

namespace
{
  enum class SweepDirection
  {
    left_to_right,
    right_to_left,
    bottom_to_top,
    top_to_bottom
  };

  enum class RecirculatingPreconditioner
  {
    none,
    jacobi,
    left_to_right_gs,
    four_direction_gs
  };

  struct RecirculatingSweepResult
  {
    std::string  preconditioner;
    double       h = 0.0;
    double       epsilon = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    int          iterations = 0;
    bool         converged = false;
    double       relative_residual = 0.0;
    double       min_value = 0.0;
    double       max_value = 0.0;
    double       undershoot = 0.0;
    double       overshoot = 0.0;
    double       l2_norm = 0.0;
    double       h1_seminorm = 0.0;
  };

  const char *preconditioner_name(
    const RecirculatingPreconditioner preconditioner)
  {
    switch (preconditioner)
      {
        case RecirculatingPreconditioner::none:
          return "none";
        case RecirculatingPreconditioner::jacobi:
          return "Jacobi";
        case RecirculatingPreconditioner::left_to_right_gs:
          return "left-to-right-GS";
        case RecirculatingPreconditioner::four_direction_gs:
          return "four-direction-GS";
      }
    return "unknown";
  }

  std::vector<int> dof_order(const FEMSpace<double, dim> &space,
                             const SweepDirection direction)
  {
    std::vector<int> order(static_cast<std::size_t>(space.n_dof()));
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(),
              order.end(),
              [&](const int left, const int right) {
                const auto &pl = space.dofInfo(left).interp_point;
                const auto &pr = space.dofInfo(right).interp_point;

                const double primary_l =
                  direction == SweepDirection::left_to_right ||
                      direction == SweepDirection::right_to_left ?
                    pl[0] :
                    pl[1];
                const double primary_r =
                  direction == SweepDirection::left_to_right ||
                      direction == SweepDirection::right_to_left ?
                    pr[0] :
                    pr[1];
                const double secondary_l =
                  direction == SweepDirection::left_to_right ||
                      direction == SweepDirection::right_to_left ?
                    pl[1] :
                    pl[0];
                const double secondary_r =
                  direction == SweepDirection::left_to_right ||
                      direction == SweepDirection::right_to_left ?
                    pr[1] :
                    pr[0];

                if (std::abs(primary_l - primary_r) > 1.0e-12)
                  {
                    const bool increasing =
                      direction == SweepDirection::left_to_right ||
                      direction == SweepDirection::bottom_to_top;
                    return increasing ? primary_l < primary_r :
                                        primary_l > primary_r;
                  }
                if (std::abs(secondary_l - secondary_r) > 1.0e-12)
                  return secondary_l < secondary_r;
                return left < right;
              });
    return order;
  }

  class DirectionalGaussSeidelPreconditioner
  {
  public:
    DirectionalGaussSeidelPreconditioner(
      const SparseMatrix<double> &matrix,
      std::vector<std::vector<int>> orders)
      : orders(std::move(orders)),
        entries(this->orders.front().size(),
                std::vector<double>(this->orders.front().size(), 0.0))
    {
      for (std::size_t i = 0; i < entries.size(); ++i)
        for (std::size_t j = 0; j < entries.size(); ++j)
          entries[i][j] = matrix.el(i, j);
    }

    void apply(const Vector<double> &src, Vector<double> &dst) const
    {
      dst.reinit(src.size(), false);
      dst = 0.0;

      for (const auto &order : orders)
        for (const int row : order)
          {
            double value = src(row);
            for (std::size_t column = 0; column < entries.size(); ++column)
              if (static_cast<int>(column) != row)
                value -= entries[row][column] * dst(column);

            const double diagonal = entries[row][row];
            dst(row) = std::abs(diagonal) > 1.0e-14 ? value / diagonal :
                                                        value;
          }
    }

  private:
    std::vector<std::vector<int>>    orders;
    std::vector<std::vector<double>> entries;
  };

  SolveInfo solve_right_preconditioned_gmres(
    const SparseMatrix<double> &matrix,
    Vector<double> &x,
    const Vector<double> &rhs,
    const std::function<void(const Vector<double> &, Vector<double> &)>
      &apply_preconditioner,
    const double tolerance,
    const int restart,
    const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);
    Vector<double> r = residual(matrix, x, rhs);
    double beta = r.l2_norm();
    double relative_residual = beta / rhs_norm;
    if (relative_residual <= tolerance)
      return {true, 0, relative_residual};

    int total_iterations = 0;
    while (total_iterations < max_iterations)
      {
        std::vector<Vector<double>> v(restart + 1, Vector<double>(rhs.size()));
        std::vector<Vector<double>> z(restart, Vector<double>(rhs.size()));
        std::vector<std::vector<double>> h(
          restart + 1, std::vector<double>(restart, 0.0));
        std::vector<double> cs(restart, 0.0);
        std::vector<double> sn(restart, 0.0);
        std::vector<double> g(restart + 1, 0.0);

        v[0] = r;
        normalize(v[0], beta);
        g[0] = beta;

        int inner_iterations = 0;
        for (; inner_iterations < restart &&
               total_iterations < max_iterations;
             ++inner_iterations, ++total_iterations)
          {
            apply_preconditioner(v[inner_iterations],
                                 z[inner_iterations]);
            matrix.vmult(v[inner_iterations + 1],
                         z[inner_iterations]);

            for (int i = 0; i <= inner_iterations; ++i)
              {
                h[i][inner_iterations] =
                  v[inner_iterations + 1].dot(v[i]);
                add_scaled(v[inner_iterations + 1],
                           -h[i][inner_iterations],
                           v[i]);
              }

            h[inner_iterations + 1][inner_iterations] =
              v[inner_iterations + 1].l2_norm();
            if (h[inner_iterations + 1][inner_iterations] != 0.0)
              normalize(v[inner_iterations + 1],
                        h[inner_iterations + 1][inner_iterations]);

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

        r = residual(matrix, x, rhs);
        beta = r.l2_norm();
        relative_residual = beta / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }

  RecirculatingSweepResult run_recirculating_sweep_case(
    const std::string &mesh_file,
    const double h,
    const double epsilon,
    const RecirculatingPreconditioner preconditioner)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);
    auto fem_space = build_space(mesh, template_element);

    ConvectionDiffusionMatrix matrix(fem_space, epsilon, true);
    matrix.algebricAccuracy() = 4;
    matrix.build();

    FEMFunction<double, dim> solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    add_supg_rhs(fem_space, rhs, epsilon, true);

    BoundaryFunction<double, dim> bottom_boundary(
      BoundaryConditionInfo::DIRICHLET, 1, &double_glazing_boundary);
    BoundaryFunction<double, dim> right_boundary(
      BoundaryConditionInfo::DIRICHLET, 2, &double_glazing_boundary);
    BoundaryFunction<double, dim> top_boundary(
      BoundaryConditionInfo::DIRICHLET, 3, &double_glazing_boundary);
    BoundaryFunction<double, dim> left_boundary(
      BoundaryConditionInfo::DIRICHLET, 4, &double_glazing_boundary);

    BoundaryConditionAdmin<double, dim> boundary_admin(fem_space);
    boundary_admin.add(bottom_boundary);
    boundary_admin.add(right_boundary);
    boundary_admin.add(top_boundary);
    boundary_admin.add(left_boundary);
    boundary_admin.apply(matrix, solution, rhs);

    for (std::size_t i = 0; i < solution.size(); ++i)
      solution(i) = 0.0;

    const std::vector<int> left_to_right =
      dof_order(fem_space, SweepDirection::left_to_right);
    DirectionalGaussSeidelPreconditioner single_sweep(
      matrix,
      {left_to_right});
    DirectionalGaussSeidelPreconditioner four_sweep(
      matrix,
      {left_to_right,
       dof_order(fem_space, SweepDirection::right_to_left),
       dof_order(fem_space, SweepDirection::bottom_to_top),
       dof_order(fem_space, SweepDirection::top_to_bottom)});

    std::function<void(const Vector<double> &, Vector<double> &)> apply =
      [](const Vector<double> &src, Vector<double> &dst) { dst = src; };

    switch (preconditioner)
      {
        case RecirculatingPreconditioner::none:
          apply = [](const Vector<double> &src, Vector<double> &dst) {
            dst = src;
          };
          break;
        case RecirculatingPreconditioner::jacobi:
          apply = [&](const Vector<double> &src, Vector<double> &dst) {
            jacobi_apply(matrix, src, dst);
          };
          break;
        case RecirculatingPreconditioner::left_to_right_gs:
          apply = [&](const Vector<double> &src, Vector<double> &dst) {
            single_sweep.apply(src, dst);
          };
          break;
        case RecirculatingPreconditioner::four_direction_gs:
          apply = [&](const Vector<double> &src, Vector<double> &dst) {
            four_sweep.apply(src, dst);
          };
          break;
      }

    const SolveInfo solve_info =
      solve_right_preconditioned_gmres(matrix,
                                       solution,
                                       rhs,
                                       apply,
                                       1.0e-8,
                                       50,
                                       5000);
    const auto diagnostics = diagnose_solution(solution);
    return {preconditioner_name(preconditioner),
            h,
            epsilon,
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            solve_info.iterations,
            solve_info.converged,
            solve_info.relative_residual,
            diagnostics.min_value,
            diagnostics.max_value,
            diagnostics.undershoot,
            diagnostics.overshoot,
            diagnostics.l2_norm,
            diagnostics.h1_seminorm};
  }

  void print_results(const std::vector<RecirculatingSweepResult> &results)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack recirculating double-glazing SUPG sweep "
                 "preconditioner test\n";
    std::cout << "Book Example 3.1.4 with Chapter 4 four-direction "
                 "Gauss-Seidel idea; status=limit is reported, not fatal.\n";
    std::cout << "h epsilon cells dofs preconditioner GMRES status relres "
                 "min max undershoot overshoot L2 H1\n";
    for (const auto &result : results)
      std::cout << result.h << ' '
                << result.epsilon << ' '
                << result.cells << ' '
                << result.dofs << ' '
                << result.preconditioner << ' '
                << result.iterations << ' '
                << (result.converged ? "ok" : "limit") << ' '
                << result.relative_residual << ' '
                << result.min_value << ' '
                << result.max_value << ' '
                << result.undershoot << ' '
                << result.overshoot << ' '
                << result.l2_norm << ' '
                << result.h1_seminorm << '\n';
  }
}

int main()
{
  try
    {
      const std::string mesh_file =
        std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p10";
      const std::vector<double> epsilons = {1.0 / 100.0,
                                            1.0 / 200.0,
                                            1.0 / 500.0};
      const std::vector<RecirculatingPreconditioner> preconditioners = {
        RecirculatingPreconditioner::none,
        RecirculatingPreconditioner::jacobi,
        RecirculatingPreconditioner::left_to_right_gs,
        RecirculatingPreconditioner::four_direction_gs};

      std::vector<RecirculatingSweepResult> results;
      for (const double epsilon : epsilons)
        for (const auto preconditioner : preconditioners)
          results.push_back(run_recirculating_sweep_case(mesh_file,
                                                         0.10,
                                                         epsilon,
                                                         preconditioner));
      print_results(results);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
