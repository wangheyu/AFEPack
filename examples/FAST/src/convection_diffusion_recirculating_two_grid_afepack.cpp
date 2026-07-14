/**
 * @file convection_diffusion_recirculating_two_grid_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：回流问题两重网格实验。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 回流问题两重网格实验：演示在粗网格校正和细网格求解之间传递信息的两重网格思路。
 * - 两重网格算法：展示粗细网格组合求解策略在模型问题中的效果。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <functional>
#include <memory>
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

  struct RecirculatingLevel
  {
    double mesh_h = 0.0;
    EasyMesh mesh;
    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof;
    BasisFunctionAdmin<double, dim, dim> basis;
    std::vector<TemplateElement<double, dim, dim>> template_element;
    FEMSpace<double, dim> space;
    std::unique_ptr<ConvectionDiffusionMatrix> matrix;
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

  struct CycleResult
  {
    int          levels = 0;
    double       finest_h = 0.0;
    double       epsilon = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    int          cycles = 0;
    int          macro_sweeps = 0;
    int          directional_passes = 0;
    double       initial_residual = 0.0;
    double       jacobi_residual = 0.0;
    double       left_to_right_residual = 0.0;
    double       four_direction_residual = 0.0;
    double       cycle_residual = 0.0;
    double       cycle_factor = 0.0;
    double       min_value = 0.0;
    double       max_value = 0.0;
    double       undershoot = 0.0;
    double       overshoot = 0.0;
    double       l2_norm = 0.0;
    double       h1_seminorm = 0.0;
  };

  struct KrylovPreconditionerResult
  {
    std::string  preconditioner;
    int          levels = 0;
    double       finest_h = 0.0;
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

  void build_space_in_place(
    EasyMesh &mesh,
    std::vector<TemplateElement<double, dim, dim>> &template_element,
    FEMSpace<double, dim> &space)
  {
    space.reinit(mesh, template_element);

    const int n_element = mesh.n_geometry(dim);
    space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
      space.element(i).reinit(space, i, 0);

    space.buildElement();
    space.buildDof();
    space.buildDofBoundaryMark();
  }

  std::unique_ptr<RecirculatingLevel>
  build_recirculating_level(const double h,
                            const std::string &mesh_file,
                            const double epsilon)
  {
    auto level = std::make_unique<RecirculatingLevel>();
    level->mesh_h = h;
    level->mesh.readData(mesh_file);
    level->template_element = make_triangle_template(level->geometry,
                                                     level->transform,
                                                     level->dof,
                                                     level->basis);
    build_space_in_place(level->mesh, level->template_element, level->space);

    level->matrix =
      std::make_unique<ConvectionDiffusionMatrix>(level->space,
                                                  epsilon,
                                                  true);
    level->matrix->algebricAccuracy() = 4;
    level->matrix->build();

    FEMFunction<double, dim> boundary_solution(level->space);
    Operator::L2Discretize(&right_hand_side, level->space, level->rhs, 4);
    add_supg_rhs(level->space, level->rhs, epsilon, true);

    BoundaryFunction<double, dim> bottom_boundary(
      BoundaryConditionInfo::DIRICHLET, 1, &double_glazing_boundary);
    BoundaryFunction<double, dim> right_boundary(
      BoundaryConditionInfo::DIRICHLET, 2, &double_glazing_boundary);
    BoundaryFunction<double, dim> top_boundary(
      BoundaryConditionInfo::DIRICHLET, 3, &double_glazing_boundary);
    BoundaryFunction<double, dim> left_boundary(
      BoundaryConditionInfo::DIRICHLET, 4, &double_glazing_boundary);

    BoundaryConditionAdmin<double, dim> boundary_admin(level->space);
    boundary_admin.add(bottom_boundary);
    boundary_admin.add(right_boundary);
    boundary_admin.add(top_boundary);
    boundary_admin.add(left_boundary);
    boundary_admin.apply(*level->matrix, boundary_solution, level->rhs);

    return level;
  }

  std::vector<double> to_std_vector_local(const Vector<double> &vector)
  {
    std::vector<double> result(vector.size(), 0.0);
    for (std::size_t i = 0; i < vector.size(); ++i)
      result[i] = vector(i);
    return result;
  }

  double dot_local(const std::vector<double> &a,
                   const std::vector<double> &b)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
      result += a[i] * b[i];
    return result;
  }

  double norm_local(const std::vector<double> &v)
  {
    return std::sqrt(dot_local(v, v));
  }

  std::vector<double> vmult_local(const SparseMatrix<double> &matrix,
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

  std::vector<double> residual_local(const SparseMatrix<double> &matrix,
                                     const std::vector<double> &x,
                                     const std::vector<double> &rhs)
  {
    std::vector<double> result = rhs;
    const auto ax = vmult_local(matrix, x);
    for (std::size_t i = 0; i < result.size(); ++i)
      result[i] -= ax[i];
    return result;
  }

  double relative_residual_local(const SparseMatrix<double> &matrix,
                                 const std::vector<double> &x,
                                 const std::vector<double> &rhs)
  {
    return norm_local(residual_local(matrix, x, rhs)) /
           std::max(norm_local(rhs), 1.0e-30);
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

                const bool horizontal =
                  direction == SweepDirection::left_to_right ||
                  direction == SweepDirection::right_to_left;
                const double primary_l = horizontal ? pl[0] : pl[1];
                const double primary_r = horizontal ? pr[0] : pr[1];
                const double secondary_l = horizontal ? pl[1] : pl[0];
                const double secondary_r = horizontal ? pr[1] : pr[0];

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

  std::vector<std::vector<int>>
  four_direction_orders(const FEMSpace<double, dim> &space)
  {
    return {dof_order(space, SweepDirection::left_to_right),
            dof_order(space, SweepDirection::right_to_left),
            dof_order(space, SweepDirection::bottom_to_top),
            dof_order(space, SweepDirection::top_to_bottom)};
  }

  void ordered_gauss_seidel_pass(const SparseMatrix<double> &matrix,
                                 const std::vector<int> &order,
                                 std::vector<double> &x,
                                 const std::vector<double> &rhs)
  {
    const auto &pattern = matrix.get_sparsity_pattern();
    const size_t *row_start = pattern.get_rowstart_indices();
    const size_t *columns = pattern.get_column_numbers();

    for (const int row : order)
      {
        double sum = 0.0;
        for (std::size_t index = row_start[row]; index < row_start[row + 1];
             ++index)
          {
            const int column = static_cast<int>(columns[index]);
            if (column != row)
              sum += matrix.global_entry(index) * x[column];
          }

        const double diagonal = matrix.diag_element(row);
        if (std::abs(diagonal) > 1.0e-14)
          x[row] = (rhs[row] - sum) / diagonal;
      }
  }

  void directional_gauss_seidel_smooth(
    const SparseMatrix<double> &matrix,
    const std::vector<std::vector<int>> &orders,
    std::vector<double> &x,
    const std::vector<double> &rhs,
    const int macro_sweeps)
  {
    for (int sweep = 0; sweep < macro_sweeps; ++sweep)
      for (const auto &order : orders)
        ordered_gauss_seidel_pass(matrix, order, x, rhs);
  }

  void weighted_jacobi_smooth_local(const SparseMatrix<double> &matrix,
                                    std::vector<double> &x,
                                    const std::vector<double> &rhs,
                                    const int sweeps,
                                    const double omega)
  {
    for (int sweep = 0; sweep < sweeps; ++sweep)
      {
        const auto r = residual_local(matrix, x, rhs);
        for (std::size_t i = 0; i < x.size(); ++i)
          {
            const double diagonal = matrix.diag_element(i);
            if (std::abs(diagonal) > 1.0e-14)
              x[i] += omega * r[i] / diagonal;
          }
      }
  }

  std::vector<double> dense_solve(const SparseMatrix<double> &matrix,
                                  const std::vector<double> &rhs)
  {
    const int n = static_cast<int>(rhs.size());
    std::vector<std::vector<double>> a(
      static_cast<std::size_t>(n),
      std::vector<double>(static_cast<std::size_t>(n), 0.0));
    std::vector<double> b = rhs;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        a[i][j] = matrix.el(i, j);

    for (int pivot = 0; pivot < n; ++pivot)
      {
        int best = pivot;
        for (int row = pivot + 1; row < n; ++row)
          if (std::abs(a[row][pivot]) > std::abs(a[best][pivot]))
            best = row;
        if (std::abs(a[best][pivot]) < 1.0e-20)
          throw std::runtime_error("singular coarse matrix in dense solve");
        if (best != pivot)
          {
            std::swap(a[best], a[pivot]);
            std::swap(b[best], b[pivot]);
          }

        for (int row = pivot + 1; row < n; ++row)
          {
            const double factor = a[row][pivot] / a[pivot][pivot];
            a[row][pivot] = 0.0;
            for (int col = pivot + 1; col < n; ++col)
              a[row][col] -= factor * a[pivot][col];
            b[row] -= factor * b[pivot];
          }
      }

    std::vector<double> x(static_cast<std::size_t>(n), 0.0);
    for (int row = n - 1; row >= 0; --row)
      {
        double sum = b[row];
        for (int col = row + 1; col < n; ++col)
          sum -= a[row][col] * x[col];
        x[row] = sum / a[row][row];
      }
    return x;
  }

  bool barycentric_coordinates_local(const Point<dim> &p,
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
  interpolation_row_local(const FEMSpace<double, dim> &coarse_space,
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
        if (!barycentric_coordinates_local(point,
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

  Transfer build_transfer_local(const FEMSpace<double, dim> &fine_space,
                                const FEMSpace<double, dim> &coarse_space)
  {
    Transfer transfer;
    transfer.coarse_size = coarse_space.n_dof();
    transfer.rows.resize(fine_space.n_dof());

    for (int i = 0; i < fine_space.n_dof(); ++i)
      transfer.rows[i] =
        interpolation_row_local(coarse_space,
                                fine_space.dofInfo(i).interp_point);

    return transfer;
  }

  std::vector<double> restrict_residual_local(
    const Transfer &transfer,
    const std::vector<double> &fine_vector)
  {
    std::vector<double> coarse_vector(transfer.coarse_size, 0.0);
    for (std::size_t i = 0; i < transfer.rows.size(); ++i)
      for (const auto &entry : transfer.rows[i])
        coarse_vector[entry.coarse_index] += entry.weight * fine_vector[i];
    return coarse_vector;
  }

  std::vector<double> prolongate_local(
    const Transfer &transfer,
    const std::vector<double> &coarse_vector)
  {
    std::vector<double> fine_vector(transfer.rows.size(), 0.0);
    for (std::size_t i = 0; i < transfer.rows.size(); ++i)
      for (const auto &entry : transfer.rows[i])
        fine_vector[i] += entry.weight * coarse_vector[entry.coarse_index];
    return fine_vector;
  }

  void v_cycle(
    const std::vector<std::unique_ptr<RecirculatingLevel>> &levels,
    const std::vector<Transfer> &transfers,
    const std::vector<std::vector<std::vector<int>>> &level_orders,
    const int level_index,
    std::vector<double> &x,
    const std::vector<double> &rhs,
    const int pre_smoothing_steps,
    const int post_smoothing_steps)
  {
    const auto &level = *levels[level_index];

    if (level_index == 0)
      {
        x = dense_solve(*level.matrix, rhs);
        return;
      }

    directional_gauss_seidel_smooth(*level.matrix,
                                    level_orders[level_index],
                                    x,
                                    rhs,
                                    pre_smoothing_steps);

    const auto fine_residual = residual_local(*level.matrix, x, rhs);
    const auto coarse_rhs =
      restrict_residual_local(transfers[level_index - 1], fine_residual);

    std::vector<double> coarse_error(
      levels[level_index - 1]->space.n_dof(), 0.0);
    v_cycle(levels,
            transfers,
            level_orders,
            level_index - 1,
            coarse_error,
            coarse_rhs,
            pre_smoothing_steps,
            post_smoothing_steps);

    const auto fine_error =
      prolongate_local(transfers[level_index - 1], coarse_error);
    for (std::size_t i = 0; i < x.size(); ++i)
      x[i] += fine_error[i];

    directional_gauss_seidel_smooth(*level.matrix,
                                    level_orders[level_index],
                                    x,
                                    rhs,
                                    post_smoothing_steps);
  }

  void copy_to_function_local(const std::vector<double> &values,
                              FEMFunction<double, dim> &solution)
  {
    if (values.size() != solution.size())
      throw std::runtime_error("solution size mismatch");
    for (std::size_t i = 0; i < values.size(); ++i)
      solution(i) = values[i];
  }

  void copy_to_vector_local(const std::vector<double> &values,
                            Vector<double> &vector)
  {
    vector.reinit(values.size(), false);
    for (std::size_t i = 0; i < values.size(); ++i)
      vector(i) = values[i];
  }

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

  std::vector<std::pair<double, std::string>> default_meshes_local()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.20, mesh_dir + "/unit_square_h0p20"},
            {0.10, mesh_dir + "/unit_square_h0p10"},
            {0.05, mesh_dir + "/unit_square_h0p05"}};
  }

  CycleResult run_cycle_case(
    std::vector<std::unique_ptr<RecirculatingLevel>> &levels,
    const std::vector<Transfer> &transfers,
    const std::vector<std::vector<std::vector<int>>> &four_orders,
    const std::vector<std::vector<std::vector<int>>> &left_to_right_orders,
    const int finest_level,
    const double epsilon)
  {
    constexpr int cycles = 6;
    constexpr int pre_smoothing_steps = 1;
    constexpr int post_smoothing_steps = 1;
    const int macro_sweeps =
      cycles * (pre_smoothing_steps + post_smoothing_steps);
    const int directional_passes = 4 * macro_sweeps;

    auto &fine = *levels[finest_level];
    const auto fine_rhs = to_std_vector_local(fine.rhs);
    std::vector<double> jacobi_solution(fine.space.n_dof(), 0.0);
    std::vector<double> single_gs_solution(fine.space.n_dof(), 0.0);
    std::vector<double> four_gs_solution(fine.space.n_dof(), 0.0);
    std::vector<double> cycle_solution(fine.space.n_dof(), 0.0);

    const double initial_residual =
      relative_residual_local(*fine.matrix, cycle_solution, fine_rhs);

    weighted_jacobi_smooth_local(*fine.matrix,
                                 jacobi_solution,
                                 fine_rhs,
                                 directional_passes,
                                 0.7);

    directional_gauss_seidel_smooth(*fine.matrix,
                                    left_to_right_orders[finest_level],
                                    single_gs_solution,
                                    fine_rhs,
                                    directional_passes);

    directional_gauss_seidel_smooth(*fine.matrix,
                                    four_orders[finest_level],
                                    four_gs_solution,
                                    fine_rhs,
                                    macro_sweeps);

    for (int cycle = 0; cycle < cycles; ++cycle)
      v_cycle(levels,
              transfers,
              four_orders,
              finest_level,
              cycle_solution,
              fine_rhs,
              pre_smoothing_steps,
              post_smoothing_steps);

    const double jacobi_residual =
      relative_residual_local(*fine.matrix, jacobi_solution, fine_rhs);
    const double single_gs_residual =
      relative_residual_local(*fine.matrix, single_gs_solution, fine_rhs);
    const double four_gs_residual =
      relative_residual_local(*fine.matrix, four_gs_solution, fine_rhs);
    const double cycle_residual =
      relative_residual_local(*fine.matrix, cycle_solution, fine_rhs);

    if (!std::isfinite(cycle_residual) ||
        cycle_residual >= initial_residual)
      throw std::runtime_error(
        "recirculating convection-diffusion V-cycle did not reduce");

    FEMFunction<double, dim> solution(fine.space);
    copy_to_function_local(cycle_solution, solution);
    const SolutionDiagnostics diagnostics = diagnose_solution(solution);

    return {finest_level + 1,
            fine.mesh_h,
            epsilon,
            static_cast<int>(fine.mesh.n_geometry(dim)),
            fine.space.n_dof(),
            cycles,
            macro_sweeps,
            directional_passes,
            initial_residual,
            jacobi_residual,
            single_gs_residual,
            four_gs_residual,
            cycle_residual,
            std::pow(cycle_residual / initial_residual,
                     1.0 / static_cast<double>(cycles)),
            diagnostics.min_value,
            diagnostics.max_value,
            diagnostics.undershoot,
            diagnostics.overshoot,
            diagnostics.l2_norm,
            diagnostics.h1_seminorm};
  }

  std::vector<CycleResult> run_default_study(const double epsilon)
  {
    std::vector<std::unique_ptr<RecirculatingLevel>> levels;
    for (const auto &mesh : default_meshes_local())
      levels.push_back(build_recirculating_level(mesh.first,
                                                 mesh.second,
                                                 epsilon));

    std::vector<Transfer> transfers;
    for (std::size_t i = 1; i < levels.size(); ++i)
      transfers.push_back(build_transfer_local(levels[i]->space,
                                               levels[i - 1]->space));

    std::vector<std::vector<std::vector<int>>> four_orders;
    std::vector<std::vector<std::vector<int>>> left_to_right_orders;
    for (const auto &level : levels)
      {
        four_orders.push_back(four_direction_orders(level->space));
        left_to_right_orders.push_back(
          {dof_order(level->space, SweepDirection::left_to_right)});
      }

    std::vector<CycleResult> results;
    for (int finest_level = 1;
         finest_level < static_cast<int>(levels.size());
         ++finest_level)
      results.push_back(run_cycle_case(levels,
                                       transfers,
                                       four_orders,
                                       left_to_right_orders,
                                       finest_level,
                                       epsilon));
    return results;
  }

  KrylovPreconditionerResult run_krylov_case(
    std::vector<std::unique_ptr<RecirculatingLevel>> &levels,
    const std::vector<Transfer> &transfers,
    const std::vector<std::vector<std::vector<int>>> &four_orders,
    const int finest_level,
    const double epsilon,
    const std::string &preconditioner)
  {
    auto &fine = *levels[finest_level];
    FEMFunction<double, dim> solution(fine.space);
    for (std::size_t i = 0; i < solution.size(); ++i)
      solution(i) = 0.0;

    std::function<void(const Vector<double> &, Vector<double> &)> apply =
      [](const Vector<double> &src, Vector<double> &dst) { dst = src; };

    if (preconditioner == "Jacobi")
      {
        apply = [&](const Vector<double> &src, Vector<double> &dst) {
          jacobi_apply(*fine.matrix, src, dst);
        };
      }
    else if (preconditioner == "four-direction-GS")
      {
        apply = [&](const Vector<double> &src, Vector<double> &dst) {
          std::vector<double> z(src.size(), 0.0);
          directional_gauss_seidel_smooth(*fine.matrix,
                                          four_orders[finest_level],
                                          z,
                                          to_std_vector_local(src),
                                          1);
          copy_to_vector_local(z, dst);
        };
      }
    else if (preconditioner == "one-V-cycle")
      {
        apply = [&](const Vector<double> &src, Vector<double> &dst) {
          std::vector<double> z(src.size(), 0.0);
          v_cycle(levels,
                  transfers,
                  four_orders,
                  finest_level,
                  z,
                  to_std_vector_local(src),
                  1,
                  1);
          copy_to_vector_local(z, dst);
        };
      }

    const SolveInfo solve_info =
      solve_right_preconditioned_gmres(*fine.matrix,
                                       solution,
                                       fine.rhs,
                                       apply,
                                       1.0e-8,
                                       50,
                                       2000);

    const SolutionDiagnostics diagnostics = diagnose_solution(solution);
    return {preconditioner,
            finest_level + 1,
            fine.mesh_h,
            epsilon,
            static_cast<int>(fine.mesh.n_geometry(dim)),
            fine.space.n_dof(),
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

  std::vector<KrylovPreconditionerResult>
  run_krylov_study(const double epsilon)
  {
    std::vector<std::unique_ptr<RecirculatingLevel>> levels;
    for (const auto &mesh : default_meshes_local())
      levels.push_back(build_recirculating_level(mesh.first,
                                                 mesh.second,
                                                 epsilon));

    std::vector<Transfer> transfers;
    for (std::size_t i = 1; i < levels.size(); ++i)
      transfers.push_back(build_transfer_local(levels[i]->space,
                                               levels[i - 1]->space));

    std::vector<std::vector<std::vector<int>>> four_orders;
    for (const auto &level : levels)
      four_orders.push_back(four_direction_orders(level->space));

    const int finest_level = static_cast<int>(levels.size()) - 1;
    std::vector<KrylovPreconditionerResult> results;
    for (const std::string &preconditioner :
         {"Jacobi", "four-direction-GS", "one-V-cycle"})
      results.push_back(run_krylov_case(levels,
                                        transfers,
                                        four_orders,
                                        finest_level,
                                        epsilon,
                                        preconditioner));
    return results;
  }

  void print_results(const std::vector<CycleResult> &results)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack recirculating double-glazing SUPG geometric "
                 "V-cycle diagnostic\n";
    std::cout << "Book Example 3.1.4; V-cycle smoother is four-direction "
                 "Gauss-Seidel, coarse solves use dense pivoted elimination.\n";
    std::cout << "levels h epsilon cells dofs cycles macro_sweeps "
                 "directional_passes r0 rJacobi rLeftGS rFourGS rV "
                 "Vfactor min max undershoot overshoot L2 H1\n";
    for (const auto &result : results)
      std::cout << result.levels << ' '
                << result.finest_h << ' '
                << result.epsilon << ' '
                << result.cells << ' '
                << result.dofs << ' '
                << result.cycles << ' '
                << result.macro_sweeps << ' '
                << result.directional_passes << ' '
                << result.initial_residual << ' '
                << result.jacobi_residual << ' '
                << result.left_to_right_residual << ' '
                << result.four_direction_residual << ' '
                << result.cycle_residual << ' '
                << result.cycle_factor << ' '
                << result.min_value << ' '
                << result.max_value << ' '
                << result.undershoot << ' '
                << result.overshoot << ' '
                << result.l2_norm << ' '
                << result.h1_seminorm << '\n';
  }

  void print_krylov_results(
    const std::vector<KrylovPreconditionerResult> &results)
  {
    std::cout << "AFEPack recirculating double-glazing right-preconditioned "
                 "GMRES with geometric V-cycle\n";
    std::cout << "levels h epsilon cells dofs preconditioner GMRES status "
                 "relres min max undershoot overshoot L2 H1\n";
    for (const auto &result : results)
      std::cout << result.levels << ' '
                << result.finest_h << ' '
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
      std::vector<CycleResult> results;
      for (const double epsilon : {1.0 / 200.0, 1.0 / 500.0})
        {
          const auto epsilon_results = run_default_study(epsilon);
          results.insert(results.end(),
                         epsilon_results.begin(),
                         epsilon_results.end());
        }
      print_results(results);

      std::vector<KrylovPreconditionerResult> krylov_results;
      for (const double epsilon : {1.0 / 200.0, 1.0 / 500.0})
        {
          const auto epsilon_results = run_krylov_study(epsilon);
          krylov_results.insert(krylov_results.end(),
                                epsilon_results.begin(),
                                epsilon_results.end());
        }
      print_krylov_results(krylov_results);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
