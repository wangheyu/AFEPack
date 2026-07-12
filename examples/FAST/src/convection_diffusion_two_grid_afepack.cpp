/**
 * @file convection_diffusion_two_grid_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：两重网格算法。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
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

#include <memory>
#include <numeric>
#include <utility>

#define main convection_diffusion_boundary_layer_hidden_main
#include "convection_diffusion_boundary_layer_afepack.cpp"
#undef main

namespace
{
  struct ConvectionDiffusionLevel
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
    unsigned int finest_dofs = 0;
    int          cycles = 0;
    int          smoothing_sweeps = 0;
    double       initial_residual = 0.0;
    double       jacobi_residual = 0.0;
    double       gs_residual = 0.0;
    double       cycle_residual = 0.0;
    double       cycle_factor = 0.0;
    double       l2_error = 0.0;
    double       h1_error = 0.0;
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

  std::unique_ptr<ConvectionDiffusionLevel>
  build_cd_level(const double h,
                 const std::string &mesh_file,
                 const double epsilon)
  {
    current_epsilon = epsilon;

    auto level = std::make_unique<ConvectionDiffusionLevel>();
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

    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &exact_solution);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(level->space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
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

  std::vector<int> dof_order_by_x_local(const FEMSpace<double, dim> &space,
                                        const bool downstream)
  {
    std::vector<int> order(static_cast<std::size_t>(space.n_dof()));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(),
              order.end(),
              [&](const int left, const int right) {
                const auto &left_point = space.dofInfo(left).interp_point;
                const auto &right_point = space.dofInfo(right).interp_point;
                if (std::abs(left_point[0] - right_point[0]) > 1.0e-12)
                  return downstream ? left_point[0] < right_point[0] :
                                      left_point[0] > right_point[0];
                if (std::abs(left_point[1] - right_point[1]) > 1.0e-12)
                  return left_point[1] < right_point[1];
                return left < right;
              });
    return order;
  }

  void ordered_gauss_seidel_smooth(const SparseMatrix<double> &matrix,
                                   const std::vector<int> &order,
                                   std::vector<double> &x,
                                   const std::vector<double> &rhs,
                                   const int sweeps)
  {
    for (int sweep = 0; sweep < sweeps; ++sweep)
      for (const int row : order)
        {
          double sum = 0.0;
          const auto &pattern = matrix.get_sparsity_pattern();
          const size_t *row_start = pattern.get_rowstart_indices();
          const size_t *columns = pattern.get_column_numbers();
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
    const std::vector<std::unique_ptr<ConvectionDiffusionLevel>> &levels,
    const std::vector<Transfer> &transfers,
    const std::vector<std::vector<int>> &downstream_orders,
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

    ordered_gauss_seidel_smooth(*level.matrix,
                                downstream_orders[level_index],
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
            downstream_orders,
            level_index - 1,
            coarse_error,
            coarse_rhs,
            pre_smoothing_steps,
            post_smoothing_steps);

    const auto fine_error =
      prolongate_local(transfers[level_index - 1], coarse_error);
    for (std::size_t i = 0; i < x.size(); ++i)
      x[i] += fine_error[i];

    ordered_gauss_seidel_smooth(*level.matrix,
                                downstream_orders[level_index],
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

  std::vector<std::pair<double, std::string>> default_meshes_local()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.20, mesh_dir + "/unit_square_h0p20"},
            {0.10, mesh_dir + "/unit_square_h0p10"},
            {0.05, mesh_dir + "/unit_square_h0p05"}};
  }

  CycleResult run_cycle_case(
    std::vector<std::unique_ptr<ConvectionDiffusionLevel>> &levels,
    const std::vector<Transfer> &transfers,
    const std::vector<std::vector<int>> &downstream_orders,
    const int finest_level,
    const double epsilon)
  {
    current_epsilon = epsilon;

    constexpr int cycles = 8;
    constexpr int pre_smoothing_steps = 2;
    constexpr int post_smoothing_steps = 2;
    const int total_smoothing_sweeps =
      cycles * (pre_smoothing_steps + post_smoothing_steps);

    auto &fine = *levels[finest_level];
    const auto fine_rhs = to_std_vector_local(fine.rhs);
    std::vector<double> jacobi_solution(fine.space.n_dof(), 0.0);
    std::vector<double> gs_solution(fine.space.n_dof(), 0.0);
    std::vector<double> cycle_solution(fine.space.n_dof(), 0.0);

    const double initial_residual =
      relative_residual_local(*fine.matrix, cycle_solution, fine_rhs);

    weighted_jacobi_smooth_local(*fine.matrix,
                                 jacobi_solution,
                                 fine_rhs,
                                 total_smoothing_sweeps,
                                 0.7);

    ordered_gauss_seidel_smooth(*fine.matrix,
                                downstream_orders[finest_level],
                                gs_solution,
                                fine_rhs,
                                total_smoothing_sweeps);

    for (int cycle = 0; cycle < cycles; ++cycle)
      v_cycle(levels,
              transfers,
              downstream_orders,
              finest_level,
              cycle_solution,
              fine_rhs,
              pre_smoothing_steps,
              post_smoothing_steps);

    const double jacobi_residual =
      relative_residual_local(*fine.matrix, jacobi_solution, fine_rhs);
    const double gs_residual =
      relative_residual_local(*fine.matrix, gs_solution, fine_rhs);
    const double cycle_residual =
      relative_residual_local(*fine.matrix, cycle_solution, fine_rhs);

    if (!std::isfinite(cycle_residual) ||
        cycle_residual >= initial_residual)
      throw std::runtime_error("convection-diffusion V-cycle did not reduce");

    FEMFunction<double, dim> solution(fine.space);
    copy_to_function_local(cycle_solution, solution);
    FunctionFunction<double> exact(&exact_solution, &exact_gradient);

    return {finest_level + 1,
            fine.mesh_h,
            fine.space.n_dof(),
            cycles,
            total_smoothing_sweeps,
            initial_residual,
            jacobi_residual,
            gs_residual,
            cycle_residual,
            std::pow(cycle_residual / initial_residual,
                     1.0 / static_cast<double>(cycles)),
            Functional::L2Error(solution, exact, 4),
            Functional::H1SemiError(solution, exact, 4)};
  }

  void run_default_study(const double epsilon)
  {
    current_epsilon = epsilon;

    std::vector<std::unique_ptr<ConvectionDiffusionLevel>> levels;
    for (const auto &mesh : default_meshes_local())
      levels.push_back(build_cd_level(mesh.first, mesh.second, epsilon));

    std::vector<Transfer> transfers;
    for (std::size_t i = 1; i < levels.size(); ++i)
      transfers.push_back(build_transfer_local(levels[i]->space,
                                               levels[i - 1]->space));

    std::vector<std::vector<int>> downstream_orders;
    for (const auto &level : levels)
      downstream_orders.push_back(dof_order_by_x_local(level->space, true));

    std::vector<CycleResult> results;
    for (int finest_level = 1;
         finest_level < static_cast<int>(levels.size());
         ++finest_level)
      results.push_back(run_cycle_case(levels,
                                       transfers,
                                       downstream_orders,
                                       finest_level,
                                       epsilon));

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack convection-diffusion SUPG geometric V-cycle diagnostic\n";
    std::cout << "book-style exponential outflow layer, beta=(1,0), epsilon="
              << epsilon << '\n';
    std::cout << std::setw(8) << "levels" << std::setw(14)
              << "h" << std::setw(10) << "dofs" << std::setw(8)
              << "cycles" << std::setw(8) << "smooth" << std::setw(14)
              << "r0" << std::setw(14) << "rJacobi" << std::setw(14)
              << "rGS" << std::setw(14) << "rV" << std::setw(14)
              << "V factor" << std::setw(14) << "L2" << std::setw(14)
              << "H1-semi" << '\n';
    for (const auto &result : results)
      {
        std::cout << std::setw(8) << result.levels << std::setw(14)
                  << result.finest_h << std::setw(10)
                  << result.finest_dofs << std::setw(8)
                  << result.cycles << std::setw(8)
                  << result.smoothing_sweeps << std::setw(14)
                  << result.initial_residual << std::setw(14)
                  << result.jacobi_residual << std::setw(14)
                  << result.gs_residual << std::setw(14)
                  << result.cycle_residual << std::setw(14)
                  << result.cycle_factor << std::setw(14)
                  << result.l2_error << std::setw(14)
                  << result.h1_error << '\n';
      }
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [epsilon]\n";
      return 2;
    }

  try
    {
      const double epsilon = argc > 1 ? std::stod(argv[1]) : 1.0e-3;
      if (!(epsilon > 0.0))
        {
          std::cerr << "epsilon must be positive\n";
          return 2;
        }

      run_default_study(epsilon);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
