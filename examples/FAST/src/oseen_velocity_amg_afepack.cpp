/**
 * @file oseen_velocity_amg_afepack.cpp
 * @brief AFEPack Oseen 线性化方程算例：速度块 AMG。
 *
 * @details
 * 本文件是 FAST 目录中的 Oseen 线性化方程 迁移算例，关注Navier-Stokes 线性化后得到的鞍点系统及其 Schur 补预条件器。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 速度块 AMG：展示速度子问题的代数多重网格预处理在鞍点系统中的使用方式。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define main stokes_velocity_amg_default_main
#include "stokes_velocity_amg_afepack.cpp"
#undef main

#include <functional>

namespace
{
  struct OseenComponentResult
  {
    unsigned int nonzeros = 0;
    SolveInfo    gmres;
    SolveInfo    jacobi_gmres;
    SolveInfo    amg_gmres;
    double       amg_l2_error = 0.0;
    double       amg_h1_error = 0.0;
  };

  struct OseenVelocityResult
  {
    double       h = 0.0;
    int          cells = 0;
    unsigned int scalar_dofs = 0;
    unsigned int velocity_dofs = 0;
    unsigned int scalar_nonzeros = 0;
    int          gmres_max_iterations = 0;
    int          gmres_total_iterations = 0;
    int          jacobi_max_iterations = 0;
    int          jacobi_total_iterations = 0;
    int          amg_max_iterations = 0;
    int          amg_total_iterations = 0;
    double       gmres_relative_residual = 0.0;
    double       jacobi_relative_residual = 0.0;
    double       amg_relative_residual = 0.0;
    double       velocity_l2_error = 0.0;
    double       velocity_h1_error = 0.0;
    bool         has_rates = false;
    double       velocity_l2_rate = 0.0;
    double       velocity_h1_rate = 0.0;
  };

  double current_oseen_viscosity = 1.0;
  double current_wind_scale = 1.0;

  double wind_dot_gradient(const double *p,
                           const GradientFunction gradient_function)
  {
    const auto gradient = gradient_function(p);
    return velocity_x_exact(p) * gradient[0] +
           velocity_y_exact(p) * gradient[1];
  }

  double oseen_velocity_x_rhs(const double *p)
  {
    return -current_oseen_viscosity * laplace_velocity_x(p) +
           current_wind_scale * wind_dot_gradient(p, &velocity_x_gradient);
  }

  double oseen_velocity_y_rhs(const double *p)
  {
    return -current_oseen_viscosity * laplace_velocity_y(p) +
           current_wind_scale * wind_dot_gradient(p, &velocity_y_gradient);
  }

  class OseenVelocityMatrix : public BilinearOperator<dim, double>
  {
  public:
    explicit OseenVelocityMatrix(FEMSpace<double, dim> &space)
      : BilinearOperator<dim, double>(space, space)
    {}

    void getElementMatrix(
      const Element<double, dim> &ele0,
      const Element<double, dim> &,
      const ActiveElementPairIterator<dim>::State =
        ActiveElementPairIterator<dim>::EQUAL) override
    {
      const double volume = ele0.templateElement().volume();
      const QuadratureInfo<dim> &quad_info =
        ele0.findQuadratureInfo(algebricAccuracy());
      const auto jacobian =
        ele0.local_to_global_jacobian(quad_info.quadraturePoint());
      const auto q_point = ele0.local_to_global(quad_info.quadraturePoint());
      const auto basis_value = ele0.basis_function_value(q_point);
      const auto basis_grad = ele0.basis_function_gradient(q_point);

      const int n_q = quad_info.n_quadraturePoint();
      const int n_dof = ele0.dof().size();
      for (int q = 0; q < n_q; ++q)
        {
          const double p[dim] = {q_point[q][0], q_point[q][1]};
          const double wind[dim] = {velocity_x_exact(p),
                                    velocity_y_exact(p)};
          const double jxw = quad_info.weight(q) * jacobian[q] * volume;
          for (int i = 0; i < n_dof; ++i)
            for (int j = 0; j < n_dof; ++j)
              {
                const double diffusion =
                  current_oseen_viscosity *
                  innerProduct(basis_grad[i][q],
                               basis_grad[j][q]);
                const double convection =
                  current_wind_scale *
                  (wind[0] * basis_grad[j][q][0] +
                   wind[1] * basis_grad[j][q][1]) *
                  basis_value[i][q];
                elementMatrix(i, j) += jxw * (diffusion + convection);
              }
        }
    }
  };

  using PreconditionerAction =
    std::function<void(const Vector<double> &, Vector<double> &)>;

  void identity_apply(const Vector<double> &src, Vector<double> &dst)
  {
    dst.reinit(src.size(), false);
    dst = src;
  }

  void jacobi_apply_local(const SparseMatrix<double> &matrix,
                          const Vector<double> &src,
                          Vector<double> &dst)
  {
    dst.reinit(src.size(), false);
    for (std::size_t i = 0; i < src.size(); ++i)
      {
        const double diagonal = matrix.diag_element(i);
        dst(i) = std::abs(diagonal) > 1.0e-14 ? src(i) / diagonal : src(i);
      }
  }

  Vector<double> residual_local(const SparseMatrix<double> &matrix,
                                const Vector<double> &x,
                                const Vector<double> &rhs)
  {
    Vector<double> ax(rhs.size());
    matrix.vmult(ax, x);

    Vector<double> r = rhs;
    r.add(-1.0, ax);
    return r;
  }

  void normalize_local(Vector<double> &v, const double norm)
  {
    if (norm == 0.0)
      return;
    v.scale(1.0 / norm);
  }

  void add_scaled_local(Vector<double> &destination,
                        const double coefficient,
                        const Vector<double> &source)
  {
    destination.add(coefficient, source);
  }

  void generate_plane_rotation_local(const double dx,
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
        const double temp = dx / dy;
        sn = 1.0 / std::sqrt(1.0 + temp * temp);
        cs = temp * sn;
      }
    else
      {
        const double temp = dy / dx;
        cs = 1.0 / std::sqrt(1.0 + temp * temp);
        sn = temp * cs;
      }
  }

  void apply_plane_rotation_local(double &dx,
                                  double &dy,
                                  const double cs,
                                  const double sn)
  {
    const double temp = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
  }

  std::vector<double> solve_upper_hessenberg_local(
    const std::vector<std::vector<double>> &h,
    const std::vector<double> &g,
    const int n)
  {
    std::vector<double> y(n, 0.0);
    for (int i = n - 1; i >= 0; --i)
      {
        double rhs = g[i];
        for (int j = i + 1; j < n; ++j)
          rhs -= h[i][j] * y[j];
        y[i] = rhs / h[i][i];
      }
    return y;
  }

  SolveInfo solve_right_preconditioned_gmres(
    const SparseMatrix<double> &matrix,
    Vector<double> &x,
    const Vector<double> &rhs,
    const PreconditionerAction &preconditioner,
    const double tolerance,
    const int restart,
    const int max_iterations)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);
    Vector<double> r = residual_local(matrix, x, rhs);
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
        normalize_local(v[0], beta);
        g[0] = beta;

        int inner_iterations = 0;
        for (; inner_iterations < restart &&
               total_iterations < max_iterations;
             ++inner_iterations, ++total_iterations)
          {
            preconditioner(v[inner_iterations], z[inner_iterations]);
            matrix.vmult(v[inner_iterations + 1], z[inner_iterations]);

            for (int i = 0; i <= inner_iterations; ++i)
              {
                h[i][inner_iterations] =
                  v[inner_iterations + 1].dot(v[i]);
                add_scaled_local(v[inner_iterations + 1],
                                 -h[i][inner_iterations],
                                 v[i]);
              }

            h[inner_iterations + 1][inner_iterations] =
              v[inner_iterations + 1].l2_norm();
            if (h[inner_iterations + 1][inner_iterations] != 0.0)
              normalize_local(v[inner_iterations + 1],
                              h[inner_iterations + 1][inner_iterations]);

            for (int i = 0; i < inner_iterations; ++i)
              apply_plane_rotation_local(h[i][inner_iterations],
                                         h[i + 1][inner_iterations],
                                         cs[i],
                                         sn[i]);

            generate_plane_rotation_local(
              h[inner_iterations][inner_iterations],
              h[inner_iterations + 1][inner_iterations],
              cs[inner_iterations],
              sn[inner_iterations]);
            apply_plane_rotation_local(
              h[inner_iterations][inner_iterations],
              h[inner_iterations + 1][inner_iterations],
              cs[inner_iterations],
              sn[inner_iterations]);
            apply_plane_rotation_local(g[inner_iterations],
                                       g[inner_iterations + 1],
                                       cs[inner_iterations],
                                       sn[inner_iterations]);

            relative_residual =
              std::abs(g[inner_iterations + 1]) / rhs_norm;
            if (relative_residual <= tolerance)
              {
                const int system_size = inner_iterations + 1;
                const auto y =
                  solve_upper_hessenberg_local(h, g, system_size);
                for (int i = 0; i < system_size; ++i)
                  add_scaled_local(x, y[i], z[i]);
                return {true, total_iterations + 1, relative_residual};
              }
          }

        const auto y =
          solve_upper_hessenberg_local(h, g, inner_iterations);
        for (int i = 0; i < inner_iterations; ++i)
          add_scaled_local(x, y[i], z[i]);

        r = residual_local(matrix, x, rhs);
        beta = r.l2_norm();
        relative_residual = beta / rhs_norm;
        if (relative_residual <= tolerance)
          return {true, total_iterations, relative_residual};
      }

    return {false, max_iterations, relative_residual};
  }

  void apply_homogeneous_oseen_boundary(
    OseenVelocityMatrix &matrix,
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

  OseenComponentResult solve_oseen_component(
    FEMSpace<double, dim> &space,
    ScalarFunction rhs_function,
    ScalarFunction exact_function,
    GradientFunction exact_gradient,
    const bool write_solution,
    const std::string &output_name)
  {
    OseenVelocityMatrix oseen_matrix(space);
    oseen_matrix.algebricAccuracy() = 5;
    oseen_matrix.build();

    StiffMatrix<dim, double> diffusion_matrix(space);
    diffusion_matrix.algebricAccuracy() = 5;
    diffusion_matrix.build();

    Vector<double> rhs;
    Operator::L2Discretize(rhs_function, space, rhs, 5);

    FEMFunction<double, dim> oseen_boundary_solution(space);
    apply_homogeneous_oseen_boundary(oseen_matrix,
                                     oseen_boundary_solution,
                                     rhs,
                                     space);

    Vector<double> diffusion_dummy_rhs(space.n_dof());
    diffusion_dummy_rhs = 0.0;
    FEMFunction<double, dim> diffusion_boundary_solution(space);
    apply_homogeneous_dirichlet_boundary(diffusion_matrix,
                                         diffusion_boundary_solution,
                                         diffusion_dummy_rhs,
                                         space);

    FEMFunction<double, dim> gmres_solution(space);
    FEMFunction<double, dim> jacobi_solution(space);
    FEMFunction<double, dim> amg_solution(space);

    const SolveInfo gmres_info =
      solve_right_preconditioned_gmres(
        oseen_matrix,
        gmres_solution,
        rhs,
        &identity_apply,
        1.0e-8,
        60,
        4000);
    if (!gmres_info.converged)
      throw std::runtime_error("Oseen velocity GMRES did not converge");

    const SolveInfo jacobi_info =
      solve_right_preconditioned_gmres(
        oseen_matrix,
        jacobi_solution,
        rhs,
        [&oseen_matrix](const Vector<double> &src, Vector<double> &dst) {
          jacobi_apply_local(oseen_matrix, src, dst);
        },
        1.0e-8,
        60,
        4000);
    if (!jacobi_info.converged)
      throw std::runtime_error("Oseen velocity Jacobi-GMRES did not converge");

    AMGPreconditioner amg_preconditioner(diffusion_matrix, 2, 2);
    amg_preconditioner.smoothStep() = 2;
    const SolveInfo amg_info =
      solve_right_preconditioned_gmres(
        oseen_matrix,
        amg_solution,
        rhs,
        [&amg_preconditioner](const Vector<double> &src, Vector<double> &dst) {
          dst.reinit(src.size(), false);
          dst = 0.0;
          amg_preconditioner.vmult(dst, src);
          dst.scale(1.0 / current_oseen_viscosity);
        },
        1.0e-8,
        60,
        4000);
    if (!amg_info.converged)
      throw std::runtime_error(
        "Oseen velocity AMG-GMRES did not converge");

    if (write_solution)
      amg_solution.writeOpenDXData(std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                                   "/" + output_name);

    FunctionFunction<double> exact(exact_function, exact_gradient);
    return {static_cast<unsigned int>(oseen_matrix.n_nonzero_elements()),
            gmres_info,
            jacobi_info,
            amg_info,
            Functional::L2Error(amg_solution, exact, 5),
            Functional::H1SemiError(amg_solution, exact, 5)};
  }

  OseenVelocityResult run_oseen_mesh_case(const double h,
                                          const std::string &mesh_file,
                                          const bool write_solution)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof(geometry);
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);
    auto velocity_space = build_space(mesh, template_element);

    const OseenComponentResult ux =
      solve_oseen_component(velocity_space,
                            &oseen_velocity_x_rhs,
                            &velocity_x_exact,
                            &velocity_x_gradient,
                            write_solution,
                            "oseen_velocity_amg_ux_h0p05.dx");
    const OseenComponentResult uy =
      solve_oseen_component(velocity_space,
                            &oseen_velocity_y_rhs,
                            &velocity_y_exact,
                            &velocity_y_gradient,
                            write_solution,
                            "oseen_velocity_amg_uy_h0p05.dx");

    return {h,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            2u * velocity_space.n_dof(),
            std::max(ux.nonzeros, uy.nonzeros),
            std::max(ux.gmres.iterations, uy.gmres.iterations),
            ux.gmres.iterations + uy.gmres.iterations,
            std::max(ux.jacobi_gmres.iterations,
                     uy.jacobi_gmres.iterations),
            ux.jacobi_gmres.iterations + uy.jacobi_gmres.iterations,
            std::max(ux.amg_gmres.iterations, uy.amg_gmres.iterations),
            ux.amg_gmres.iterations + uy.amg_gmres.iterations,
            std::max(ux.gmres.relative_residual,
                     uy.gmres.relative_residual),
            std::max(ux.jacobi_gmres.relative_residual,
                     uy.jacobi_gmres.relative_residual),
            std::max(ux.amg_gmres.relative_residual,
                     uy.amg_gmres.relative_residual),
            std::sqrt(ux.amg_l2_error * ux.amg_l2_error +
                      uy.amg_l2_error * uy.amg_l2_error),
            std::sqrt(ux.amg_h1_error * ux.amg_h1_error +
                      uy.amg_h1_error * uy.amg_h1_error)};
  }

  void add_oseen_rates(std::vector<OseenVelocityResult> &results)
  {
    for (std::size_t i = 1; i < results.size(); ++i)
      {
        results[i].has_rates = true;
        results[i].velocity_l2_rate =
          convergence_rate(results[i - 1].h,
                           results[i].h,
                           results[i - 1].velocity_l2_error,
                           results[i].velocity_l2_error);
        results[i].velocity_h1_rate =
          convergence_rate(results[i - 1].h,
                           results[i].h,
                           results[i - 1].velocity_h1_error,
                           results[i].velocity_h1_error);
      }
  }

  void run_default_oseen_velocity_study(const std::string &mesh_dir)
  {
    current_oseen_viscosity = 1.0;
    current_wind_scale = 1.0;

    const std::vector<std::pair<double, std::string>> mesh_cases = {
      {0.20, mesh_dir + "/unit_square_h0p20"},
      {0.10, mesh_dir + "/unit_square_h0p10"},
      {0.05, mesh_dir + "/unit_square_h0p05"}};

    std::vector<OseenVelocityResult> results;
    for (std::size_t i = 0; i < mesh_cases.size(); ++i)
      results.push_back(run_oseen_mesh_case(mesh_cases[i].first,
                                            mesh_cases[i].second,
                                            i + 1 == mesh_cases.size()));

    add_oseen_rates(results);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Oseen velocity block AMG diagnostic\n";
    std::cout << "P2 scalar velocity blocks, frozen manufactured wind, "
                 "homogeneous Dirichlet.\n";
    std::cout << "h cells scalar_dofs velocity_dofs scalar_nnz "
                 "GMRESmax GMRESsum JGMRESmax JGMRESsum "
                 "AMG-GMRESmax AMG-GMRESsum GMRES_res JGMRES_res "
                 "AMG_res velocity_L2 velocity_H1 "
                 "velocity_L2_rate velocity_H1_rate\n";

    for (const auto &result : results)
      std::cout << result.h << ' '
                << result.cells << ' '
                << result.scalar_dofs << ' '
                << result.velocity_dofs << ' '
                << result.scalar_nonzeros << ' '
                << result.gmres_max_iterations << ' '
                << result.gmres_total_iterations << ' '
                << result.jacobi_max_iterations << ' '
                << result.jacobi_total_iterations << ' '
                << result.amg_max_iterations << ' '
                << result.amg_total_iterations << ' '
                << result.gmres_relative_residual << ' '
                << result.jacobi_relative_residual << ' '
                << result.amg_relative_residual << ' '
                << result.velocity_l2_error << ' '
                << result.velocity_h1_error << ' '
                << std::setw(9)
                << format_rate(result.has_rates,
                               result.velocity_l2_rate) << ' '
                << std::setw(9)
                << format_rate(result.has_rates,
                               result.velocity_h1_rate) << '\n';

    if (!(results.back().velocity_l2_error < results.front().velocity_l2_error &&
          results.back().velocity_h1_error < results.front().velocity_h1_error))
      throw std::runtime_error("Oseen velocity errors did not decrease");

    for (const auto &result : results)
      if (result.amg_max_iterations > result.gmres_max_iterations)
        throw std::runtime_error(
          "AMG-preconditioned Oseen solve should improve over GMRES");

    std::cout << "Wrote "
              << std::string(AFEPACK_DEFAULT_OUTPUT_DIR)
              << "/oseen_velocity_amg_{ux,uy}_h0p05.dx\n";
  }

  double mesh_h_from_base(const std::string &mesh_base)
  {
    if (mesh_base.find("h0p20") != std::string::npos)
      return 0.20;
    if (mesh_base.find("h0p10") != std::string::npos)
      return 0.10;
    if (mesh_base.find("h0p05") != std::string::npos)
      return 0.05;
    if (mesh_base.find("h0p025") != std::string::npos)
      return 0.025;
    return 0.0;
  }

  std::vector<double> parse_viscosity_list(const std::string &csv)
  {
    std::vector<double> values;
    std::stringstream input(csv);
    std::string token;
    while (std::getline(input, token, ','))
      {
        if (token.empty())
          continue;
        const double value = std::stod(token);
        if (value <= 0.0)
          throw std::runtime_error("viscosities must be positive");
        values.push_back(value);
      }

    if (values.empty())
      throw std::runtime_error("empty viscosity list");
    return values;
  }

  void run_oseen_velocity_viscosity_sweep(const std::string &mesh_base,
                                          const std::vector<double> &viscosities)
  {
    current_wind_scale = 1.0;
    const double h = mesh_h_from_base(mesh_base);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Oseen velocity block viscosity sweep\n";
    std::cout << "mesh nu formal_Re cells scalar_dofs GMRESmax JGMRESmax "
                 "AMG-GMRESmax GMRES_res JGMRES_res AMG_res "
                 "velocity_L2 velocity_H1\n";

    for (const double nu : viscosities)
      {
        current_oseen_viscosity = nu;
        const OseenVelocityResult result =
          run_oseen_mesh_case(h, mesh_base, false);

        std::cout << mesh_base << ' '
                  << nu << ' '
                  << 1.0 / nu << ' '
                  << result.cells << ' '
                  << result.scalar_dofs << ' '
                  << result.gmres_max_iterations << ' '
                  << result.jacobi_max_iterations << ' '
                  << result.amg_max_iterations << ' '
                  << result.gmres_relative_residual << ' '
                  << result.jacobi_relative_residual << ' '
                  << result.amg_relative_residual << ' '
                  << result.velocity_l2_error << ' '
                  << result.velocity_h1_error << '\n';
      }
  }
}

int main(int argc, char **argv)
{
  if (argc > 4)
    {
      std::cerr << "usage: " << argv[0]
                << " [mesh-directory]\n"
                << "       " << argv[0]
                << " --viscosity-sweep [mesh-base] [nu_csv]\n";
      return 2;
    }

  try
    {
      if (argc > 1 && std::string(argv[1]) == "--viscosity-sweep")
        {
          const std::string mesh_base =
            argc > 2 ? argv[2]
                     : std::string(AFEPACK_DEFAULT_MESH_DIR) +
                         "/unit_square_h0p10";
          const std::string viscosity_csv =
            argc > 3 ? argv[3] : "1,0.5,0.25,0.125,0.0625";
          run_oseen_velocity_viscosity_sweep(
            mesh_base, parse_viscosity_list(viscosity_csv));
        }
      else
        {
          const std::string mesh_dir =
            argc > 1 ? argv[1] : std::string(AFEPACK_DEFAULT_MESH_DIR);
          run_default_oseen_velocity_study(mesh_dir);
        }
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
