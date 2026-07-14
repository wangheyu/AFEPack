/**
 * @file poisson_book_examples_afepack.cpp
 * @brief AFEPack Poisson 方程算例：教材例题合集。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 教材例题合集：集中组织教材中的 Poisson 基础算例，用统一入口生成多组对照结果。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
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
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
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

  enum class BookProblem
  {
    square_constant_source,
    lshape_constant_source,
    square_harmonic,
    lshape_singular
  };

  struct ErrorNorms
  {
    double l2 = -1.0;
    double h1_seminorm = -1.0;
  };

  struct PointDiagnostic
  {
    double value = 0.0;
    double distance = 0.0;
  };

  struct PoissonBookResult
  {
    BookProblem  problem;
    std::string  mesh_file;
    std::string  output_file;
    double       h = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    PointDiagnostic center;
    double       max_nodal_value = 0.0;
    ErrorNorms   error;
  };

  BookProblem active_problem = BookProblem::square_constant_source;

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

  void book_coordinates(const double *p, double &x, double &y)
  {
    x = 2.0 * p[0] - 1.0;
    y = 2.0 * p[1] - 1.0;
  }

  double square_harmonic_exact_from_book_coordinates(const double x,
                                                     const double y)
  {
    const double numerator = 2.0 * (1.0 + y);
    const double denominator =
      (3.0 + x) * (3.0 + x) + (1.0 + y) * (1.0 + y);
    return numerator / denominator;
  }

  std::vector<double>
  square_harmonic_gradient_in_unit_coordinates(const double *p)
  {
    double x = 0.0;
    double y = 0.0;
    book_coordinates(p, x, y);

    const double a = 3.0 + x;
    const double b = 1.0 + y;
    const double denominator = a * a + b * b;
    const double denominator_square = denominator * denominator;

    const double du_dx_book = -4.0 * a * b / denominator_square;
    const double du_dy_book =
      2.0 / denominator - 4.0 * b * b / denominator_square;

    return {2.0 * du_dx_book, 2.0 * du_dy_book};
  }

  double lshape_singular_exact_from_book_coordinates(const double x,
                                                     const double y)
  {
    const double radius = std::hypot(x, y);
    if (radius < 1.0e-14)
      return 0.0;

    const double theta = std::atan2(y, x);
    return std::pow(radius, 2.0 / 3.0) *
           std::sin((2.0 * theta + pi) / 3.0);
  }

  std::vector<double>
  lshape_singular_gradient_in_unit_coordinates(const double *p)
  {
    double x = 0.0;
    double y = 0.0;
    book_coordinates(p, x, y);

    const double radius = std::hypot(x, y);
    if (radius < 1.0e-14)
      return {0.0, 0.0};

    const double alpha = 2.0 / 3.0;
    const double theta = std::atan2(y, x);
    const double psi = (2.0 * theta + pi) / 3.0;
    const double scale = alpha * std::pow(radius, alpha - 1.0);

    const double du_dx_book =
      scale * (std::sin(psi) * std::cos(theta) -
               std::cos(psi) * std::sin(theta));
    const double du_dy_book =
      scale * (std::sin(psi) * std::sin(theta) +
               std::cos(psi) * std::cos(theta));

    return {2.0 * du_dx_book, 2.0 * du_dy_book};
  }

  double exact_solution(const double *p)
  {
    double x = 0.0;
    double y = 0.0;
    book_coordinates(p, x, y);

    switch (active_problem)
      {
        case BookProblem::square_constant_source:
        case BookProblem::lshape_constant_source:
          return 0.0;
        case BookProblem::square_harmonic:
          return square_harmonic_exact_from_book_coordinates(x, y);
        case BookProblem::lshape_singular:
          return lshape_singular_exact_from_book_coordinates(x, y);
      }

    return 0.0;
  }

  std::vector<double> exact_gradient(const double *p)
  {
    switch (active_problem)
      {
        case BookProblem::square_constant_source:
        case BookProblem::lshape_constant_source:
          return {0.0, 0.0};
        case BookProblem::square_harmonic:
          return square_harmonic_gradient_in_unit_coordinates(p);
        case BookProblem::lshape_singular:
          return lshape_singular_gradient_in_unit_coordinates(p);
      }

    return {0.0, 0.0};
  }

  double right_hand_side(const double *p)
  {
    (void)p;
    switch (active_problem)
      {
        case BookProblem::square_constant_source:
        case BookProblem::lshape_constant_source:
          return 4.0;
        case BookProblem::square_harmonic:
        case BookProblem::lshape_singular:
          return 0.0;
      }

    return 0.0;
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
    BoundaryFunction<double, dim> boundary5(
      BoundaryConditionInfo::DIRICHLET, 5, &exact_solution);
    BoundaryFunction<double, dim> boundary6(
      BoundaryConditionInfo::DIRICHLET, 6, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.add(boundary5);
    boundary_admin.add(boundary6);
    boundary_admin.apply(stiff_matrix, solution, rhs);
  }

  PointDiagnostic nearest_dof_value(const FEMFunction<double, dim> &solution,
                                    const double x,
                                    const double y)
  {
    const auto &space = solution.femSpace();

    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto &element = space.element(e);
        const auto &dofs = element.dof();
        if (dofs.size() != 3)
          continue;

        const auto &p0 = space.dofInfo(dofs[0]).interp_point;
        const auto &p1 = space.dofInfo(dofs[1]).interp_point;
        const auto &p2 = space.dofInfo(dofs[2]).interp_point;

        const double denominator =
          (p1[1] - p2[1]) * (p0[0] - p2[0]) +
          (p2[0] - p1[0]) * (p0[1] - p2[1]);
        if (std::abs(denominator) < 1.0e-14)
          continue;

        const double lambda0 =
          ((p1[1] - p2[1]) * (x - p2[0]) +
           (p2[0] - p1[0]) * (y - p2[1])) /
          denominator;
        const double lambda1 =
          ((p2[1] - p0[1]) * (x - p2[0]) +
           (p0[0] - p2[0]) * (y - p2[1])) /
          denominator;
        const double lambda2 = 1.0 - lambda0 - lambda1;

        if (lambda0 >= -1.0e-12 && lambda1 >= -1.0e-12 &&
            lambda2 >= -1.0e-12)
          {
            const double value =
              lambda0 * solution(dofs[0]) + lambda1 * solution(dofs[1]) +
              lambda2 * solution(dofs[2]);
            return {value, 0.0};
          }
      }

    unsigned int nearest = 0;
    double nearest_distance_square = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < solution.size(); ++i)
      {
        const auto &point = space.dofInfo(static_cast<int>(i)).interp_point;
        const double dx = point[0] - x;
        const double dy = point[1] - y;
        const double distance_square = dx * dx + dy * dy;
        if (distance_square < nearest_distance_square)
          {
            nearest = i;
            nearest_distance_square = distance_square;
          }
      }

    return {solution(nearest), std::sqrt(nearest_distance_square)};
  }

  double max_nodal_value(const FEMFunction<double, dim> &solution)
  {
    double maximum = -std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < solution.size(); ++i)
      maximum = std::max(maximum, static_cast<double>(solution(i)));
    return maximum;
  }

  bool has_exact_error(const BookProblem problem)
  {
    return problem == BookProblem::square_harmonic ||
           problem == BookProblem::lshape_singular;
  }

  ErrorNorms compute_error_norms(const FEMFunction<double, dim> &solution)
  {
    double l2_square = 0.0;
    double h1_square = 0.0;

    const auto &space = solution.femSpace();
    for (int e = 0; e < space.n_element(); ++e)
      {
        const auto &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info = element.findQuadratureInfo(5);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto values = solution.value(q_point, element);
        const auto gradients = solution.gradient(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double point[dim] = {q_point[q][0], q_point[q][1]};
            const double value_error = exact_solution(point) - values[q];
            const auto exact_grad = exact_gradient(point);
            const double grad_error_x = exact_grad[0] - gradients[q][0];
            const double grad_error_y = exact_grad[1] - gradients[q][1];
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;

            l2_square += jxw * value_error * value_error;
            h1_square +=
              jxw * (grad_error_x * grad_error_x +
                     grad_error_y * grad_error_y);
          }
      }

    return {2.0 * std::sqrt(l2_square), std::sqrt(h1_square)};
  }

  PoissonBookResult solve_problem(const BookProblem problem,
                                  const std::string &mesh_file,
                                  const std::string &output_file,
                                  const double h,
                                  const bool write_solution)
  {
    active_problem = problem;

    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> geometry;
    CoordTransform<dim, dim> transform;
    TemplateDOF<dim> dof(geometry);
    BasisFunctionAdmin<double, dim, dim> basis(dof);
    auto template_element =
      make_triangle_template(geometry, transform, dof, basis);
    auto fem_space = build_space(mesh, template_element);

    StiffMatrix<dim, double> stiff_matrix(fem_space);
    stiff_matrix.algebricAccuracy() = 4;
    stiff_matrix.build();

    FEMFunction<double, dim> solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    apply_dirichlet_boundary(stiff_matrix, solution, rhs, fem_space);

    AMGSolver solver(stiff_matrix);
    solver.solve(solution, rhs, 1.0e-10, 500);

    if (write_solution)
      solution.writeOpenDXData(output_file);

    PoissonBookResult result;
    result.problem = problem;
    result.mesh_file = mesh_file;
    result.output_file = write_solution ? output_file : "";
    result.h = h;
    result.cells = mesh.n_geometry(dim);
    result.dofs = fem_space.n_dof();
    result.center = nearest_dof_value(solution, 0.5, 0.5);
    result.max_nodal_value = max_nodal_value(solution);
    if (has_exact_error(problem))
      result.error = compute_error_norms(solution);

    return result;
  }

  double convergence_rate(const double coarse_error,
                          const double fine_error,
                          const double coarse_h,
                          const double fine_h)
  {
    if (coarse_error <= 0.0 || fine_error <= 0.0)
      return std::numeric_limits<double>::quiet_NaN();
    return std::log(coarse_error / fine_error) /
           std::log(coarse_h / fine_h);
  }

  std::string mesh_dir()
  {
    return AFEPACK_DEFAULT_MESH_DIR;
  }

  std::string output_dir()
  {
    return AFEPACK_DEFAULT_OUTPUT_DIR;
  }

  std::vector<std::pair<double, std::string>> square_meshes()
  {
    const std::string dir = mesh_dir();
    return {{0.20, dir + "/unit_square_h0p20"},
            {0.10, dir + "/unit_square_h0p10"},
            {0.05, dir + "/unit_square_h0p05"},
            {0.025, dir + "/unit_square_h0p025"}};
  }

  std::vector<std::pair<double, std::string>> lshape_meshes()
  {
    const std::string dir = mesh_dir();
    return {{0.20, dir + "/lshape_h0p20"},
            {0.10, dir + "/lshape_h0p10"},
            {0.05, dir + "/lshape_h0p05"}};
  }

  void print_square_constant_source()
  {
    constexpr double exact_center = 0.294685413126;
    std::vector<PoissonBookResult> results;
    const auto meshes = square_meshes();
    for (std::size_t i = 0; i < meshes.size(); ++i)
      {
        const bool write_solution = i + 1 == meshes.size();
        results.push_back(solve_problem(
          BookProblem::square_constant_source,
          meshes[i].second,
          output_dir() + "/poisson_book_example_1_1_1.dx",
          meshes[i].first,
          write_solution));
      }

    std::cout << "Example 1.1.1 square, f=1, zero boundary\n";
    std::cout << std::setw(10) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(16)
              << "uh(0,0)" << std::setw(16) << "abs error"
              << std::setw(14) << "rate" << std::setw(14)
              << "dof dist" << '\n';

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        const double error =
          std::abs(results[i].center.value - exact_center);
        std::cout << std::setw(10) << results[i].h
                  << std::setw(10) << results[i].cells
                  << std::setw(10) << results[i].dofs
                  << std::setw(16) << results[i].center.value
                  << std::setw(16) << error;
        if (i == 0)
          std::cout << std::setw(14) << "-";
        else
          std::cout << std::setw(14)
                    << convergence_rate(
                         std::abs(results[i - 1].center.value -
                                  exact_center),
                         error,
                         results[i - 1].h,
                         results[i].h);
        std::cout << std::setw(14) << results[i].center.distance << '\n';
      }
    std::cout << "Exact u(0,0): " << exact_center << '\n';
    std::cout << "Wrote " << results.back().output_file << "\n\n";
  }

  void print_lshape_constant_source()
  {
    std::vector<PoissonBookResult> results;
    const auto meshes = lshape_meshes();
    for (std::size_t i = 0; i < meshes.size(); ++i)
      {
        const bool write_solution = i + 1 == meshes.size();
        results.push_back(solve_problem(
          BookProblem::lshape_constant_source,
          meshes[i].second,
          output_dir() + "/poisson_book_example_1_1_2.dx",
          meshes[i].first,
          write_solution));
      }

    std::cout << "Example 1.1.2 L-shaped domain, f=1, zero boundary\n";
    std::cout << std::setw(10) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(16)
              << "max uh" << std::setw(16) << "successive"
              << std::setw(14) << "rate" << '\n';

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        std::cout << std::setw(10) << results[i].h
                  << std::setw(10) << results[i].cells
                  << std::setw(10) << results[i].dofs
                  << std::setw(16) << results[i].max_nodal_value;
        if (i == 0)
          {
            std::cout << std::setw(16) << "-" << std::setw(14) << "-";
          }
        else
          {
            const double current_difference =
              std::abs(results[i].max_nodal_value -
                       results[i - 1].max_nodal_value);
            std::cout << std::setw(16) << current_difference;
            if (i == 1)
              std::cout << std::setw(14) << "-";
            else
              {
                const double previous_difference =
                  std::abs(results[i - 1].max_nodal_value -
                           results[i - 2].max_nodal_value);
                std::cout << std::setw(14)
                          << convergence_rate(previous_difference,
                                              current_difference,
                                              results[i - 2].h,
                                              results[i - 1].h);
              }
          }
        std::cout << '\n';
      }
    std::cout << "Wrote " << results.back().output_file << "\n\n";
  }

  void print_square_harmonic()
  {
    std::vector<PoissonBookResult> results;
    const auto meshes = square_meshes();
    for (std::size_t i = 0; i < meshes.size(); ++i)
      {
        const bool write_solution = i + 1 == meshes.size();
        results.push_back(solve_problem(
          BookProblem::square_harmonic,
          meshes[i].second,
          output_dir() + "/poisson_book_example_1_1_3.dx",
          meshes[i].first,
          write_solution));
      }

    std::cout << "Example 1.1.3 square harmonic analytic solution\n";
    std::cout << std::setw(10) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(16)
              << "L2" << std::setw(16) << "H1 semi"
              << std::setw(14) << "H1 rate" << '\n';

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        std::cout << std::setw(10) << results[i].h
                  << std::setw(10) << results[i].cells
                  << std::setw(10) << results[i].dofs
                  << std::setw(16) << results[i].error.l2
                  << std::setw(16) << results[i].error.h1_seminorm;
        if (i == 0)
          std::cout << std::setw(14) << "-";
        else
          std::cout << std::setw(14)
                    << convergence_rate(results[i - 1].error.h1_seminorm,
                                        results[i].error.h1_seminorm,
                                        results[i - 1].h,
                                        results[i].h);
        std::cout << '\n';
      }
    std::cout << "Wrote " << results.back().output_file << "\n\n";
  }

  void print_lshape_singular()
  {
    std::vector<PoissonBookResult> results;
    const auto meshes = lshape_meshes();
    for (std::size_t i = 0; i < meshes.size(); ++i)
      {
        const bool write_solution = i + 1 == meshes.size();
        results.push_back(solve_problem(
          BookProblem::lshape_singular,
          meshes[i].second,
          output_dir() + "/poisson_book_example_1_1_4.dx",
          meshes[i].first,
          write_solution));
      }

    std::cout << "Example 1.1.4 L-shaped singular analytic solution\n";
    std::cout << std::setw(10) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(16)
              << "L2" << std::setw(16) << "H1 semi"
              << std::setw(14) << "H1 rate" << '\n';

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        std::cout << std::setw(10) << results[i].h
                  << std::setw(10) << results[i].cells
                  << std::setw(10) << results[i].dofs
                  << std::setw(16) << results[i].error.l2
                  << std::setw(16) << results[i].error.h1_seminorm;
        if (i == 0)
          std::cout << std::setw(14) << "-";
        else
          std::cout << std::setw(14)
                    << convergence_rate(results[i - 1].error.h1_seminorm,
                                        results[i].error.h1_seminorm,
                                        results[i - 1].h,
                                        results[i].h);
        std::cout << '\n';
      }
    std::cout << "Expected H1 rate is limited by the r^(2/3) singularity.\n";
    std::cout << "Wrote " << results.back().output_file << "\n\n";
  }
}

int main()
{
  try
    {
      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Poisson book reference examples 1.1.1-1.1.4\n";
      std::cout << "Book coordinates (X,Y)=2(x,y)-(1,1); RHS is scaled by 4.\n\n";

      print_square_constant_source();
      print_lshape_constant_source();
      print_square_harmonic();
      print_lshape_singular();
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
