/**
 * @file poisson_mixed_boundary_afepack.cpp
 * @brief AFEPack Poisson 方程算例：混合边界条件。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 混合边界条件：同时处理 Dirichlet 与 Neumann 边界，展示边界积分和约束装配。
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
  constexpr int    right_boundary_mark = 2;

  struct MixedBoundaryResult
  {
    std::string  mesh_file;
    std::string  output_file;
    double       h = 0.0;
    int          cells = 0;
    unsigned int dofs = 0;
    double       l2_error = 0.0;
    double       h1_error = 0.0;
    double       residual = 0.0;
  };

  double exact_solution(const double *p)
  {
    return std::exp(p[0]) * std::sin(pi * p[1]);
  }

  std::vector<double> exact_gradient(const double *p)
  {
    return {std::exp(p[0]) * std::sin(pi * p[1]),
            pi * std::exp(p[0]) * std::cos(pi * p[1])};
  }

  double right_hand_side(const double *p)
  {
    return (pi * pi - 1.0) * exact_solution(p);
  }

  double neumann_data(const double *p)
  {
    return std::exp(1.0) * std::sin(pi * p[1]);
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

  std::vector<TemplateElement<double, dim, dim>> make_triangle_template(
    TemplateGeometry<dim> &geometry,
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

  double distance(const Point<dim> &a, const Point<dim> &b)
  {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    return std::sqrt(dx * dx + dy * dy);
  }

  void add_neumann_boundary_rhs(const EasyMesh &mesh,
                                const FEMSpace<double, dim> &space,
                                Vector<double> &rhs)
  {
    const double gauss_points[2] = {
      0.5 * (1.0 - 1.0 / std::sqrt(3.0)),
      0.5 * (1.0 + 1.0 / std::sqrt(3.0))};
    const double gauss_weights[2] = {0.5, 0.5};

    for (unsigned int side = 0; side < mesh.n_geometry(1); ++side)
      {
        const GeometryBM &edge = mesh.geometry(1, side);
        if (edge.boundaryMark() != right_boundary_mark)
          continue;

        const int vertex0 = edge.vertex(0);
        const int vertex1 = edge.vertex(1);
        const Point<dim> &p0 = mesh.point(vertex0);
        const Point<dim> &p1 = mesh.point(vertex1);
        const double length = distance(p0, p1);

        const auto &dof0 = space.geometryDof(0, vertex0);
        const auto &dof1 = space.geometryDof(0, vertex1);
        if (dof0.size() != 1 || dof1.size() != 1)
          throw std::runtime_error("expected one P1 dof on each boundary vertex");

        for (int q = 0; q < 2; ++q)
          {
            const double t = gauss_points[q];
            const double p[dim] = {(1.0 - t) * p0[0] + t * p1[0],
                                   (1.0 - t) * p0[1] + t * p1[1]};
            const double value = neumann_data(p) * gauss_weights[q] * length;
            rhs(dof0[0]) += value * (1.0 - t);
            rhs(dof1[0]) += value * t;
          }
      }
  }

  void apply_dirichlet_boundary(StiffMatrix<dim, double> &matrix,
                                FEMFunction<double, dim> &solution,
                                Vector<double> &rhs,
                                const FEMSpace<double, dim> &space)
  {
    BoundaryFunction<double, dim> bottom(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryFunction<double, dim> top(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
    BoundaryFunction<double, dim> left(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(bottom);
    boundary_admin.add(top);
    boundary_admin.add(left);
    boundary_admin.apply(matrix, solution, rhs);
  }

  double relative_residual(const StiffMatrix<dim, double> &matrix,
                           const Vector<double> &solution,
                           const Vector<double> &rhs)
  {
    Vector<double> ax(rhs.size());
    matrix.vmult(ax, solution);
    Vector<double> residual = rhs;
    residual.add(-1.0, ax);
    return residual.l2_norm() / std::max(rhs.l2_norm(), 1.0e-30);
  }

  MixedBoundaryResult solve_mixed_boundary_poisson(
    const std::string &mesh_file,
    const std::string &output_file,
    const bool write_solution,
    const double h)
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

    StiffMatrix<dim, double> stiff_matrix(fem_space);
    stiff_matrix.algebricAccuracy() = 4;
    stiff_matrix.build();

    FEMFunction<double, dim> solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    add_neumann_boundary_rhs(mesh, fem_space, rhs);
    apply_dirichlet_boundary(stiff_matrix, solution, rhs, fem_space);

    AMGSolver solver(stiff_matrix);
    solver.solve(solution, rhs, 1.0e-10, 400);

    if (write_solution)
      solution.writeOpenDXData(output_file);

    FunctionFunction<double> exact(&exact_solution, &exact_gradient);
    return {mesh_file,
            output_file,
            h,
            static_cast<int>(mesh.n_geometry(dim)),
            fem_space.n_dof(),
            Functional::L2Error(solution, exact, 4),
            Functional::H1SemiError(solution, exact, 4),
            relative_residual(stiff_matrix, solution, rhs)};
  }

  std::vector<std::pair<double, std::string>> default_meshes()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.200, mesh_dir + "/unit_square_h0p20"},
            {0.100, mesh_dir + "/unit_square_h0p10"},
            {0.050, mesh_dir + "/unit_square_h0p05"},
            {0.025, mesh_dir + "/unit_square_h0p025"}};
  }

  double convergence_rate(const double coarse_error,
                          const double fine_error,
                          const double coarse_h,
                          const double fine_h)
  {
    return std::log(coarse_error / fine_error) / std::log(coarse_h / fine_h);
  }

  void run_default_convergence_study()
  {
    const auto meshes = default_meshes();
    std::vector<MixedBoundaryResult> results;
    results.reserve(meshes.size());

    for (std::size_t i = 0; i < meshes.size(); ++i)
      {
        const bool write_solution = i + 1 == meshes.size();
        const std::string output_file =
          write_solution ?
            std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
              "/poisson_mixed_boundary_afepack_h0p025.dx" :
            "";
        results.push_back(solve_mixed_boundary_poisson(meshes[i].second,
                                                       output_file,
                                                       write_solution,
                                                       meshes[i].first));
      }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson P1 mixed Dirichlet/Neumann boundary\n";
    std::cout << "Exact solution exp(x) sin(pi y); Neumann boundary x=1.\n";
    std::cout << std::setw(12) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(16) << "L2"
              << std::setw(14) << "L2 rate" << std::setw(16) << "H1 semi"
              << std::setw(14) << "H1 rate" << std::setw(16)
              << "residual" << '\n';

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        std::cout << std::setw(12) << results[i].h << std::setw(10)
                  << results[i].cells << std::setw(10) << results[i].dofs
                  << std::setw(16) << results[i].l2_error;
        if (i == 0)
          std::cout << std::setw(14) << "-";
        else
          std::cout << std::setw(14)
                    << convergence_rate(results[i - 1].l2_error,
                                        results[i].l2_error,
                                        results[i - 1].h,
                                        results[i].h);

        std::cout << std::setw(16) << results[i].h1_error;
        if (i == 0)
          std::cout << std::setw(14) << "-";
        else
          std::cout << std::setw(14)
                    << convergence_rate(results[i - 1].h1_error,
                                        results[i].h1_error,
                                        results[i - 1].h,
                                        results[i].h);
        std::cout << std::setw(16) << results[i].residual << '\n';
      }

    std::cout << "Wrote " << results.back().output_file << '\n';
  }

  void print_single_result(const MixedBoundaryResult &result)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson P1 mixed Dirichlet/Neumann boundary\n";
    std::cout << "mesh: " << result.mesh_file << '\n';
    std::cout << "cells: " << result.cells << '\n';
    std::cout << "dofs: " << result.dofs << '\n';
    std::cout << "L2 error: " << result.l2_error << '\n';
    std::cout << "H1 seminorm error: " << result.h1_error << '\n';
    std::cout << "relative residual: " << result.residual << '\n';
    if (!result.output_file.empty())
      std::cout << "Wrote " << result.output_file << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output.dx]]\n";
      std::cerr << "       no arguments runs the default convergence study\n";
      return 2;
    }

  try
    {
      if (argc == 1)
        {
          run_default_convergence_study();
          return 0;
        }

      const std::string mesh_file = argv[1];
      const std::string output_file =
        argc > 2 ? argv[2] : std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                              "/poisson_mixed_boundary_afepack.dx";
      print_single_result(
        solve_mixed_boundary_poisson(mesh_file, output_file, true, 0.0));
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
