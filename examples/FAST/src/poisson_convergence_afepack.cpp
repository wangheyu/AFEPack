/**
 * @file poisson_convergence_afepack.cpp
 * @brief AFEPack Poisson 方程算例：收敛阶验证。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 收敛阶验证：在多级网格上计算误差和经验收敛阶，验证有限元空间的理论阶数。
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

  struct PoissonResult
  {
    std::string  mesh_file;
    std::string  output_file;
    double       h;
    int          cells;
    unsigned int dofs;
    double       l2_error;
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

  PoissonResult solve_poisson(const std::string &mesh_file,
                              const std::string &output_file,
                              const bool         write_solution,
                              const double       target_h)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    ensure_template_path();

    TemplateGeometry<dim> triangle_template_geometry;
    triangle_template_geometry.readData("triangle.tmp_geo");

    CoordTransform<dim, dim> triangle_coord_transform;
    triangle_coord_transform.readData("triangle.crd_trs");

    TemplateDOF<dim> triangle_template_dof(triangle_template_geometry);
    triangle_template_dof.readData("triangle.1.tmp_dof");

    BasisFunctionAdmin<double, dim, dim> triangle_basis_function(
      triangle_template_dof);
    triangle_basis_function.readData("triangle.1.bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(triangle_template_geometry,
                               triangle_template_dof,
                               triangle_coord_transform,
                               triangle_basis_function);

    FEMSpace<double, dim> fem_space;
    fem_space.reinit(mesh, template_element);

    const int n_element = mesh.n_geometry(dim);
    fem_space.element().resize(n_element);
    for (int i = 0; i < n_element; ++i)
      fem_space.element(i).reinit(fem_space, i, 0);

    fem_space.buildElement();
    fem_space.buildDof();
    fem_space.buildDofBoundaryMark();

    StiffMatrix<dim, double> stiff_matrix(fem_space);
    stiff_matrix.algebricAccuracy() = 4;
    stiff_matrix.build();

    FEMFunction<double, dim> solution(fem_space);
    Vector<double>           rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);

    BoundaryFunction<double, dim> boundary1(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
    BoundaryFunction<double, dim> boundary2(
      BoundaryConditionInfo::DIRICHLET, 2, &exact_solution);
    BoundaryFunction<double, dim> boundary3(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
    BoundaryFunction<double, dim> boundary4(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);

    BoundaryConditionAdmin<double, dim> boundary_admin(fem_space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.apply(stiff_matrix, solution, rhs);

    AMGSolver solver(stiff_matrix);
    solver.solve(solution, rhs, 1.0e-10, 400);

    if (write_solution)
      solution.writeOpenDXData(output_file);

    return {mesh_file,
            output_file,
            target_h,
            n_element,
            fem_space.n_dof(),
            Functional::L2Error(solution,
                                FunctionFunction<double>(&exact_solution),
                                4)};
  }

  std::vector<std::pair<double, std::string>> default_convergence_meshes()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.200, mesh_dir + "/unit_square_h0p20"},
            {0.100, mesh_dir + "/unit_square_h0p10"},
            {0.050, mesh_dir + "/unit_square_h0p05"},
            {0.025, mesh_dir + "/unit_square_h0p025"}};
  }

  double l2_rate(const PoissonResult &coarse, const PoissonResult &fine)
  {
    return std::log(coarse.l2_error / fine.l2_error) /
           std::log(coarse.h / fine.h);
  }

  void print_single_result(const PoissonResult &result)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson P1 manufactured solution\n";
    std::cout << "mesh: " << result.mesh_file << '\n';
    std::cout << "cells: " << result.cells << '\n';
    std::cout << "dofs: " << result.dofs << '\n';
    std::cout << "L2 error: " << result.l2_error << '\n';
    if (!result.output_file.empty())
      std::cout << "Wrote " << result.output_file << '\n';
  }

  void run_default_convergence_study()
  {
    const auto meshes = default_convergence_meshes();
    std::vector<PoissonResult> results;
    results.reserve(meshes.size());

    for (std::size_t i = 0; i < meshes.size(); ++i)
      {
        const bool write_solution = i + 1 == meshes.size();
        const std::string output_file =
          write_solution ?
            std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
              "/poisson_afepack_h0p025.dx" :
            "";
        results.push_back(solve_poisson(meshes[i].second,
                                         output_file,
                                         write_solution,
                                         meshes[i].first));
      }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson P1 convergence study\n";
    std::cout << std::setw(12) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(16) << "L2"
              << std::setw(14) << "rate" << '\n';

    for (std::size_t i = 0; i < results.size(); ++i)
      {
        std::cout << std::setw(12) << results[i].h << std::setw(10)
                  << results[i].cells << std::setw(10) << results[i].dofs
                  << std::setw(16) << results[i].l2_error;
        if (i == 0)
          std::cout << std::setw(14) << "-";
        else
          std::cout << std::setw(14) << l2_rate(results[i - 1], results[i]);
        std::cout << '\n';
      }

    std::cout << "Wrote " << results.back().output_file << '\n';
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

      const std::string mesh_file   = argv[1];
      const std::string output_file =
        argc > 2 ? argv[2] : std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                              "/poisson_afepack.dx";
      print_single_result(solve_poisson(mesh_file, output_file, true, 0.0));
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
