/**
 * @file poisson_amg_sweep_afepack.cpp
 * @brief AFEPack Poisson 方程算例：AMG 参数扫描。
 *
 * @details
 * 本文件是 FAST 目录中的 Poisson 方程 迁移算例，关注二阶椭圆模型问题的有限元离散、误差评估和线性求解器对比。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - AMG 参数扫描：批量改变代数多重网格设置并比较残差收敛行为。
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

namespace
{
  constexpr int    dim = 2;
  constexpr double pi  = 3.141592653589793238462643383279502884;

  struct AMGConfig
  {
    unsigned int smooth_steps = 0;
    unsigned int cycles = 0;
  };

  struct AMGResult
  {
    double       h = 0.0;
    unsigned int cells = 0;
    unsigned int dofs = 0;
    unsigned int nonzeros = 0;
    unsigned int smooth_steps = 0;
    unsigned int cycles = 0;
    double       relative_residual = 0.0;
    double       residual_factor = 0.0;
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

    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(boundary1);
    boundary_admin.add(boundary2);
    boundary_admin.add(boundary3);
    boundary_admin.add(boundary4);
    boundary_admin.apply(stiff_matrix, solution, rhs);
  }

  double relative_residual(const SparseMatrix<double> &matrix,
                           const Vector<double> &x,
                           const Vector<double> &rhs)
  {
    const double rhs_norm = std::max(rhs.l2_norm(), 1.0e-30);

    Vector<double> ax(rhs.size());
    matrix.vmult(ax, x);

    Vector<double> residual = rhs;
    residual.add(-1.0, ax);
    return residual.l2_norm() / rhs_norm;
  }

  std::vector<std::pair<double, std::string>> default_meshes()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.200, mesh_dir + "/unit_square_h0p20"},
            {0.100, mesh_dir + "/unit_square_h0p10"},
            {0.050, mesh_dir + "/unit_square_h0p05"},
            {0.025, mesh_dir + "/unit_square_h0p025"}};
  }

  std::vector<AMGResult> run_mesh_file(const double h,
                                       const std::string &mesh_file,
                                       const std::vector<AMGConfig> &configs)
  {
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

    FEMFunction<double, dim> boundary_solution(fem_space);
    Vector<double> rhs;
    Operator::L2Discretize(&right_hand_side, fem_space, rhs, 4);
    apply_dirichlet_boundary(stiff_matrix, boundary_solution, rhs, fem_space);

    std::vector<AMGResult> results;
    results.reserve(configs.size());

    for (const auto &config : configs)
      {
        if (config.smooth_steps == 0 || config.cycles == 0)
          throw std::runtime_error("AMG smooth steps and cycles must be positive");

        FEMFunction<double, dim> solution(fem_space);
        const double initial_residual =
          relative_residual(stiff_matrix, solution, rhs);

        AMGSolver solver(stiff_matrix, 1.0e-12, config.smooth_steps);
        solver.solve(solution, rhs, 0.0, config.cycles - 1, 1);

        const double final_residual =
          relative_residual(stiff_matrix, solution, rhs);
        if (!std::isfinite(final_residual) ||
            final_residual >= initial_residual)
          throw std::runtime_error("AMG V-cycle did not reduce the residual");

        const double residual_factor =
          std::pow(final_residual / initial_residual,
                   1.0 / static_cast<double>(config.cycles));
        const double l2_error =
          Functional::L2Error(solution, FunctionFunction<double>(&exact_solution), 4);

        results.push_back({h,
                           mesh.n_geometry(dim),
                           fem_space.n_dof(),
                           static_cast<unsigned int>(
                             stiff_matrix.n_nonzero_elements()),
                           config.smooth_steps,
                           config.cycles,
                           final_residual,
                           residual_factor,
                           l2_error});
      }

    return results;
  }

  void print_results(const std::string &title,
                     const std::vector<AMGResult> &results)
  {
    std::cout << title << '\n';
    std::cout << std::setw(10) << "h" << std::setw(10) << "cells"
              << std::setw(10) << "dofs" << std::setw(10) << "nnz"
              << std::setw(10) << "smooth" << std::setw(10) << "cycles"
              << std::setw(14) << "rel res" << std::setw(14)
              << "factor" << std::setw(14) << "L2 error" << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(10) << result.h << std::setw(10)
                  << result.cells << std::setw(10) << result.dofs
                  << std::setw(10) << result.nonzeros << std::setw(10)
                  << result.smooth_steps << std::setw(10) << result.cycles
                  << std::setw(14) << result.relative_residual
                  << std::setw(14) << result.residual_factor
                  << std::setw(14) << result.l2_error << '\n';
      }
  }

  void run_default_study()
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack Poisson AMG V-cycle sweep\n";
    std::cout << "P1 triangles, manufactured solution, homogeneous Dirichlet boundary.\n";

    std::vector<AMGResult> mesh_results;
    for (const auto &mesh : default_meshes())
      {
        const auto partial =
          run_mesh_file(mesh.first, mesh.second, {{2, 4}});
        mesh_results.insert(mesh_results.end(), partial.begin(), partial.end());
      }
    print_results("Mesh sweep with 2 Gauss-Seidel smoothing steps and 4 cycles:",
                  mesh_results);

    const std::string finest_mesh =
      std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p025";
    const std::vector<AMGConfig> smoother_configs = {
      {1, 1}, {1, 2}, {1, 4}, {1, 8},
      {2, 1}, {2, 2}, {2, 4}, {2, 8},
      {3, 1}, {3, 2}, {3, 4}, {3, 8}};
    const auto smoother_results =
      run_mesh_file(0.025, finest_mesh, smoother_configs);
    print_results("Finest-grid smoother/cycle sweep:", smoother_results);

    const double strongest_residual =
      smoother_results.back().relative_residual;
    if (strongest_residual > 1.0e-8)
      throw std::runtime_error("strongest AMG sweep did not reach the expected residual");
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
