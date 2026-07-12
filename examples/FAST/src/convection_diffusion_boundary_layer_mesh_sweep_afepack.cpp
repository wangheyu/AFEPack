/**
 * @file convection_diffusion_boundary_layer_mesh_sweep_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：边界层算例。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 边界层算例：展示靠近边界的薄层结构对网格尺度、扩散系数和离散稳定性的要求。
 * - 网格尺度扫描：在多组网格输入上重复求解，比较自由度增长对误差和迭代次数的影响。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#define main convection_diffusion_boundary_layer_default_main
#include "convection_diffusion_boundary_layer_afepack.cpp"
#undef main

namespace
{
  struct MeshCase
  {
    double      h = 0.0;
    std::string mesh_file;
  };

  struct MeshSweepResult
  {
    double       h = 0.0;
    double       epsilon = 0.0;
    double       cell_peclet = 0.0;
    bool         supg = false;
    int          cells = 0;
    unsigned int dofs = 0;
    int          iterations = 0;
    double       relative_residual = 0.0;
    double       l2_error = 0.0;
    double       h1_error = 0.0;
  };

  std::vector<MeshCase> default_meshes()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{0.200, mesh_dir + "/unit_square_h0p20"},
            {0.100, mesh_dir + "/unit_square_h0p10"},
            {0.050, mesh_dir + "/unit_square_h0p05"},
            {0.025, mesh_dir + "/unit_square_h0p025"}};
  }

  MeshSweepResult run_mesh_case(const MeshCase &mesh,
                                const double epsilon,
                                const bool use_supg)
  {
    const CaseResult result =
      run_case(mesh.mesh_file, "", epsilon, use_supg, false);

    return {mesh.h,
            epsilon,
            beta_norm() * mesh.h / (2.0 * epsilon),
            result.supg,
            result.cells,
            result.dofs,
            result.iterations,
            result.relative_residual,
            result.l2_error,
            result.h1_error};
  }

  void print_mesh_sweep(const std::vector<MeshSweepResult> &results)
  {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "AFEPack convection-diffusion boundary-layer mesh sweep\n";
    std::cout << "epsilon=1e-3, exact exponential outflow layer, P1 triangles.\n";
    std::cout << std::setw(10) << "h" << ' ' << std::setw(12)
              << "method" << ' ' << std::setw(12) << "Pe_h" << ' '
              << std::setw(10) << "cells" << ' ' << std::setw(10)
              << "dofs" << ' ' << std::setw(8) << "iter" << ' '
              << std::setw(14) << "relres" << ' ' << std::setw(14)
              << "L2" << ' ' << std::setw(14) << "H1-semi" << '\n';

    for (const auto &result : results)
      {
        std::cout << std::setw(10) << result.h << ' ' << std::setw(12)
                  << (result.supg ? "SUPG" : "Galerkin") << ' '
                  << std::setw(12) << result.cell_peclet << ' '
                  << std::setw(10) << result.cells << ' ' << std::setw(10)
                  << result.dofs << ' ' << std::setw(8)
                  << result.iterations << ' ' << std::setw(14)
                  << result.relative_residual << ' ' << std::setw(14)
                  << result.l2_error << ' ' << std::setw(14)
                  << result.h1_error << '\n';
      }
  }
}

int main()
{
  try
    {
      const double epsilon = 1.0e-3;
      std::vector<MeshSweepResult> results;
      for (const auto &mesh : default_meshes())
        {
          results.push_back(run_mesh_case(mesh, epsilon, false));
          results.push_back(run_mesh_case(mesh, epsilon, true));
        }

      print_mesh_sweep(results);

      if (results.empty() || results.back().l2_error > 1.0e-1)
        throw std::runtime_error("SUPG finest-grid L2 error is too large");
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
