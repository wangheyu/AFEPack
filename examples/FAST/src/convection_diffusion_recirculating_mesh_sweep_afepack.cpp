/**
 * @file convection_diffusion_recirculating_mesh_sweep_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：网格尺度扫描。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 网格尺度扫描：在多组网格输入上重复求解，比较自由度增长对误差和迭代次数的影响。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#define main convection_diffusion_recirculating_hidden_main
#include "convection_diffusion_recirculating_layer_afepack.cpp"
#undef main

namespace
{
  struct RecirculatingMeshCase
  {
    std::string name;
    std::string mesh_file;
    double      h = 0.0;
  };

  struct RecirculatingMeshSweepRow
  {
    RecirculatingMeshCase mesh_case;
    CaseResult            result;
  };

  std::vector<RecirculatingMeshCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{"unit_square_h0p20",
             mesh_dir + "/unit_square_h0p20",
             0.20},
            {"unit_square_h0p10",
             mesh_dir + "/unit_square_h0p10",
             0.10}};
  }

  void print_result_row(const RecirculatingMeshSweepRow &row)
  {
    const auto &result = row.result;
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << result.epsilon << ' '
              << (result.supg ? "SUPG" : "Galerkin") << ' '
              << result.cells << ' '
              << result.dofs << ' '
              << result.iterations << ' '
              << result.relative_residual << ' '
              << result.min_value << ' '
              << result.max_value << ' '
              << result.undershoot << ' '
              << result.overshoot << ' '
              << result.l2_norm << ' '
              << result.h1_seminorm << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [epsilon]\n";
      return 2;
    }

  const double epsilon = argc > 1 ? std::stod(argv[1]) : 1.0 / 200.0;
  if (!(epsilon > 0.0))
    {
      std::cerr << "epsilon must be positive\n";
      return 2;
    }

  try
    {
      std::vector<RecirculatingMeshSweepRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        {
          rows.push_back({mesh_case,
                          run_case(mesh_case.mesh_file,
                                   "",
                                   epsilon,
                                   false,
                                   false)});
          rows.push_back({mesh_case,
                          run_case(mesh_case.mesh_file,
                                   "",
                                   epsilon,
                                   true,
                                   false)});
        }

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack convection-diffusion recirculating "
                   "double-glazing mesh sweep (book Example 3.1.4)\n";
      std::cout << "mesh h epsilon method cells dofs GMRES relres min max "
                   "undershoot overshoot L2 H1\n";
      for (const auto &row : rows)
        print_result_row(row);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
