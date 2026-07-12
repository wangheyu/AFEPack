/**
 * @file convection_diffusion_recirculating_grid_diffusion_afepack.cpp
 * @brief AFEPack 对流扩散方程算例：回流问题网格扩散实验。
 *
 * @details
 * 本文件是 FAST 目录中的 对流扩散方程 迁移算例，关注标量输运问题中扩散占优、对流占优以及边界层/内层结构的离散表现。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 回流问题网格扩散实验：比较回流问题中网格分辨率与扩散参数变化对数值解的影响。
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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define main convection_diffusion_recirculating_hidden_main
#include "convection_diffusion_recirculating_layer_afepack.cpp"
#undef main

namespace
{
  struct RecirculatingGridDiffusionCase
  {
    std::string name;
    std::string mesh_file;
    double      h = 0.0;
  };

  struct RecirculatingGridDiffusionRow
  {
    RecirculatingGridDiffusionCase mesh_case;
    CaseResult                     result;
  };

  std::vector<RecirculatingGridDiffusionCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    return {{"unit_square_h0p20",
             mesh_dir + "/unit_square_h0p20",
             0.20},
            {"unit_square_h0p10",
             mesh_dir + "/unit_square_h0p10",
             0.10}};
  }

  std::vector<double> parse_epsilons(const std::string &text)
  {
    std::vector<double> epsilons;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ','))
      {
        if (token.empty())
          continue;
        const double epsilon = std::stod(token);
        if (!(epsilon > 0.0))
          throw std::runtime_error("epsilon values must be positive");
        epsilons.push_back(epsilon);
      }

    if (epsilons.empty())
      throw std::runtime_error("empty epsilon list");
    return epsilons;
  }

  void print_result_row(const RecirculatingGridDiffusionRow &row)
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
      std::cerr << "usage: " << argv[0] << " [epsilons]\n";
      return 2;
    }

  try
    {
      const auto epsilons =
        parse_epsilons(argc > 1 ? argv[1] : "0.01,0.005,0.002");

      std::vector<RecirculatingGridDiffusionRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        for (const double epsilon : epsilons)
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
                   "double-glazing grid/diffusion sweep "
                   "(book Example 3.1.4)\n";
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
