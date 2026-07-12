/**
 * @file stokes_colliding_schur_fine_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：碰撞流算例。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 碰撞流算例：设置相向流动边界条件，观察压力场和速度场在中心区域的耦合。
 * - Schur 补求解框架：显式组织速度块、压力块和 Schur 补近似以求解不可压流鞍点系统。
 * - 细网格算例：使用更高分辨率网格考察结果收敛性和计算代价。
 *
 * 网格与数据：主要依赖 meshes/afepack 中的 EasyMesh 输入；执行 make -C examples/FAST meshes 可重新生成网格文件。
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

#define main stokes_colliding_schur_hidden_main
#include "stokes_colliding_schur_afepack.cpp"
#undef main

namespace
{
  struct FineCollidingCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
    double      h = 0.0;
  };

  struct FineCollidingRow
  {
    FineCollidingCase mesh_case;
    StokesResult      result;
  };

  FineCollidingCase default_mesh_case()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {"unit_square_h0p05",
            mesh_dir + "/unit_square_h0p05",
            output_dir + "/stokes_colliding_schur_fine_h0p05",
            0.05};
  }

  std::vector<PreconditionerKind> default_preconditioners()
  {
    return {PreconditionerKind::block_diagonal,
            PreconditionerKind::schur_mass_pcg,
            PreconditionerKind::schur_mass_amg1,
            PreconditionerKind::schur_mass_amg2};
  }

  void print_result_row(const FineCollidingRow &row)
  {
    const auto &result = row.result;
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << result.preconditioner << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.system_nonzeros << ' '
              << result.gmres_iterations << ' '
              << result.gmres_relative_residual << ' '
              << result.velocity_l2_error << ' '
              << result.velocity_h1_error << ' '
              << result.pressure_l2_error << ' '
              << result.divergence_l2_norm << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc != 1)
    {
      std::cerr << "usage: " << argv[0] << '\n';
      return 2;
    }

  try
    {
      constexpr double gmres_tolerance = 1.0e-8;
      constexpr int    max_gmres_iterations = 5000;

      const auto mesh_case = default_mesh_case();
      std::vector<FineCollidingRow> rows;
      for (const auto preconditioner : default_preconditioners())
        rows.push_back({mesh_case,
                        solve_stokes(mesh_case.mesh_file,
                                     mesh_case.output_prefix,
                                     preconditioner,
                                     gmres_tolerance,
                                     max_gmres_iterations)});

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes colliding-flow fine-grid Schur/AMG "
                   "check (book Example 5.1.4)\n";
      std::cout << "mesh h preconditioner cells velocity_dofs pressure_dofs "
                   "nonzeros GMRES relres velocity_L2_error "
                   "velocity_H1_error pressure_L2_mod_constant "
                   "divergence_L2\n";
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
