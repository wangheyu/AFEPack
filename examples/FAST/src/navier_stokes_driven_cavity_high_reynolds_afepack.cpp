/**
 * @file navier_stokes_driven_cavity_high_reynolds_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：驱动方腔流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 驱动方腔流：经典顶盖驱动方腔基准问题，用于验证速度-压力耦合和涡结构。
 * - 高雷诺数实验：在更强非线性和薄边界层条件下检验求解器稳健性。
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

#define main navier_stokes_driven_cavity_schur_hidden_main
#include "navier_stokes_driven_cavity_schur_afepack.cpp"
#undef main

namespace
{
  struct HighReynoldsCavityCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
    double      h = 0.0;
  };

  struct HighReynoldsCavityRow
  {
    HighReynoldsCavityCase mesh_case;
    CavityResult           result;
  };

  HighReynoldsCavityCase default_mesh_case()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {"unit_square_h0p10",
            mesh_dir + "/unit_square_h0p10",
            output_dir + "/navier_stokes_driven_cavity_high_reynolds_h0p10",
            0.10};
  }

  void print_result_row(const HighReynoldsCavityRow &row)
  {
    const auto &result = row.result;
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.system_nonzeros << ' '
              << result.picard_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.velocity_l2_norm << ' '
              << result.velocity_h1_seminorm << ' '
              << result.max_speed << ' '
              << result.pressure_min << ' '
              << result.pressure_max << ' '
              << result.pressure_mean << ' '
              << result.divergence_l2_norm << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 2)
    {
      std::cerr << "usage: " << argv[0] << " [viscosities]\n";
      return 2;
    }

  const auto viscosities =
    parse_viscosities(argc > 1 ? argv[1] : "0.02,0.01,0.005,0.002,0.001");

  try
    {
      const auto mesh_case = default_mesh_case();
      const auto results = solve_cavity_continuation(mesh_case.mesh_file,
                                                     mesh_case.output_prefix,
                                                     viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes regularized driven cavity "
                   "high-Reynolds Schur/AMG continuation "
                   "(book Example 7.1.3 / Figure 7.2 diagnostic)\n";
      std::cout << "preconditioner: tri(Aamg2,Mp-pcg)\n";
      std::cout << "mesh h nu Re cells velocity_dofs pressure_dofs "
                   "nonzeros Picard GMRES final_update final_gmres "
                   "velocity_L2 velocity_H1 max_speed pressure_min "
                   "pressure_max pressure_mean divergence_L2\n";
      for (const auto &result : results)
        print_result_row({mesh_case, result});
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
