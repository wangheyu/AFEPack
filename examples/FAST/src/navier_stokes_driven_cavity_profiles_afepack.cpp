/**
 * @file navier_stokes_driven_cavity_profiles_afepack.cpp
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
 * - 剖面数据输出：沿指定截线输出速度或压力剖面，便于与基准数据和图形对比。
 *
 * 网格与数据：主要依赖 单位方形或规则通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#define main navier_stokes_driven_cavity_schur_hidden_main
#include "navier_stokes_driven_cavity_schur_afepack.cpp"
#undef main

namespace
{
  double value_at_coordinate(const std::vector<double> &coordinates,
                             const std::vector<double> &values,
                             const double target)
  {
    for (std::size_t i = 0; i < coordinates.size() && i < values.size(); ++i)
      if (std::abs(coordinates[i] - target) < 1.0e-12)
        return values[i];
    throw std::runtime_error("requested centerline coordinate was not sampled");
  }

  void print_summary_row(const CavityResult &result)
  {
    const auto ux_minmax =
      std::minmax_element(result.vertical_centerline_ux.begin(),
                          result.vertical_centerline_ux.end());
    const auto uy_minmax =
      std::minmax_element(result.horizontal_centerline_uy.begin(),
                          result.horizontal_centerline_uy.end());

    std::cout << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.picard_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.velocity_l2_norm << ' '
              << result.velocity_h1_seminorm << ' '
              << result.max_speed << ' '
              << result.divergence_l2_norm << ' '
              << *ux_minmax.first << ' '
              << *ux_minmax.second << ' '
              << *uy_minmax.first << ' '
              << *uy_minmax.second << ' '
              << value_at_coordinate(result.vertical_centerline_y,
                                     result.vertical_centerline_ux,
                                     0.0) << ' '
              << value_at_coordinate(result.horizontal_centerline_x,
                                     result.horizontal_centerline_uy,
                                     0.0) << '\n';
  }

  void print_profile_samples(const CavityResult &result)
  {
    for (std::size_t i = 0; i < result.vertical_centerline_y.size(); ++i)
      std::cout << result.viscosity << ' '
                << result.reynolds_number << ' '
                << "vertical_x0 "
                << result.vertical_centerline_y[i] << ' '
                << "u_x "
                << result.vertical_centerline_ux[i] << '\n';

    for (std::size_t i = 0; i < result.horizontal_centerline_x.size(); ++i)
      std::cout << result.viscosity << ' '
                << result.reynolds_number << ' '
                << "horizontal_y0 "
                << result.horizontal_centerline_x[i] << ' '
                << "u_y "
                << result.horizontal_centerline_uy[i] << '\n';
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
    parse_viscosities(argc > 1 ? argv[1] : "0.1,0.02");

  try
    {
      const std::string mesh_file =
        std::string(AFEPACK_DEFAULT_MESH_DIR) + "/unit_square_h0p10";
      const std::string output_prefix =
        std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
        "/navier_stokes_driven_cavity_profiles_afepack";

      const auto results =
        solve_cavity_continuation(mesh_file, output_prefix, viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes regularized driven cavity "
                   "centerline profile diagnostic (book Example 7.1.3)\n";
      std::cout << "mesh: unit_square_h0p10\n";
      std::cout << "preconditioner: tri(Aamg2,Mp-pcg)\n";
      std::cout << "summary nu Re cells velocity_dofs pressure_dofs "
                   "Picard GMRES final_update final_gmres velocity_L2 "
                   "velocity_H1 max_speed divergence_L2 ux_vertical_min "
                   "ux_vertical_max uy_horizontal_min uy_horizontal_max "
                   "center_ux center_uy\n";
      for (const auto &result : results)
        print_summary_row(result);

      std::cout << "profile_samples nu Re line coordinate component value\n";
      for (const auto &result : results)
        print_profile_samples(result);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
