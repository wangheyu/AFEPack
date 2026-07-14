/**
 * @file navier_stokes_step_flow_outflow_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：后台阶流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 后台阶流：在带突扩几何的通道中模拟回流区和再附长度等典型流动量。
 * - 出流边界对比：比较不同出流边界处理对后台阶流压力和回流长度的影响。
 *
 * 网格与数据：主要依赖 后台阶通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define main navier_stokes_step_flow_schur_hidden_main
#include "navier_stokes_step_flow_schur_afepack.cpp"
#undef main

namespace
{
  struct StepFlowMeshCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
    double      outlet_x = 0.0;
  };

  struct OutflowComparisonRow
  {
    StepFlowResult original;
    StepFlowResult long_domain;
    double         max_profile_diff = 0.0;
    double         rms_profile_diff = 0.0;
  };

  std::vector<StepFlowMeshCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {{"step_channel_h0p25",
             mesh_dir + "/step_channel_h0p25",
             output_dir + "/navier_stokes_step_flow_outflow_original",
             5.0},
            {"step_channel_long_h0p25",
             mesh_dir + "/step_channel_long_h0p25",
             output_dir + "/navier_stokes_step_flow_outflow_long",
             10.0}};
  }

  std::vector<StepFlowResult>
  solve_case(const StepFlowMeshCase &mesh_case,
             const std::vector<double> &viscosities)
  {
    return solve_step_flow_continuation(mesh_case.mesh_file,
                                        mesh_case.output_prefix,
                                        viscosities);
  }

  OutflowComparisonRow
  compare_profiles(const StepFlowResult &original,
                   const StepFlowResult &long_domain)
  {
    if (original.profile_ux.empty() || long_domain.profile_ux.empty())
      throw std::runtime_error("missing ux profile data");

    if (original.profile_ux.size() != long_domain.profile_ux.size() ||
        original.profile_y.size() != long_domain.profile_y.size() ||
        original.profile_ux.size() != original.profile_y.size())
      throw std::runtime_error("incompatible ux profile sizes");

    double max_diff = 0.0;
    double sum_squares = 0.0;
    for (std::size_t i = 0; i < original.profile_ux.size(); ++i)
      {
        if (std::abs(original.profile_y[i] - long_domain.profile_y[i]) >
            1.0e-12)
          throw std::runtime_error("incompatible ux profile sampling points");

        const double diff =
          std::abs(original.profile_ux[i] - long_domain.profile_ux[i]);
        max_diff = std::max(max_diff, diff);
        sum_squares += diff * diff;
      }

    OutflowComparisonRow row;
    row.original = original;
    row.long_domain = long_domain;
    row.max_profile_diff = max_diff;
    row.rms_profile_diff =
      std::sqrt(sum_squares / static_cast<double>(original.profile_ux.size()));
    return row;
  }

  void print_result_row(const OutflowComparisonRow &row)
  {
    std::cout << row.original.viscosity << ' '
              << row.original.reynolds_number << ' '
              << row.original.picard_iterations << ' '
              << row.long_domain.picard_iterations << ' '
              << row.original.total_gmres_iterations << ' '
              << row.long_domain.total_gmres_iterations << ' '
              << row.original.flux_imbalance << ' '
              << row.long_domain.flux_imbalance << ' '
              << row.original.reattachment_x << ' '
              << row.long_domain.reattachment_x << ' '
              << row.max_profile_diff << ' '
              << row.rms_profile_diff << '\n';
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
    parse_viscosities(argc > 1 ? argv[1] : "0.05,0.02");

  try
    {
      const auto mesh_cases = default_mesh_cases();
      const auto original_results = solve_case(mesh_cases[0], viscosities);
      const auto long_results = solve_case(mesh_cases[1], viscosities);

      if (original_results.size() != long_results.size())
        throw std::runtime_error("incompatible continuation result sizes");

      std::vector<OutflowComparisonRow> rows;
      rows.reserve(original_results.size());
      for (std::size_t i = 0; i < original_results.size(); ++i)
        rows.push_back(compare_profiles(original_results[i], long_results[i]));

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes expansion step-flow "
                   "outflow-length comparison (book Exercise 7.4)\n";
      std::cout << "preconditioner: tri(Aamg2,Mp-pcg)\n";
      std::cout << "profile sample: x=5, "
                   "y=-0.75,-0.5,-0.25,0,0.25,0.5,0.75\n";
      std::cout << "original outlet x: " << mesh_cases[0].outlet_x
                << " long outlet x: " << mesh_cases[1].outlet_x << '\n';
      std::cout << "nu Re Picard_x5 Picard_x10 GMRES_x5 GMRES_x10 "
                   "flux_x5 flux_x10 reattach_x5 reattach_x10 "
                   "max_profile_diff rms_profile_diff\n";
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
