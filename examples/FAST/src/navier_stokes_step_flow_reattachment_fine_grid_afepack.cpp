/**
 * @file navier_stokes_step_flow_reattachment_fine_grid_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：细网格再附长度实验。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 细网格再附长度实验：在后台阶流中用细网格评估再附点位置和剪切量。
 * - 后台阶流：在带突扩几何的通道中模拟回流区和再附长度等典型流动量。
 * - 细网格算例：使用更高分辨率网格考察结果收敛性和计算代价。
 *
 * 网格与数据：主要依赖 后台阶通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <cmath>
#include <exception>
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
  struct FineGridStepFlowCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
    double      h = 0.0;
  };

  struct FineGridReattachmentRow
  {
    FineGridStepFlowCase mesh_case;
    StepFlowResult       result;
  };

  std::vector<FineGridStepFlowCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;
    return {{"step_channel_h0p50",
             mesh_dir + "/step_channel_h0p50",
             output_dir + "/navier_stokes_step_flow_reattach_fine_h0p50",
             0.50},
            {"step_channel_h0p25",
             mesh_dir + "/step_channel_h0p25",
             output_dir + "/navier_stokes_step_flow_reattach_fine_h0p25",
             0.25},
            {"step_channel_h0p125",
             mesh_dir + "/step_channel_h0p125",
             output_dir + "/navier_stokes_step_flow_reattach_fine_h0p125",
             0.125}};
  }

  void print_result_row(const FineGridReattachmentRow &row)
  {
    const auto &result = row.result;
    std::cout << row.mesh_case.name << ' '
              << row.mesh_case.h << ' '
              << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.picard_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << ' '
              << result.wall_shear_samples << ' '
              << result.wall_shear_min << ' '
              << result.wall_shear_max << ' '
              << result.reattachment_found << ' '
              << result.reattachment_x << '\n';
  }

  const FineGridReattachmentRow &
  row_for(const std::vector<FineGridReattachmentRow> &rows,
          const std::string                         &mesh_name,
          const double                               viscosity)
  {
    for (const auto &row : rows)
      if (row.mesh_case.name == mesh_name &&
          std::abs(row.result.viscosity - viscosity) < 1.0e-12)
        return row;
    throw std::runtime_error("missing fine-grid reattachment row");
  }

  void print_convergence_rows(
    const std::vector<FineGridReattachmentRow> &rows,
    const std::vector<double>                  &viscosities)
  {
    std::cout << "convergence nu Re x_h0p50 x_h0p25 x_h0p125 "
                 "delta_0p50_to_0p25 delta_0p25_to_0p125\n";
    for (const double viscosity : viscosities)
      {
        const auto &coarse = row_for(rows, "step_channel_h0p50", viscosity);
        const auto &middle = row_for(rows, "step_channel_h0p25", viscosity);
        const auto &fine = row_for(rows, "step_channel_h0p125", viscosity);

        std::cout << viscosity << ' '
                  << middle.result.reynolds_number << ' '
                  << coarse.result.reattachment_x << ' '
                  << middle.result.reattachment_x << ' '
                  << fine.result.reattachment_x << ' '
                  << middle.result.reattachment_x -
                       coarse.result.reattachment_x << ' '
                  << fine.result.reattachment_x -
                       middle.result.reattachment_x << '\n';
      }
  }

  void print_reattachment_richardson_rows(
    const std::vector<FineGridReattachmentRow> &rows,
    const std::vector<double>                  &viscosities)
  {
    std::cout << "reattachment_richardson nu Re delta_coarse_to_mid "
                 "delta_mid_to_fine abs_delta_ratio observed_order_abs "
                 "sign_change richardson_p1 fine_to_p1_distance\n";

    for (const double viscosity : viscosities)
      {
        const auto &coarse = row_for(rows, "step_channel_h0p50", viscosity);
        const auto &middle = row_for(rows, "step_channel_h0p25", viscosity);
        const auto &fine = row_for(rows, "step_channel_h0p125", viscosity);

        const double delta_coarse_to_mid =
          middle.result.reattachment_x - coarse.result.reattachment_x;
        const double delta_mid_to_fine =
          fine.result.reattachment_x - middle.result.reattachment_x;
        const double abs_coarse_delta = std::abs(delta_coarse_to_mid);
        const double abs_fine_delta = std::abs(delta_mid_to_fine);
        const bool   has_two_deltas =
          abs_coarse_delta > 1.0e-14 && abs_fine_delta > 1.0e-14;
        const int sign_change =
          (delta_coarse_to_mid * delta_mid_to_fine < 0.0) ? 1 : 0;
        const double richardson_p1 =
          fine.result.reattachment_x + delta_mid_to_fine;

        std::cout << viscosity << ' '
                  << middle.result.reynolds_number << ' '
                  << delta_coarse_to_mid << ' '
                  << delta_mid_to_fine << ' ';
        if (has_two_deltas)
          std::cout << abs_fine_delta / abs_coarse_delta << ' '
                    << std::log(abs_coarse_delta / abs_fine_delta) /
                         std::log(2.0);
        else
          std::cout << "NA NA";
        std::cout << ' ' << sign_change << ' '
                  << richardson_p1 << ' '
                  << std::abs(richardson_p1 - fine.result.reattachment_x)
                  << '\n';
      }
  }

  double positive_ratio_or_negative(const double numerator,
                                    const double denominator)
  {
    if (std::abs(denominator) <= 1.0e-14)
      return -1.0;
    return numerator / denominator;
  }

  void print_conservation_rows(
    const std::vector<FineGridReattachmentRow> &rows,
    const std::vector<double>                  &viscosities)
  {
    std::cout << "conservation_convergence nu Re abs_flux_h0p50 "
                 "abs_flux_h0p25 abs_flux_h0p125 "
                 "flux_ratio_0p50_to_0p25 flux_ratio_0p25_to_0p125 "
                 "div_h0p50 div_h0p25 div_h0p125 "
                 "div_ratio_0p50_to_0p25 div_ratio_0p25_to_0p125\n";

    for (const double viscosity : viscosities)
      {
        const auto &coarse = row_for(rows, "step_channel_h0p50", viscosity);
        const auto &middle = row_for(rows, "step_channel_h0p25", viscosity);
        const auto &fine = row_for(rows, "step_channel_h0p125", viscosity);

        const double flux_coarse = std::abs(coarse.result.flux_imbalance);
        const double flux_middle = std::abs(middle.result.flux_imbalance);
        const double flux_fine = std::abs(fine.result.flux_imbalance);
        const double div_coarse = coarse.result.divergence_l2_norm;
        const double div_middle = middle.result.divergence_l2_norm;
        const double div_fine = fine.result.divergence_l2_norm;

        std::cout << viscosity << ' '
                  << middle.result.reynolds_number << ' '
                  << flux_coarse << ' '
                  << flux_middle << ' '
                  << flux_fine << ' '
                  << positive_ratio_or_negative(flux_middle, flux_coarse)
                  << ' '
                  << positive_ratio_or_negative(flux_fine, flux_middle)
                  << ' '
                  << div_coarse << ' '
                  << div_middle << ' '
                  << div_fine << ' '
                  << positive_ratio_or_negative(div_middle, div_coarse)
                  << ' '
                  << positive_ratio_or_negative(div_fine, div_middle)
                  << '\n';
      }
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
      std::vector<FineGridReattachmentRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        {
          const auto results = solve_step_flow_continuation(mesh_case.mesh_file,
                                                            mesh_case.output_prefix,
                                                            viscosities);
          for (const auto &result : results)
            rows.push_back({mesh_case, result});
        }

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes expansion step-flow "
                   "fine-grid reattachment check "
                   "(book Exercise 7.5)\n";
      std::cout << "preconditioner: tri(Aamg2,Mp-pcg)\n";
      std::cout << "mesh h nu Re cells velocity_dofs pressure_dofs "
                   "Picard GMRES final_update final_gmres flux_imbalance "
                   "divergence_L2 shear_samples shear_min shear_max "
                   "reattachment_found reattachment_x\n";
      for (const auto &row : rows)
        print_result_row(row);
      print_convergence_rows(rows, viscosities);
      print_reattachment_richardson_rows(rows, viscosities);
      print_conservation_rows(rows, viscosities);
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
