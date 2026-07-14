/**
 * @file navier_stokes_blasius_mesh_compare_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：Blasius 网格对比。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - Blasius 网格对比：比较不同 Blasius 网格设置下的阻力、通量和残差表现。
 * - Blasius 边界层代理：使用不可压 Navier-Stokes 方程构造平板边界层的数值代理模型。
 *
 * 网格与数据：主要依赖 平板边界层系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
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

#define main navier_stokes_blasius_hidden_main
#include "navier_stokes_blasius_afepack.cpp"
#undef main

namespace
{
  struct BlasiusMeshCase
  {
    std::string name;
    std::string mesh_file;
    std::string output_prefix;
  };

  struct BlasiusMeshCompareRow
  {
    BlasiusMeshCase mesh_case;
    BlasiusResult   result;
  };

  std::vector<BlasiusMeshCase> default_mesh_cases()
  {
    const std::string mesh_dir = AFEPACK_DEFAULT_MESH_DIR;
    const std::string output_dir = AFEPACK_DEFAULT_OUTPUT_DIR;

    const std::string generated_mesh_file =
      mesh_dir + "/generated_blasius_plate_compare_stretched";
    const GeneratedMesh stretched_mesh = generate_stretched_blasius_mesh();
    write_generated_mesh(stretched_mesh, generated_mesh_file);

    return {{"easymesh_h0p25",
             mesh_dir + "/blasius_plate_h0p25",
             output_dir + "/navier_stokes_blasius_mesh_compare_easymesh"},
            {"generated_stretched",
             generated_mesh_file,
             output_dir + "/navier_stokes_blasius_mesh_compare_stretched"}};
  }

  double average_x1_thickness(const BlasiusResult &result)
  {
    return average_valid_thickness(result.upper_layer_x1,
                                   result.lower_layer_x1);
  }

  double average_x4_thickness(const BlasiusResult &result)
  {
    return average_valid_thickness(result.upper_layer_x4,
                                   result.lower_layer_x4);
  }

  void print_result_row(const BlasiusMeshCompareRow &row)
  {
    const auto &result = row.result;
    const double delta_x1 = average_x1_thickness(result);
    const double delta_x4 = average_x4_thickness(result);

    std::cout << row.mesh_case.name << ' '
              << result.viscosity << ' '
              << result.reynolds_number << ' '
              << result.cells << ' '
              << result.velocity_dofs << ' '
              << result.pressure_dofs << ' '
              << result.system_nonzeros << ' '
              << result.nonlinear_iterations << ' '
              << result.total_gmres_iterations << ' '
              << result.final_relative_update << ' '
              << result.final_gmres_relative_residual << ' '
              << result.final_nonlinear_relative_residual << ' '
              << result.velocity_l2_norm << ' '
              << result.velocity_h1_seminorm << ' '
              << result.max_speed << ' '
              << result.flux_imbalance << ' '
              << result.divergence_l2_norm << ' '
              << delta_x1 << ' '
              << delta_x4 << ' '
              << normalized_layer_thickness(delta_x1,
                                            result.viscosity,
                                            1.0)
              << ' '
              << normalized_layer_thickness(delta_x4,
                                            result.viscosity,
                                            4.0)
              << '\n';
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
      std::vector<BlasiusMeshCompareRow> rows;
      for (const auto &mesh_case : default_mesh_cases())
        {
          const auto results = solve_blasius_continuation(mesh_case.mesh_file,
                                                         mesh_case.output_prefix,
                                                         viscosities);
          for (const auto &result : results)
            rows.push_back({mesh_case, result});
        }

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes Blasius thin-plate proxy "
                   "mesh comparison (book Example 7.1.4 / Exercise 7.2)\n";
      std::cout << "mesh nu Re cells velocity_dofs pressure_dofs nonzeros "
                   "Nonlinear GMRES final_update final_gmres final_nonlinear "
                   "velocity_L2 "
                   "velocity_H1 max_speed "
                   "flux_imbalance divergence_L2 delta_x1_avg delta_x4_avg "
                   "delta_x1_over_sqrt_nu delta_x4_over_sqrt_4nu\n";
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
