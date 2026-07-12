/**
 * @file stokes_contracting_channel_profiles_afepack.cpp
 * @brief AFEPack 不可压 Stokes 方程算例：收缩通道流。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Stokes 方程 迁移算例，关注低雷诺数流动的速度-压力混合有限元离散与鞍点线性系统。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - 收缩通道流：在收缩几何中考察入口、出口和壁面边界对流动结构的影响。
 * - 剖面数据输出：沿指定截线输出速度或压力剖面，便于与基准数据和图形对比。
 *
 * 网格与数据：主要依赖 收缩通道网格、收缩通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>

#include <lapacke.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef AFEPACK_DEFAULT_MESH_DIR
#  define AFEPACK_DEFAULT_MESH_DIR "build/meshes/afepack"
#endif

#ifndef AFEPACK_DEFAULT_OUTPUT_DIR
#  define AFEPACK_DEFAULT_OUTPUT_DIR "build"
#endif

namespace
{
  constexpr int    dim = 2;
  constexpr double viscosity = 1.0;
  constexpr double boundary_tolerance = 1.0e-12;
  constexpr double stabilization_parameter = 1.0;
  constexpr double outlet_x = 4.0;

  enum class ContractingChannelGeometry
  {
    linear_proxy,
    quadratic_bow_tie_proxy,
    ifiss_bow_tie_half
  };

  ContractingChannelGeometry active_channel_geometry =
    ContractingChannelGeometry::linear_proxy;

  struct VerticalProfileSample
  {
    double y = 0.0;
    double ux = 0.0;
  };

  struct VerticalProfile
  {
    double x = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;
    double flow_rate = 0.0;
    double ux_min = 0.0;
    double ux_center = 0.0;
    double ux_max = 0.0;
    std::vector<VerticalProfileSample> samples;
  };

  struct ContractingChannelResult
  {
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    double       velocity_l2_norm = 0.0;
    double       velocity_h1_seminorm = 0.0;
    double       max_speed = 0.0;
    double       inlet_rate = 0.0;
    double       outlet_rate = 0.0;
    double       flux_imbalance = 0.0;
    double       pressure_min = 0.0;
    double       pressure_max = 0.0;
    double       pressure_mean = 0.0;
    double       divergence_l2_norm = 0.0;
    double       max_cell_divergence_integral = 0.0;
    double       cell_mean_divergence_l2_norm = 0.0;
    std::vector<VerticalProfile> profiles;
  };

  struct CellDivergenceDiagnostic
  {
    double max_integral = 0.0;
    double mean_l2_norm = 0.0;
  };

  bool near(const double value, const double target)
  {
    return std::abs(value - target) < boundary_tolerance;
  }

  double channel_half_width(const double x)
  {
    switch (active_channel_geometry)
      {
      case ContractingChannelGeometry::linear_proxy:
        return 1.0 - 0.125 * x;
      case ContractingChannelGeometry::quadratic_bow_tie_proxy:
        {
          const double t = 1.0 - x / outlet_x;
          return 0.5 + 0.5 * t * t;
        }
      case ContractingChannelGeometry::ifiss_bow_tie_half:
        return 0.5 * (1.0 - 3.0 * x / 16.0);
      }

    return 1.0 - 0.125 * x;
  }

  double channel_y_min(const double x)
  {
    return active_channel_geometry ==
             ContractingChannelGeometry::ifiss_bow_tie_half ?
             0.0 :
             -channel_half_width(x);
  }

  double channel_y_max(const double x)
  {
    return active_channel_geometry ==
             ContractingChannelGeometry::ifiss_bow_tie_half ?
             1.0 - 3.0 * x / 16.0 :
             channel_half_width(x);
  }

  double channel_area()
  {
    switch (active_channel_geometry)
      {
      case ContractingChannelGeometry::linear_proxy:
        return 6.0;
      case ContractingChannelGeometry::quadratic_bow_tie_proxy:
        return outlet_x + outlet_x / 3.0;
      case ContractingChannelGeometry::ifiss_bow_tie_half:
        return 2.5;
      }

    return 6.0;
  }

  double reference_flow_rate()
  {
    return active_channel_geometry ==
             ContractingChannelGeometry::ifiss_bow_tie_half ?
             2.0 / 3.0 :
             4.0 / 3.0;
  }

  double profile_sampling_density()
  {
    return active_channel_geometry ==
             ContractingChannelGeometry::ifiss_bow_tie_half ?
             64.0 :
             8.0;
  }

  void set_contracting_channel_geometry(
    const ContractingChannelGeometry geometry)
  {
    active_channel_geometry = geometry;
  }

  double inflow_velocity_x(const double *p)
  {
    if (active_channel_geometry ==
        ContractingChannelGeometry::ifiss_bow_tie_half)
      return 4.0 * p[1] * (1.0 - p[1]);
    return 1.0 - p[1] * p[1];
  }

  double velocity_boundary_value(const int component, const double *p)
  {
    if (component == 0 && near(p[0], 0.0))
      return inflow_velocity_x(p);
    return 0.0;
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
  make_triangle_template_by_name(
    const std::string &template_name,
    TemplateGeometry<dim> &geometry,
    CoordTransform<dim, dim> &transform,
    TemplateDOF<dim> &dof,
    BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle." + template_name + ".tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle." + template_name + ".bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(geometry, dof, transform, basis);
    return template_element;
  }

  std::vector<TemplateElement<double, dim, dim>>
  make_triangle_template(const int order,
                         TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    return make_triangle_template_by_name(std::to_string(order),
                                          geometry,
                                          transform,
                                          dof,
                                          basis);
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

  int velocity_index(const int component,
                     const int dof,
                     const int n_velocity_dofs)
  {
    return component * n_velocity_dofs + dof;
  }

  int pressure_index(const int dof, const int n_velocity_dofs)
  {
    return 2 * n_velocity_dofs + dof;
  }

  void add_entry(std::vector<double> &matrix,
                 const int n,
                 const int row,
                 const int col,
                 const double value)
  {
    matrix[static_cast<std::size_t>(row) * n + col] += value;
  }

  void apply_dirichlet_constraint(std::vector<double> &matrix,
                                  std::vector<double> &rhs,
                                  const int n,
                                  const int index,
                                  const double value)
  {
    for (int row = 0; row < n; ++row)
      if (row != index)
        rhs[row] -= matrix[static_cast<std::size_t>(row) * n + index] *
                    value;

    for (int j = 0; j < n; ++j)
      {
        matrix[static_cast<std::size_t>(index) * n + j] = 0.0;
        matrix[static_cast<std::size_t>(j) * n + index] = 0.0;
      }
    matrix[static_cast<std::size_t>(index) * n + index] = 1.0;
    rhs[index] = value;
  }

  bool is_velocity_dirichlet_dof(const FEMSpace<double, dim> &velocity_space,
                                 const int dof)
  {
    if (velocity_space.dofBoundaryMark(dof) == 0)
      return false;

    const auto &point = velocity_space.dofInfo(dof).interp_point;
    const double x = point[0];
    const double y = point[1];
    const double y_min = channel_y_min(outlet_x);
    const double y_max = channel_y_max(outlet_x);

    if (near(x, outlet_x) && !near(y, y_min) && !near(y, y_max))
      return false;

    return true;
  }

  void solve_dense_system(std::vector<double> &matrix,
                          std::vector<double> &rhs,
                          const int n)
  {
    std::vector<int> pivots(n);
    const int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,
                                   n,
                                   1,
                                   matrix.data(),
                                   n,
                                   pivots.data(),
                                   rhs.data(),
                                   1);
    if (info != 0)
      throw std::runtime_error("LAPACKE_dgesv failed with info=" +
                               std::to_string(info));
  }

  int nearest_pressure_dof(const FEMSpace<double, dim> &pressure_space,
                           const double x,
                           const double y)
  {
    int    nearest = 0;
    double nearest_distance_square = std::numeric_limits<double>::max();
    for (int i = 0; i < pressure_space.n_dof(); ++i)
      {
        const auto &point = pressure_space.dofInfo(i).interp_point;
        const double dx = point[0] - x;
        const double dy = point[1] - y;
        const double distance_square = dx * dx + dy * dy;
        if (distance_square < nearest_distance_square)
          {
            nearest = i;
            nearest_distance_square = distance_square;
          }
      }
    return nearest;
  }

  double pressure_integral(const FEMFunction<double, dim> &pressure)
  {
    const auto &space = pressure.femSpace();
    double integral = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(4);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto values = pressure.value(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          integral += quad_info.weight(q) * jacobian[q] * volume * values[q];
      }

    return integral;
  }

  double divergence_l2_norm(const FEMFunction<double, dim> &u_x,
                            const FEMFunction<double, dim> &u_y)
  {
    const auto &space = u_x.femSpace();
    double norm_square = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(5);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto grad_x = u_x.gradient(q_point, element);
        const auto grad_y = u_y.gradient(q_point, element);

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double div_u = grad_x[q][0] + grad_y[q][1];
            norm_square += jxw * div_u * div_u;
          }
      }

    return std::sqrt(norm_square);
  }

  CellDivergenceDiagnostic cell_divergence_diagnostic(
    const FEMFunction<double, dim> &u_x,
    const FEMFunction<double, dim> &u_y)
  {
    const auto &space = u_x.femSpace();
    CellDivergenceDiagnostic diagnostic;
    double cell_mean_norm_square = 0.0;

    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(5);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());
        const auto q_point =
          element.local_to_global(quad_info.quadraturePoint());
        const auto grad_x = u_x.gradient(q_point, element);
        const auto grad_y = u_y.gradient(q_point, element);

        double cell_measure = 0.0;
        double divergence_integral = 0.0;
        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double div_u = grad_x[q][0] + grad_y[q][1];
            cell_measure += jxw;
            divergence_integral += jxw * div_u;
          }

        diagnostic.max_integral =
          std::max(diagnostic.max_integral, std::abs(divergence_integral));
        if (cell_measure > 0.0)
          cell_mean_norm_square +=
            divergence_integral * divergence_integral / cell_measure;
      }

    diagnostic.mean_l2_norm = std::sqrt(cell_mean_norm_square);
    return diagnostic;
  }

  double max_speed(const FEMFunction<double, dim> &u_x,
                   const FEMFunction<double, dim> &u_y)
  {
    double maximum = 0.0;
    for (std::size_t i = 0; i < u_x.size(); ++i)
      {
        const double speed =
          std::sqrt(u_x(i) * u_x(i) + u_y(i) * u_y(i));
        maximum = std::max(maximum, speed);
      }
    return maximum;
  }

  double integrate_vertical_samples(
    const std::vector<VerticalProfileSample> &samples)
  {
    if (samples.size() < 2)
      return 0.0;

    if (samples.size() % 2 == 1)
      {
        double rate = 0.0;
        bool   valid_simpson_rule = true;
        for (std::size_t i = 0; i + 2 < samples.size(); i += 2)
          {
            const double dy0 = samples[i + 1].y - samples[i].y;
            const double dy1 = samples[i + 2].y - samples[i + 1].y;
            if (dy0 <= 0.0 ||
                std::abs(dy0 - dy1) >
                  64.0 * boundary_tolerance * std::max(1.0, std::abs(dy0)))
              {
                valid_simpson_rule = false;
                break;
              }
            rate += (dy0 + dy1) / 6.0 *
                    (samples[i].ux + 4.0 * samples[i + 1].ux +
                     samples[i + 2].ux);
          }

        if (valid_simpson_rule)
          return rate;
      }

    double rate = 0.0;
    for (std::size_t i = 1; i < samples.size(); ++i)
      {
        const double dy = samples[i].y - samples[i - 1].y;
        rate += 0.5 * dy * (samples[i].ux + samples[i - 1].ux);
      }
    return rate;
  }

  bool barycentric_coordinates(const Point<dim> &p,
                               const Point<dim> &a,
                               const Point<dim> &b,
                               const Point<dim> &c,
                               double &lambda_a,
                               double &lambda_b,
                               double &lambda_c)
  {
    const double determinant =
      (b[1] - c[1]) * (a[0] - c[0]) +
      (c[0] - b[0]) * (a[1] - c[1]);
    if (std::abs(determinant) < 1.0e-30)
      return false;

    lambda_a = ((b[1] - c[1]) * (p[0] - c[0]) +
                (c[0] - b[0]) * (p[1] - c[1])) /
               determinant;
    lambda_b = ((c[1] - a[1]) * (p[0] - c[0]) +
                (a[0] - c[0]) * (p[1] - c[1])) /
               determinant;
    lambda_c = 1.0 - lambda_a - lambda_b;
    return true;
  }

  bool point_in_triangle(const Point<dim> &point,
                         const Point<dim> &a,
                         const Point<dim> &b,
                         const Point<dim> &c)
  {
    constexpr double tolerance = 1.0e-9;

    double lambda_a = 0.0;
    double lambda_b = 0.0;
    double lambda_c = 0.0;
    if (!barycentric_coordinates(point,
                                 a,
                                 b,
                                 c,
                                 lambda_a,
                                 lambda_b,
                                 lambda_c))
      return false;

    return lambda_a >= -tolerance && lambda_b >= -tolerance &&
           lambda_c >= -tolerance && lambda_a <= 1.0 + tolerance &&
           lambda_b <= 1.0 + tolerance && lambda_c <= 1.0 + tolerance;
  }

  bool evaluate_in_element(const FEMFunction<double, dim> &u_x,
                           const Point<dim> &point,
                           double &value)
  {
    const auto &space = u_x.femSpace();
    const EasyMesh &mesh = static_cast<const EasyMesh &>(space.mesh());

    for (unsigned int cell = 0; cell < mesh.n_geometry(dim); ++cell)
      {
        const GeometryBM &triangle = mesh.geometry(dim, cell);
        if (triangle.n_vertex() != 3)
          continue;

        const Point<dim> &a = mesh.point(triangle.vertex(0));
        const Point<dim> &b = mesh.point(triangle.vertex(1));
        const Point<dim> &c = mesh.point(triangle.vertex(2));
        if (!point_in_triangle(point, a, b, c))
          continue;

        value = u_x.value(point, space.element(static_cast<int>(cell)));
        return true;
      }

    return false;
  }

  VerticalProfile vertical_line_profile(const FEMFunction<double, dim> &u_x,
                                        const double x)
  {
    VerticalProfile profile;
    profile.x = x;
    profile.y_min = channel_y_min(x);
    profile.y_max = channel_y_max(x);

    const int intervals =
      std::max(2, static_cast<int>(std::lround(
                    (profile.y_max - profile.y_min) *
                    profile_sampling_density())));
    for (int i = 0; i <= intervals; ++i)
      {
        Point<dim> point;
        point[0] = x;
        point[1] = profile.y_min +
                   (profile.y_max - profile.y_min) * i / intervals;

        double value = 0.0;
        if (evaluate_in_element(u_x, point, value))
          profile.samples.push_back({point[1], value});
      }

    profile.flow_rate = integrate_vertical_samples(profile.samples);

    if (!profile.samples.empty())
      {
        profile.ux_min = profile.samples.front().ux;
        profile.ux_max = profile.samples.front().ux;
        profile.ux_center = profile.samples.front().ux;

        double nearest_center_distance = std::numeric_limits<double>::max();
        const double center = 0.5 * (profile.y_min + profile.y_max);
        for (const auto &sample : profile.samples)
          {
            profile.ux_min = std::min(profile.ux_min, sample.ux);
            profile.ux_max = std::max(profile.ux_max, sample.ux);
            const double distance = std::abs(sample.y - center);
            if (distance < nearest_center_distance)
              {
                nearest_center_distance = distance;
                profile.ux_center = sample.ux;
              }
          }
      }

    return profile;
  }

  std::vector<VerticalProfile> contracting_channel_profiles(
    const FEMFunction<double, dim> &u_x)
  {
    return {vertical_line_profile(u_x, 0.0),
            vertical_line_profile(u_x, 1.0),
            vertical_line_profile(u_x, 4.0)};
  }

  ContractingChannelResult solve_contracting_channel(
    const std::string &mesh_file,
    const std::string &output_prefix,
    const int velocity_order,
    const int pressure_order,
    const bool stabilized,
    const bool discontinuous_pressure = false,
    const bool pin_pressure = true)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> velocity_geometry;
    CoordTransform<dim, dim> velocity_transform;
    TemplateDOF<dim> velocity_dof;
    BasisFunctionAdmin<double, dim, dim> velocity_basis(velocity_dof);
    auto velocity_template = make_triangle_template(velocity_order,
                                                    velocity_geometry,
                                                    velocity_transform,
                                                    velocity_dof,
                                                    velocity_basis);

    TemplateGeometry<dim> pressure_geometry;
    CoordTransform<dim, dim> pressure_transform;
    TemplateDOF<dim> pressure_dof;
    BasisFunctionAdmin<double, dim, dim> pressure_basis(pressure_dof);
    const std::string pressure_template_name =
      discontinuous_pressure ? std::to_string(pressure_order) + ".D" :
                               std::to_string(pressure_order);
    auto pressure_template =
      make_triangle_template_by_name(pressure_template_name,
                                     pressure_geometry,
                                     pressure_transform,
                                     pressure_dof,
                                     pressure_basis);

    auto velocity_space = build_space(mesh, velocity_template);
    auto pressure_space = build_space(mesh, pressure_template);

    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;
    const int n_system = n_total;

    std::vector<double> matrix(static_cast<std::size_t>(n_system) * n_system,
                               0.0);
    std::vector<double> rhs(n_system, 0.0);

    for (int e = 0; e < velocity_space.n_element(); ++e)
      {
        const auto &velocity_element = velocity_space.element(e);
        const auto &pressure_element = pressure_space.element(e);
        const double volume = velocity_element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          velocity_element.findQuadratureInfo(5);
        const auto jacobian =
          velocity_element.local_to_global_jacobian(
            quad_info.quadraturePoint());
        const auto q_point =
          velocity_element.local_to_global(quad_info.quadraturePoint());
        const auto velocity_grads =
          velocity_element.basis_function_gradient(q_point);
        const auto pressure_values =
          pressure_element.basis_function_value(q_point);

        const auto &velocity_element_dofs = velocity_element.dof();
        const auto &pressure_element_dofs = pressure_element.dof();
        const int n_local_velocity_dofs = velocity_element_dofs.size();
        const int n_local_pressure_dofs = pressure_element_dofs.size();

        double element_measure = 0.0;
        std::vector<double> pressure_mean(n_local_pressure_dofs, 0.0);
        if (stabilized)
          {
            for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
              {
                const double jxw =
                  quad_info.weight(q) * jacobian[q] * volume;
                element_measure += jxw;
                for (int i = 0; i < n_local_pressure_dofs; ++i)
                  pressure_mean[i] += jxw * pressure_values[i][q];
              }
            for (double &mean : pressure_mean)
              mean /= element_measure;
          }

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    const double stiffness =
                      viscosity * jxw *
                      innerProduct(velocity_grads[i][q],
                                   velocity_grads[j][q]);
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_system,
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                stiffness);
                  }

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj = pressure_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_system,
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                pressure_index(pj, n_velocity_dofs),
                                -jxw * pressure_values[j][q] *
                                  velocity_grads[i][q][component]);
                  }
              }

            for (int i = 0; i < n_local_pressure_dofs; ++i)
              {
                const int pi_dof = pressure_element_dofs[i];
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_system,
                                pressure_index(pi_dof, n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                -jxw * pressure_values[i][q] *
                                  velocity_grads[j][q][component]);
                  }

                if (stabilized)
                  for (int j = 0; j < n_local_pressure_dofs; ++j)
                    {
                      const int pj = pressure_element_dofs[j];
                      const double stabilized_i =
                        pressure_values[i][q] - pressure_mean[i];
                      const double stabilized_j =
                        pressure_values[j][q] - pressure_mean[j];
                      add_entry(matrix,
                                n_system,
                                pressure_index(pi_dof, n_velocity_dofs),
                                pressure_index(pj, n_velocity_dofs),
                                -stabilization_parameter * jxw *
                                  stabilized_i * stabilized_j);
                    }
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (is_velocity_dirichlet_dof(velocity_space, i))
        {
          const auto &point = velocity_space.dofInfo(i).interp_point;
          const double p[dim] = {point[0], point[1]};
          for (int component = 0; component < dim; ++component)
            apply_dirichlet_constraint(matrix,
                                       rhs,
                                       n_system,
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs),
                                       velocity_boundary_value(component, p));
        }

    if (pin_pressure)
      {
        const int pressure_pin = nearest_pressure_dof(pressure_space,
                                                      outlet_x,
                                                      0.0);
        apply_dirichlet_constraint(
          matrix,
          rhs,
          n_system,
          pressure_index(pressure_pin, n_velocity_dofs),
          0.0);
      }

    solve_dense_system(matrix, rhs, n_system);

    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);
    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        u_x(i) = rhs[velocity_index(0, i, n_velocity_dofs)];
        u_y(i) = rhs[velocity_index(1, i, n_velocity_dofs)];
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      p_h(i) = rhs[pressure_index(i, n_velocity_dofs)];

    u_x.writeOpenDXData(output_prefix + "_ux.dx");
    u_y.writeOpenDXData(output_prefix + "_uy.dx");
    p_h.writeOpenDXData(output_prefix + "_p.dx");

    const double ux_l2 = Functional::L2Norm(u_x, 5);
    const double uy_l2 = Functional::L2Norm(u_y, 5);
    const double ux_h1 = Functional::H1Seminorm(u_x, 5);
    const double uy_h1 = Functional::H1Seminorm(u_y, 5);

    double pressure_min = std::numeric_limits<double>::max();
    double pressure_max = -std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < p_h.size(); ++i)
      {
        pressure_min = std::min(pressure_min, p_h(i));
        pressure_max = std::max(pressure_max, p_h(i));
      }

    const auto profiles = contracting_channel_profiles(u_x);
    const double inlet_rate = profiles[0].flow_rate;
    const double outlet_rate = profiles[2].flow_rate;
    const auto cell_divergence =
      cell_divergence_diagnostic(u_x, u_y);

    return {static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            max_speed(u_x, u_y),
            inlet_rate,
            outlet_rate,
            outlet_rate - inlet_rate,
            pressure_min,
            pressure_max,
            pressure_integral(p_h) / channel_area(),
            divergence_l2_norm(u_x, u_y),
            cell_divergence.max_integral,
            cell_divergence.mean_l2_norm,
            profiles};
  }

  void print_profile_summary(const std::string &method,
                             const ContractingChannelResult &result)
  {
    for (const auto &profile : result.profiles)
      std::cout << method << ' '
                << profile.x << ' '
                << profile.y_min << ' '
                << profile.y_max << ' '
                << profile.samples.size() << ' '
                << profile.flow_rate << ' '
                << profile.flow_rate - reference_flow_rate() << ' '
                << profile.ux_min << ' '
                << profile.ux_center << ' '
                << profile.ux_max << '\n';
  }

  void print_profile_samples(const std::string &method,
                             const ContractingChannelResult &result)
  {
    for (const auto &profile : result.profiles)
      for (const auto &sample : profile.samples)
        std::cout << method << ' '
                  << profile.x << ' '
                  << sample.y << ' '
                  << sample.ux << '\n';
  }
}

int main(int argc, char **argv)
{
  if (argc > 3)
    {
      std::cerr << "usage: " << argv[0] << " [mesh-file [output-prefix]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) +
                 "/contracting_channel_h0p25";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/stokes_contracting_channel_profiles_afepack";

  try
    {
      const auto taylor_hood =
        solve_contracting_channel(mesh_file,
                                  output_prefix + "_p2p1",
                                  2,
                                  1,
                                  false);
      const auto equal_order =
        solve_contracting_channel(mesh_file,
                                  output_prefix + "_p1p1_lps",
                                  1,
                                  1,
                                  true);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Stokes contracting-channel profile diagnostic\n";
      std::cout << "Book Exercise 5.7 proxy on the Exercise 1.5 bow-tie "
                   "contraction: u_x profiles at x=0, x=1, and x=4.\n";
      std::cout << "reference_inlet_flow_rate "
                << reference_flow_rate() << '\n';
      std::cout << "summary method x y_min y_max samples flow_rate "
                   "flow_rate_minus_inlet ux_min ux_center ux_max\n";
      print_profile_summary("P2P1_TaylorHood", taylor_hood);
      print_profile_summary("P1P1_LPS", equal_order);

      std::cout << "profile_samples method x y ux\n";
      print_profile_samples("P2P1_TaylorHood", taylor_hood);
      print_profile_samples("P1P1_LPS", equal_order);

      std::cout << "diagnostics method cells velocity_dofs pressure_dofs "
                   "velocity_L2 velocity_H1 max_speed flux_imbalance "
                   "pressure_min pressure_max pressure_mean divergence_L2\n";
      std::cout << "P2P1_TaylorHood "
                << taylor_hood.cells << ' '
                << taylor_hood.velocity_dofs << ' '
                << taylor_hood.pressure_dofs << ' '
                << taylor_hood.velocity_l2_norm << ' '
                << taylor_hood.velocity_h1_seminorm << ' '
                << taylor_hood.max_speed << ' '
                << taylor_hood.flux_imbalance << ' '
                << taylor_hood.pressure_min << ' '
                << taylor_hood.pressure_max << ' '
                << taylor_hood.pressure_mean << ' '
                << taylor_hood.divergence_l2_norm << '\n';
      std::cout << "P1P1_LPS "
                << equal_order.cells << ' '
                << equal_order.velocity_dofs << ' '
                << equal_order.pressure_dofs << ' '
                << equal_order.velocity_l2_norm << ' '
                << equal_order.velocity_h1_seminorm << ' '
                << equal_order.max_speed << ' '
                << equal_order.flux_imbalance << ' '
                << equal_order.pressure_min << ' '
                << equal_order.pressure_max << ' '
                << equal_order.pressure_mean << ' '
                << equal_order.divergence_l2_norm << '\n';
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
