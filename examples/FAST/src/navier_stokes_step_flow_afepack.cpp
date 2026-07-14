/**
 * @file navier_stokes_step_flow_afepack.cpp
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
 *
 * 网格与数据：主要依赖 后台阶通道系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
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
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
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
  constexpr int    max_picard_iterations = 160;
  constexpr double picard_tolerance = 1.0e-8;
  constexpr double boundary_tolerance = 1.0e-12;
  constexpr double exercise_7_4_profile_x = 5.0;

  struct StepFlowResult
  {
    double       viscosity = 0.0;
    double       reynolds_number = 0.0;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    int          picard_iterations = 0;
    double       final_relative_update = 0.0;
    double       velocity_l2_norm = 0.0;
    double       velocity_h1_seminorm = 0.0;
    double       max_speed = 0.0;
    double       inflow_rate = 0.0;
    double       outflow_rate = 0.0;
    double       flux_imbalance = 0.0;
    double       pressure_min = 0.0;
    double       pressure_max = 0.0;
    double       pressure_mean = 0.0;
    double       divergence_l2_norm = 0.0;
    int          wall_shear_samples = 0;
    int          reattachment_found = 0;
    double       reattachment_x = -1.0;
    double       wall_shear_min = 0.0;
    double       wall_shear_max = 0.0;
    double       outlet_x = 5.0;
    double       profile_x = exercise_7_4_profile_x;
    std::vector<double> profile_y;
    std::vector<double> profile_ux;
  };

  struct WallShearDiagnostic
  {
    int    samples = 0;
    int    reattachment_found = 0;
    double reattachment_x = -1.0;
    double shear_min = 0.0;
    double shear_max = 0.0;
  };

  struct VerticalProfile
  {
    double x = exercise_7_4_profile_x;
    std::vector<double> y;
    std::vector<double> ux;
  };

  bool near(const double value, const double target)
  {
    return std::abs(value - target) < boundary_tolerance;
  }

  double inflow_velocity_x(const double *p)
  {
    return 4.0 * p[1] * (1.0 - p[1]);
  }

  double velocity_boundary_value(const int component, const double *p)
  {
    if (component == 0 && near(p[0], -1.0))
      return inflow_velocity_x(p);
    return 0.0;
  }

  double rhs_x(const double *p)
  {
    (void)p;
    return 0.0;
  }

  double rhs_y(const double *p)
  {
    (void)p;
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
  make_triangle_template(const int order,
                         TemplateGeometry<dim> &geometry,
                         CoordTransform<dim, dim> &transform,
                         TemplateDOF<dim> &dof,
                         BasisFunctionAdmin<double, dim, dim> &basis)
  {
    ensure_template_path();

    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof.reinit(geometry);
    dof.readData("triangle." + std::to_string(order) + ".tmp_dof");
    basis.reinit(dof);
    basis.readData("triangle." + std::to_string(order) + ".bas_fun");

    std::vector<TemplateElement<double, dim, dim>> template_element(1);
    template_element[0].reinit(geometry, dof, transform, basis);
    return template_element;
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

  double mesh_outlet_x(const EasyMesh &mesh)
  {
    double outlet_x = -std::numeric_limits<double>::max();
    for (int i = 0; i < mesh.n_geometry(0); ++i)
      outlet_x = std::max(outlet_x, mesh.point(i)[0]);
    return outlet_x;
  }

  bool is_velocity_dirichlet_dof(const FEMSpace<double, dim> &velocity_space,
                                 const int dof,
                                 const double outlet_x)
  {
    if (velocity_space.dofBoundaryMark(dof) == 0)
      return false;

    const auto &point = velocity_space.dofInfo(dof).interp_point;
    const double x = point[0];
    const double y = point[1];

    if (near(x, outlet_x) && !near(y, -1.0) && !near(y, 1.0))
      return false;

    return true;
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

  void assemble_picard_system(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    const FEMFunction<double, dim> &previous_x,
    const FEMFunction<double, dim> &previous_y,
    const double viscosity,
    const double outlet_x,
    std::vector<double> &matrix,
    std::vector<double> &rhs)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    matrix.assign(static_cast<std::size_t>(n_total) * n_total, 0.0);
    rhs.assign(n_total, 0.0);

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
        const auto velocity_values =
          velocity_element.basis_function_value(q_point);
        const auto velocity_grads =
          velocity_element.basis_function_gradient(q_point);
        const auto pressure_values =
          pressure_element.basis_function_value(q_point);
        const auto previous_x_values =
          previous_x.value(q_point, velocity_element);
        const auto previous_y_values =
          previous_y.value(q_point, velocity_element);

        const auto &velocity_element_dofs = velocity_element.dof();
        const auto &pressure_element_dofs = pressure_element.dof();
        const int n_local_velocity_dofs = velocity_element_dofs.size();
        const int n_local_pressure_dofs = pressure_element_dofs.size();

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          {
            const double jxw =
              quad_info.weight(q) * jacobian[q] * volume;
            const double p[dim] = {q_point[q][0], q_point[q][1]};
            const double f[dim] = {rhs_x(p), rhs_y(p)};
            const double wind[dim] = {previous_x_values[q],
                                      previous_y_values[q]};

            for (int i = 0; i < n_local_velocity_dofs; ++i)
              {
                const int vi = velocity_element_dofs[i];
                for (int j = 0; j < n_local_velocity_dofs; ++j)
                  {
                    const int vj = velocity_element_dofs[j];
                    const double diffusion =
                      viscosity * innerProduct(velocity_grads[i][q],
                                               velocity_grads[j][q]);
                    const double convection =
                      (wind[0] * velocity_grads[j][q][0] +
                       wind[1] * velocity_grads[j][q][1]) *
                      velocity_values[i][q];
                    const double velocity_block =
                      jxw * (diffusion + convection);

                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_total,
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                velocity_block);
                  }

                for (int component = 0; component < dim; ++component)
                  rhs[velocity_index(component, vi, n_velocity_dofs)] +=
                    jxw * f[component] * velocity_values[i][q];

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj = pressure_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      add_entry(matrix,
                                n_total,
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
                                n_total,
                                pressure_index(pi_dof, n_velocity_dofs),
                                velocity_index(component,
                                               vj,
                                               n_velocity_dofs),
                                -jxw * pressure_values[i][q] *
                                  velocity_grads[j][q][component]);
                  }
              }
          }
      }

    for (int i = 0; i < n_velocity_dofs; ++i)
      if (is_velocity_dirichlet_dof(velocity_space, i, outlet_x))
        {
          const auto &point = velocity_space.dofInfo(i).interp_point;
          const double p[dim] = {point[0], point[1]};
          for (int component = 0; component < dim; ++component)
            apply_dirichlet_constraint(matrix,
                                       rhs,
                                       n_total,
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs),
                                       velocity_boundary_value(component, p));
        }

    const int pressure_pin =
      nearest_pressure_dof(pressure_space, outlet_x, 0.0);
    apply_dirichlet_constraint(matrix,
                               rhs,
                               n_total,
                               pressure_index(pressure_pin, n_velocity_dofs),
                               0.0);
  }

  void assign_solution(const std::vector<double> &solution,
                       FEMFunction<double, dim> &u_x,
                       FEMFunction<double, dim> &u_y,
                       FEMFunction<double, dim> &p_h)
  {
    const int n_velocity_dofs = u_x.femSpace().n_dof();
    const int n_pressure_dofs = p_h.femSpace().n_dof();

    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        u_x(i) = solution[velocity_index(0, i, n_velocity_dofs)];
        u_y(i) = solution[velocity_index(1, i, n_velocity_dofs)];
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      p_h(i) = solution[pressure_index(i, n_velocity_dofs)];
  }

  void copy_velocity(FEMFunction<double, dim> &destination_x,
                     FEMFunction<double, dim> &destination_y,
                     const FEMFunction<double, dim> &source_x,
                     const FEMFunction<double, dim> &source_y)
  {
    for (std::size_t i = 0; i < destination_x.size(); ++i)
      {
        destination_x(i) = source_x(i);
        destination_y(i) = source_y(i);
      }
  }

  void set_zero(FEMFunction<double, dim> &function)
  {
    for (std::size_t i = 0; i < function.size(); ++i)
      function(i) = 0.0;
  }

  double relative_velocity_update(const FEMFunction<double, dim> &old_x,
                                  const FEMFunction<double, dim> &old_y,
                                  const FEMFunction<double, dim> &new_x,
                                  const FEMFunction<double, dim> &new_y)
  {
    double update_square = 0.0;
    double solution_square = 0.0;

    for (std::size_t i = 0; i < new_x.size(); ++i)
      {
        const double dx = new_x(i) - old_x(i);
        const double dy = new_y(i) - old_y(i);
        update_square += dx * dx + dy * dy;
        solution_square += new_x(i) * new_x(i) + new_y(i) * new_y(i);
      }

    return std::sqrt(update_square) /
           std::max(std::sqrt(solution_square), 1.0e-30);
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

  double domain_area(const FEMSpace<double, dim> &space)
  {
    double area = 0.0;
    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        const double volume = element.templateElement().volume();
        const QuadratureInfo<dim> &quad_info =
          element.findQuadratureInfo(2);
        const auto jacobian =
          element.local_to_global_jacobian(quad_info.quadraturePoint());

        for (int q = 0; q < quad_info.n_quadraturePoint(); ++q)
          area += quad_info.weight(q) * jacobian[q] * volume;
      }
    return area;
  }

  double oriented_area2(const Point<dim> &a,
                        const Point<dim> &b,
                        const Point<dim> &c)
  {
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
  }

  bool point_in_triangle(const Point<dim> &point,
                         const std::vector<Point<dim>> &vertices)
  {
    if (vertices.size() != 3)
      return false;

    const double denominator =
      oriented_area2(vertices[0], vertices[1], vertices[2]);
    if (std::abs(denominator) < 1.0e-30)
      return false;

    const double lambda0 =
      oriented_area2(point, vertices[1], vertices[2]) / denominator;
    const double lambda1 =
      oriented_area2(vertices[0], point, vertices[2]) / denominator;
    const double lambda2 =
      oriented_area2(vertices[0], vertices[1], point) / denominator;
    const double tolerance = 1.0e-10;
    return lambda0 >= -tolerance && lambda1 >= -tolerance &&
           lambda2 >= -tolerance && lambda0 <= 1.0 + tolerance &&
           lambda1 <= 1.0 + tolerance && lambda2 <= 1.0 + tolerance;
  }

  bool value_at_point(const FEMFunction<double, dim> &function,
                      const Point<dim> &point,
                      double &value)
  {
    const auto &space = function.femSpace();
    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        std::vector<Point<dim>> vertices;
        element.buildVertexArray(vertices);
        if (!point_in_triangle(point, vertices))
          continue;

        const std::vector<Point<dim>> points = {point};
        const auto values = function.value(points, element);
        value = values[0];
        return true;
      }
    return false;
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

  VerticalProfile horizontal_velocity_profile(
    const FEMFunction<double, dim> &u_x,
    const double requested_x,
    const double outlet_x)
  {
    VerticalProfile profile;
    profile.x = requested_x;
    profile.y = {-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75};

    if (requested_x > outlet_x + boundary_tolerance)
      return profile;

    const double x =
      near(requested_x, outlet_x) ? outlet_x - 1.0e-10 : requested_x;
    for (const double y : profile.y)
      {
        double value = 0.0;
        const Point<dim> point(x, y);
        if (!value_at_point(u_x, point, value))
          {
            profile.ux.clear();
            return profile;
          }
        profile.ux.push_back(value);
      }

    return profile;
  }

  WallShearDiagnostic lower_wall_shear_diagnostic(
    const FEMFunction<double, dim> &u_x,
    const double outlet_x)
  {
    struct Sample
    {
      double x = 0.0;
      double shear = 0.0;
    };

    std::vector<Sample> samples;
    const auto &space = u_x.femSpace();
    for (int e = 0; e < space.n_element(); ++e)
      {
        const Element<double, dim> &element = space.element(e);
        std::vector<Point<dim>> vertices;
        element.buildVertexArray(vertices);

        std::vector<int> lower_vertices;
        for (std::size_t i = 0; i < vertices.size(); ++i)
          if (near(vertices[i][1], -1.0) &&
              vertices[i][0] > -boundary_tolerance &&
              vertices[i][0] < outlet_x + boundary_tolerance)
            lower_vertices.push_back(static_cast<int>(i));

        if (lower_vertices.size() < 2)
          continue;

        int    first = lower_vertices[0];
        int    second = lower_vertices[1];
        double edge_length = 0.0;
        for (std::size_t i = 0; i < lower_vertices.size(); ++i)
          for (std::size_t j = i + 1; j < lower_vertices.size(); ++j)
            {
              const int a = lower_vertices[i];
              const int b = lower_vertices[j];
              const double length = std::abs(vertices[a][0] - vertices[b][0]);
              if (length > edge_length)
                {
                  edge_length = length;
                  first = a;
                  second = b;
                }
            }

        if (edge_length <= boundary_tolerance)
          continue;

        const double x_mid =
          0.5 * (vertices[first][0] + vertices[second][0]);
        const Point<dim> point(x_mid, -1.0 + 1.0e-10);
        const std::vector<Point<dim>> q_point = {point};
        const auto gradients = u_x.gradient(q_point, element);
        samples.push_back({x_mid, gradients[0][1]});
      }

    WallShearDiagnostic diagnostic;
    diagnostic.samples = static_cast<int>(samples.size());
    if (samples.empty())
      return diagnostic;

    std::sort(samples.begin(),
              samples.end(),
              [](const Sample &left, const Sample &right) {
                return left.x < right.x;
              });

    diagnostic.shear_min = samples.front().shear;
    diagnostic.shear_max = samples.front().shear;
    for (const auto &sample : samples)
      {
        diagnostic.shear_min =
          std::min(diagnostic.shear_min, sample.shear);
        diagnostic.shear_max =
          std::max(diagnostic.shear_max, sample.shear);
      }

    auto interpolate_crossing = [](const Sample &left,
                                   const Sample &right) {
      const double delta_shear = right.shear - left.shear;
      if (std::abs(delta_shear) < 1.0e-30)
        return 0.5 * (left.x + right.x);
      const double theta = -left.shear / delta_shear;
      const double clamped_theta = std::max(0.0, std::min(1.0, theta));
      return left.x + clamped_theta * (right.x - left.x);
    };

    for (std::size_t i = 1; i < samples.size(); ++i)
      if (samples[i - 1].shear <= 0.0 && samples[i].shear >= 0.0)
        {
          diagnostic.reattachment_found = 1;
          diagnostic.reattachment_x =
            interpolate_crossing(samples[i - 1], samples[i]);
          return diagnostic;
        }

    for (std::size_t i = 1; i < samples.size(); ++i)
      if ((samples[i - 1].shear <= 0.0 && samples[i].shear >= 0.0) ||
          (samples[i - 1].shear >= 0.0 && samples[i].shear <= 0.0))
        {
          diagnostic.reattachment_found = 1;
          diagnostic.reattachment_x =
            interpolate_crossing(samples[i - 1], samples[i]);
          return diagnostic;
        }

    return diagnostic;
  }

  double vertical_line_flow_rate(const FEMFunction<double, dim> &u_x,
                                 const double x,
                                 const double y_min,
                                 const double y_max)
  {
    std::vector<std::pair<double, double>> samples;
    const auto &space = u_x.femSpace();
    for (std::size_t i = 0; i < u_x.size(); ++i)
      {
        const auto &point = space.dofInfo(static_cast<int>(i)).interp_point;
        if (near(point[0], x) &&
            point[1] > y_min - boundary_tolerance &&
            point[1] < y_max + boundary_tolerance)
          samples.emplace_back(point[1], u_x(i));
      }

    std::sort(samples.begin(), samples.end());
    if (samples.size() < 2)
      return 0.0;

    if (samples.size() % 2 == 1)
      {
        double rate = 0.0;
        bool   valid_simpson_rule = true;
        for (std::size_t i = 0; i + 2 < samples.size(); i += 2)
          {
            const double dy0 = samples[i + 1].first - samples[i].first;
            const double dy1 = samples[i + 2].first - samples[i + 1].first;
            if (dy0 <= 0.0 ||
                std::abs(dy0 - dy1) >
                  64.0 * boundary_tolerance * std::max(1.0, std::abs(dy0)))
              {
                valid_simpson_rule = false;
                break;
              }
            rate += (dy0 + dy1) / 6.0 *
                    (samples[i].second + 4.0 * samples[i + 1].second +
                     samples[i + 2].second);
          }

        if (valid_simpson_rule)
          return rate;
      }

    double rate = 0.0;
    for (std::size_t i = 1; i < samples.size(); ++i)
      {
        const double dy = samples[i].first - samples[i - 1].first;
        rate += 0.5 * dy * (samples[i].second + samples[i - 1].second);
      }
    return rate;
  }

  std::string viscosity_label(const double viscosity)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << viscosity;
    std::string label = stream.str();
    while (!label.empty() && label.back() == '0')
      label.pop_back();
    if (!label.empty() && label.back() == '.')
      label.pop_back();
    std::replace(label.begin(), label.end(), '.', 'p');
    return label;
  }

  std::vector<double> default_viscosities()
  {
    return {1.0, 0.5, 0.25, 0.125, 0.05, 0.02};
  }

  std::vector<double> parse_viscosities(const std::string &text)
  {
    if (text.empty())
      return default_viscosities();

    std::vector<double> viscosities;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ','))
      {
        const double viscosity = std::stod(token);
        if (!(viscosity > 0.0))
          throw std::runtime_error("viscosity values must be positive");
        viscosities.push_back(viscosity);
      }

    if (viscosities.empty())
      throw std::runtime_error("viscosity list is empty");

    return viscosities;
  }

  StepFlowResult solve_viscosity_step(
    const EasyMesh &mesh,
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    const double viscosity,
    const double outlet_x,
    const double area,
    const std::string &output_prefix,
    FEMFunction<double, dim> &old_x,
    FEMFunction<double, dim> &old_y,
    FEMFunction<double, dim> &u_x,
    FEMFunction<double, dim> &u_y,
    FEMFunction<double, dim> &p_h)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    std::vector<double> matrix;
    std::vector<double> linear_solution;
    double final_relative_update = 1.0;
    int picard_iteration = 0;

    for (; picard_iteration < max_picard_iterations; ++picard_iteration)
      {
        assemble_picard_system(velocity_space,
                               pressure_space,
                               old_x,
                               old_y,
                               viscosity,
                               outlet_x,
                               matrix,
                               linear_solution);

        if (static_cast<int>(linear_solution.size()) != n_total)
          throw std::runtime_error("unexpected linear system size");

        solve_dense_system(matrix, linear_solution, n_total);
        assign_solution(linear_solution, u_x, u_y, p_h);

        final_relative_update =
          relative_velocity_update(old_x, old_y, u_x, u_y);
        ++picard_iteration;

        copy_velocity(old_x, old_y, u_x, u_y);

        if (final_relative_update < picard_tolerance)
          break;
      }

    if (final_relative_update >= picard_tolerance)
      throw std::runtime_error("Picard iteration did not converge for nu=" +
                               std::to_string(viscosity) +
                               "; final relative update=" +
                               std::to_string(final_relative_update));

    const std::string step_prefix =
      output_prefix + "_nu" + viscosity_label(viscosity);
    u_x.writeOpenDXData(step_prefix + "_ux.dx");
    u_y.writeOpenDXData(step_prefix + "_uy.dx");
    p_h.writeOpenDXData(step_prefix + "_p.dx");

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

    const double inflow_rate = vertical_line_flow_rate(u_x, -1.0, 0.0, 1.0);
    const double outflow_rate =
      vertical_line_flow_rate(u_x, outlet_x, -1.0, 1.0);
    const auto wall_shear = lower_wall_shear_diagnostic(u_x, outlet_x);
    const auto profile =
      horizontal_velocity_profile(u_x, exercise_7_4_profile_x, outlet_x);

    return {viscosity,
            2.0 / viscosity,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            picard_iteration,
            final_relative_update,
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            max_speed(u_x, u_y),
            inflow_rate,
            outflow_rate,
            outflow_rate - inflow_rate,
            pressure_min,
            pressure_max,
            pressure_integral(p_h) / area,
            divergence_l2_norm(u_x, u_y),
            wall_shear.samples,
            wall_shear.reattachment_found,
            wall_shear.reattachment_x,
            wall_shear.shear_min,
            wall_shear.shear_max,
            outlet_x,
            profile.x,
            profile.y,
            profile.ux};
  }

  std::vector<StepFlowResult> solve_step_flow_continuation(
    const std::string &mesh_file,
    const std::string &output_prefix,
    const std::vector<double> &viscosities)
  {
    EasyMesh mesh;
    mesh.readData(mesh_file);

    TemplateGeometry<dim> velocity_geometry;
    CoordTransform<dim, dim> velocity_transform;
    TemplateDOF<dim> velocity_dof;
    BasisFunctionAdmin<double, dim, dim> velocity_basis(velocity_dof);
    auto velocity_template = make_triangle_template(2,
                                                    velocity_geometry,
                                                    velocity_transform,
                                                    velocity_dof,
                                                    velocity_basis);

    TemplateGeometry<dim> pressure_geometry;
    CoordTransform<dim, dim> pressure_transform;
    TemplateDOF<dim> pressure_dof;
    BasisFunctionAdmin<double, dim, dim> pressure_basis(pressure_dof);
    auto pressure_template = make_triangle_template(1,
                                                    pressure_geometry,
                                                    pressure_transform,
                                                    pressure_dof,
                                                    pressure_basis);

    auto velocity_space = build_space(mesh, velocity_template);
    auto pressure_space = build_space(mesh, pressure_template);
    const double outlet_x = mesh_outlet_x(mesh);
    const double area = domain_area(pressure_space);

    FEMFunction<double, dim> old_x(velocity_space);
    FEMFunction<double, dim> old_y(velocity_space);
    FEMFunction<double, dim> u_x(velocity_space);
    FEMFunction<double, dim> u_y(velocity_space);
    FEMFunction<double, dim> p_h(pressure_space);

    set_zero(old_x);
    set_zero(old_y);
    set_zero(u_x);
    set_zero(u_y);
    set_zero(p_h);

    std::vector<StepFlowResult> results;
    results.reserve(viscosities.size());

    for (const double viscosity : viscosities)
      results.push_back(solve_viscosity_step(mesh,
                                             velocity_space,
                                             pressure_space,
                                             viscosity,
                                             outlet_x,
                                             area,
                                             output_prefix,
                                             old_x,
                                             old_y,
                                             u_x,
                                             u_y,
                                             p_h));

    return results;
  }
}

int main(int argc, char **argv)
{
  if (argc > 4)
    {
      std::cerr << "usage: " << argv[0]
                << " [mesh-file [output-prefix [viscosities]]]\n";
      return 2;
    }

  const std::string mesh_file =
    argc > 1 ? argv[1] :
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/step_channel_h0p25";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/navier_stokes_step_flow_afepack";
  const auto viscosities =
    parse_viscosities(argc > 3 ? argv[3] : std::string());

  try
    {
      const auto results =
        solve_step_flow_continuation(mesh_file, output_prefix, viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes expansion step-flow "
                   "continuation (book Example 7.1.2)\n";
      for (const auto &result : results)
        {
          std::cout << "nu: " << result.viscosity
                    << " book Re=2/nu: " << result.reynolds_number << '\n';
          std::cout << "cells: " << result.cells << '\n';
          std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
          std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
          std::cout << "outlet x: " << result.outlet_x << '\n';
          std::cout << "Picard iterations: "
                    << result.picard_iterations << '\n';
          std::cout << "final relative velocity update: "
                    << result.final_relative_update << '\n';
          std::cout << "velocity L2 norm: "
                    << result.velocity_l2_norm << '\n';
          std::cout << "velocity H1-semi norm: "
                    << result.velocity_h1_seminorm << '\n';
          std::cout << "max nodal speed: " << result.max_speed << '\n';
          std::cout << "inflow/outflow rates: " << result.inflow_rate << " "
                    << result.outflow_rate << '\n';
          std::cout << "outflow minus inflow: "
                    << result.flux_imbalance << '\n';
          std::cout << "pressure min/max: " << result.pressure_min << " "
                    << result.pressure_max << '\n';
          std::cout << "pressure mean: " << result.pressure_mean << '\n';
          std::cout << "divergence L2 norm: "
                    << result.divergence_l2_norm << '\n';
          std::cout << "lower-wall shear samples: "
                    << result.wall_shear_samples << '\n';
          std::cout << "lower-wall shear min/max: "
                    << result.wall_shear_min << " "
                    << result.wall_shear_max << '\n';
          std::cout << "lower-wall reattachment found/x: "
                    << result.reattachment_found << " "
                    << result.reattachment_x << '\n';
          if (!result.profile_ux.empty())
            {
              std::cout << "ux profile x=5 y:";
              for (const double y : result.profile_y)
                std::cout << " " << y;
              std::cout << '\n';
              std::cout << "ux profile x=5 values:";
              for (const double ux : result.profile_ux)
                std::cout << " " << ux;
              std::cout << '\n';
            }
        }
      std::cout << "Wrote " << output_prefix << "_nu*_ux.dx, "
                << output_prefix << "_nu*_uy.dx, and "
                << output_prefix << "_nu*_p.dx\n";
    }
  catch (const std::exception &exc)
    {
      std::cerr << "Exception: " << exc.what() << '\n';
      return 1;
    }

  return 0;
}
