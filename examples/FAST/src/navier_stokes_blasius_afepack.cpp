/**
 * @file navier_stokes_blasius_afepack.cpp
 * @brief AFEPack 不可压 Navier-Stokes 方程算例：Blasius 边界层代理。
 *
 * @details
 * 本文件是 FAST 目录中的 不可压 Navier-Stokes 方程 迁移算例，关注非线性不可压流动问题的速度-压力耦合、线性化以及块预条件求解。
 * 算例默认使用 Makefile 提供的 AFEPACK_DEFAULT_MESH_DIR 与 AFEPACK_DEFAULT_OUTPUT_DIR，
 * 因而可直接读取 examples/FAST/build/meshes/afepack 下由 EasyMesh 生成的网格，
 * 并把数据文件、剖面文件或可视化结果写入 examples/FAST/build。
 *
 * 主要演示内容：
 * - Blasius 边界层代理：使用不可压 Navier-Stokes 方程构造平板边界层的数值代理模型。
 *
 * 网格与数据：主要依赖 平板边界层系列网格；执行 make -C examples/FAST meshes 可重新生成网格文件。
 * 编译与运行：执行 make -C examples/FAST 目标名 编译单个算例，
 * 或执行 make -C examples/FAST run-目标名 在配置好的库路径和模板路径下运行。
 *
 * @note 这是从 ~/Projects/FEM/examples/afepack 复制到本仓库 examples/FAST/src 的迁移副本；
 *       注释只服务于本项目中的 FAST 示例文档化，原始 FEM 算例未被修改。
 */

#include <AFEPack/AMGSolver.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

#include "navier_stokes_schur_common.h"

#ifdef FEM_HAVE_EIGEN
#  include <Eigen/SparseCore>
#  include <Eigen/SparseLU>
#endif

#ifdef FEM_HAVE_UMFPACK
#  include <Eigen/UmfPackSupport>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
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
  using afepack_examples::navier_stokes_schur::RowSparseMatrix;
  using afepack_examples::navier_stokes_schur::SchurMassPreconditioner;
  using afepack_examples::navier_stokes_schur::SolveInfo;
  using afepack_examples::navier_stokes_schur::add_scaled;
  using afepack_examples::navier_stokes_schur::apply_dirichlet_constraint;
  using afepack_examples::navier_stokes_schur::apply_homogeneous_constraint;
  using afepack_examples::navier_stokes_schur::norm;
  using afepack_examples::navier_stokes_schur::solve_right_preconditioned_gmres;

  constexpr int    dim = 2;
  constexpr int    default_max_picard_iterations = 120;
  constexpr double picard_tolerance = 1.0e-8;
  constexpr double boundary_tolerance = 1.0e-12;
  constexpr double blasius_domain_area = 11.6;
  constexpr double free_stream_speed = 1.0;
  constexpr double plate_length = 5.0;
  constexpr double plate_half_thickness = 0.04;
  double           active_plate_half_thickness = plate_half_thickness;
  double           active_blasius_domain_area = blasius_domain_area;

  struct Node
  {
    double x = 0.0;
    double y = 0.0;
    int    boundary_mark = 0;
  };

  struct ElementData
  {
    std::array<int, 3> node = {};
    std::array<int, 3> neighbor = {{-1, -1, -1}};
    std::array<int, 3> side = {{-1, -1, -1}};
  };

  struct SideData
  {
    int a = -1;
    int b = -1;
    int ea = -1;
    int eb = -1;
    int boundary_mark = 0;
  };

  struct GeneratedMesh
  {
    std::vector<Node>        nodes;
    std::vector<ElementData> elements;
    std::vector<SideData>    sides;
    double                   min_edge = 0.0;
    double                   max_edge = 0.0;
  };

  struct BlasiusResult
  {
    double       viscosity = 0.0;
    double       reynolds_number = 0.0;
    int          cells = 0;
    unsigned int velocity_dofs = 0;
    unsigned int pressure_dofs = 0;
    unsigned int system_nonzeros = 0;
    int          nonlinear_iterations = 0;
    int          total_gmres_iterations = 0;
    double       final_relative_update = 0.0;
    double       final_gmres_relative_residual = 0.0;
    double       final_nonlinear_relative_residual = 0.0;
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
    double       upper_layer_x1 = -1.0;
    double       lower_layer_x1 = -1.0;
    double       upper_layer_x4 = -1.0;
    double       lower_layer_x4 = -1.0;
  };

  enum class BlasiusLinearSolver
  {
    schur_gmres,
    eigen_sparse_lu,
    eigen_umfpack
  };

  enum class BlasiusNonlinearSolver
  {
    picard,
    newton
  };

  bool near(const double value, const double target)
  {
    return std::abs(value - target) < boundary_tolerance;
  }

  BlasiusLinearSolver blasius_linear_solver()
  {
    const char *environment = std::getenv("BLASIUS_LINEAR_SOLVER");
    if (environment == nullptr)
      return BlasiusLinearSolver::schur_gmres;

    const std::string value(environment);
    if (value == "schur_gmres")
      return BlasiusLinearSolver::schur_gmres;
    if (value == "eigen_sparse_lu")
      return BlasiusLinearSolver::eigen_sparse_lu;
    if (value == "eigen_umfpack")
      return BlasiusLinearSolver::eigen_umfpack;
    throw std::runtime_error(
      "BLASIUS_LINEAR_SOLVER must be schur_gmres, eigen_sparse_lu, or "
      "eigen_umfpack");
  }

  const char *blasius_linear_solver_name(const BlasiusLinearSolver solver)
  {
    if (solver == BlasiusLinearSolver::schur_gmres)
      return "schur_gmres";
    if (solver == BlasiusLinearSolver::eigen_sparse_lu)
      return "eigen_sparse_lu";
    return "eigen_umfpack";
  }

  BlasiusNonlinearSolver blasius_nonlinear_solver()
  {
    const char *environment = std::getenv("BLASIUS_NONLINEAR_SOLVER");
    if (environment == nullptr)
      return BlasiusNonlinearSolver::picard;

    const std::string value(environment);
    if (value == "picard")
      return BlasiusNonlinearSolver::picard;
    if (value == "newton")
      return BlasiusNonlinearSolver::newton;
    throw std::runtime_error(
      "BLASIUS_NONLINEAR_SOLVER must be picard or newton");
  }

  const char *
  blasius_nonlinear_solver_name(const BlasiusNonlinearSolver solver)
  {
    return solver == BlasiusNonlinearSolver::picard ? "picard" : "newton";
  }

#ifdef FEM_HAVE_EIGEN
  using EigenSparseMatrix = Eigen::SparseMatrix<double>;

  EigenSparseMatrix
  make_eigen_sparse_matrix(const RowSparseMatrix &matrix)
  {
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(matrix.nonzeros());
    for (int row = 0; row < matrix.size(); ++row)
      for (const auto &[column, value] : matrix.row_entries(row))
        entries.emplace_back(row, column, value);

    EigenSparseMatrix eigen_matrix(matrix.size(), matrix.size());
    eigen_matrix.setFromTriplets(entries.begin(), entries.end());
    eigen_matrix.makeCompressed();
    return eigen_matrix;
  }

  SolveInfo finish_eigen_direct_solve(
    const RowSparseMatrix &matrix,
    std::vector<double> &solution,
    const std::vector<double> &rhs,
    const Eigen::VectorXd &eigen_solution)
  {
    if (!eigen_solution.allFinite())
      throw std::runtime_error("Eigen direct solve produced non-finite values");

    solution.assign(eigen_solution.data(),
                    eigen_solution.data() + eigen_solution.size());
    std::vector<double> residual = rhs;
    add_scaled(residual, -1.0, matrix.vmult(solution));
    const double relative_residual =
      norm(residual) / std::max(norm(rhs), 1.0e-30);
    return {std::isfinite(relative_residual) &&
              relative_residual <= 1.0e-10,
            0,
            relative_residual};
  }
#endif

  SolveInfo solve_eigen_sparse_lu(const RowSparseMatrix &matrix,
                                  std::vector<double> &solution,
                                  const std::vector<double> &rhs)
  {
#ifdef FEM_HAVE_EIGEN
    const EigenSparseMatrix eigen_matrix = make_eigen_sparse_matrix(matrix);
    Eigen::SparseLU<EigenSparseMatrix, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(eigen_matrix);
    solver.factorize(eigen_matrix);
    if (solver.info() != Eigen::Success)
      throw std::runtime_error("Eigen SparseLU factorization failed");

    Eigen::Map<const Eigen::VectorXd> eigen_rhs(rhs.data(), rhs.size());
    const Eigen::VectorXd eigen_solution = solver.solve(eigen_rhs);
    if (solver.info() != Eigen::Success)
      throw std::runtime_error("Eigen SparseLU solve failed");
    return finish_eigen_direct_solve(matrix,
                                     solution,
                                     rhs,
                                     eigen_solution);
#else
    static_cast<void>(matrix);
    static_cast<void>(solution);
    static_cast<void>(rhs);
    throw std::runtime_error(
      "BLASIUS_LINEAR_SOLVER=eigen_sparse_lu requires Eigen3 at configure time");
#endif
  }

  class EigenUmfpackSolverCache
  {
  public:
    SolveInfo solve(const RowSparseMatrix &matrix,
                    std::vector<double> &solution,
                    const std::vector<double> &rhs)
    {
#ifdef FEM_HAVE_UMFPACK
      const EigenSparseMatrix eigen_matrix = make_eigen_sparse_matrix(matrix);
      if (!same_pattern(eigen_matrix))
        {
          solver.analyzePattern(eigen_matrix);
          if (solver.info() != Eigen::Success)
            throw std::runtime_error("Eigen UMFPACK analysis failed");
          remember_pattern(eigen_matrix);
        }
      solver.factorize(eigen_matrix);
      if (solver.info() != Eigen::Success)
        throw std::runtime_error("Eigen UMFPACK factorization failed");

      Eigen::Map<const Eigen::VectorXd> eigen_rhs(rhs.data(), rhs.size());
      const Eigen::VectorXd eigen_solution = solver.solve(eigen_rhs);
      if (solver.info() != Eigen::Success)
        throw std::runtime_error("Eigen UMFPACK solve failed");
      return finish_eigen_direct_solve(matrix,
                                       solution,
                                       rhs,
                                       eigen_solution);
#else
      static_cast<void>(matrix);
      static_cast<void>(solution);
      static_cast<void>(rhs);
      throw std::runtime_error(
        "BLASIUS_LINEAR_SOLVER=eigen_umfpack requires Eigen3 and UMFPACK "
        "at configure time");
#endif
    }

  private:
#ifdef FEM_HAVE_UMFPACK
    bool same_pattern(const EigenSparseMatrix &matrix) const
    {
      if (!pattern_initialized ||
          outer_indices.size() !=
            static_cast<std::size_t>(matrix.outerSize() + 1) ||
          inner_indices.size() !=
            static_cast<std::size_t>(matrix.nonZeros()))
        return false;

      return std::equal(outer_indices.begin(),
                        outer_indices.end(),
                        matrix.outerIndexPtr()) &&
             std::equal(inner_indices.begin(),
                        inner_indices.end(),
                        matrix.innerIndexPtr());
    }

    void remember_pattern(const EigenSparseMatrix &matrix)
    {
      outer_indices.assign(matrix.outerIndexPtr(),
                           matrix.outerIndexPtr() + matrix.outerSize() + 1);
      inner_indices.assign(matrix.innerIndexPtr(),
                           matrix.innerIndexPtr() + matrix.nonZeros());
      pattern_initialized = true;
    }

    Eigen::UmfPackLU<EigenSparseMatrix> solver;
    std::vector<int> outer_indices;
    std::vector<int> inner_indices;
    bool pattern_initialized = false;
#endif
  };

  double blasius_gmres_tolerance()
  {
    const char *environment = std::getenv("BLASIUS_GMRES_TOLERANCE");
    if (environment == nullptr)
      return 1.0e-9;
    const double tolerance = std::stod(environment);
    if (!(tolerance > 0.0))
      throw std::runtime_error(
        "BLASIUS_GMRES_TOLERANCE must be positive");
    return tolerance;
  }

  int blasius_gmres_max_iterations()
  {
    const char *environment = std::getenv("BLASIUS_GMRES_MAX_ITERATIONS");
    if (environment == nullptr)
      return 6000;
    const int max_iterations = std::stoi(environment);
    if (max_iterations < 1)
      throw std::runtime_error(
        "BLASIUS_GMRES_MAX_ITERATIONS must be positive");
    return max_iterations;
  }

  double blasius_picard_relaxation()
  {
    const char *environment = std::getenv("BLASIUS_PICARD_RELAXATION");
    if (environment == nullptr)
      return 1.0;
    const double relaxation = std::stod(environment);
    if (!(relaxation > 0.0 && relaxation <= 1.0))
      throw std::runtime_error(
        "BLASIUS_PICARD_RELAXATION must be in (0,1]");
    return relaxation;
  }

  int blasius_picard_max_iterations()
  {
    const char *environment = std::getenv("BLASIUS_PICARD_MAX_ITERATIONS");
    if (environment == nullptr)
      return default_max_picard_iterations;
    const int max_iterations = std::stoi(environment);
    if (max_iterations < 1)
      throw std::runtime_error(
        "BLASIUS_PICARD_MAX_ITERATIONS must be positive");
    return max_iterations;
  }

  void set_blasius_plate_geometry(const double half_thickness,
                                  const double domain_area)
  {
    if (half_thickness < 0.0)
      throw std::runtime_error("Blasius plate half-thickness must be nonnegative");
    if (!(domain_area > 0.0))
      throw std::runtime_error("Blasius domain area must be positive");
    active_plate_half_thickness = half_thickness;
    active_blasius_domain_area = domain_area;
  }

  double edge_length(const Node &a, const Node &b)
  {
    return std::hypot(a.x - b.x, a.y - b.y);
  }

  int blasius_boundary_mark(const double x, const double y)
  {
    if (near(x, -1.0))
      return 1;
    if (near(y, -1.0))
      return 2;
    if (near(x, 5.0))
      return 3;
    if (near(y, 1.0))
      return 4;

    const double half_thickness = active_plate_half_thickness;
    if (half_thickness <= boundary_tolerance)
      return x > -boundary_tolerance &&
             x < 5.0 + boundary_tolerance && near(y, 0.0) ? 6 : 0;

    if (x > -boundary_tolerance && x < 5.0 + boundary_tolerance &&
        near(std::abs(y), half_thickness))
      return 6;
    if (near(x, 0.0) &&
        std::abs(y) < half_thickness + boundary_tolerance)
      return 6;
    return 0;
  }

  std::pair<long long, long long>
  coordinate_key(const double x, const double y)
  {
    constexpr double scale = 1.0e12;
    return {static_cast<long long>(std::llround(scale * x)),
            static_cast<long long>(std::llround(scale * y))};
  }

  int add_generated_node(
    GeneratedMesh &mesh,
    std::map<std::pair<long long, long long>, int> &node_id,
    const double x,
    const double y)
  {
    const auto key = coordinate_key(x, y);
    const auto it = node_id.find(key);
    if (it != node_id.end())
      return it->second;

    const int id = static_cast<int>(mesh.nodes.size());
    mesh.nodes.push_back({x, y, blasius_boundary_mark(x, y)});
    node_id.emplace(key, id);
    return id;
  }

  void add_rectangular_patch(
    GeneratedMesh &mesh,
    std::map<std::pair<long long, long long>, int> &node_id,
    const std::vector<double> &x_coordinates,
    const std::vector<double> &y_coordinates)
  {
    for (std::size_t j = 0; j + 1 < y_coordinates.size(); ++j)
      for (std::size_t i = 0; i + 1 < x_coordinates.size(); ++i)
        {
          const int n00 = add_generated_node(mesh,
                                             node_id,
                                             x_coordinates[i],
                                             y_coordinates[j]);
          const int n10 = add_generated_node(mesh,
                                             node_id,
                                             x_coordinates[i + 1],
                                             y_coordinates[j]);
          const int n11 = add_generated_node(mesh,
                                             node_id,
                                             x_coordinates[i + 1],
                                             y_coordinates[j + 1]);
          const int n01 = add_generated_node(mesh,
                                             node_id,
                                             x_coordinates[i],
                                             y_coordinates[j + 1]);

          ElementData lower;
          lower.node = {{n00, n10, n11}};
          mesh.elements.push_back(lower);

          ElementData upper;
          upper.node = {{n00, n11, n01}};
          mesh.elements.push_back(upper);
        }
  }

  GeneratedMesh generate_stretched_blasius_mesh()
  {
    const std::vector<double> upstream_x = {-1.0, -0.5, 0.0};
    const std::vector<double> plate_x =
      {0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
    const std::vector<double> lower_y =
      {-1.0, -0.7, -0.45, -0.25, -0.14, -0.08,
       -plate_half_thickness};
    const std::vector<double> upper_y =
      {plate_half_thickness, 0.08, 0.14, 0.25, 0.45, 0.7, 1.0};
    const std::vector<double> full_y =
      {-1.0, -0.7, -0.45, -0.25, -0.14, -0.08,
       -plate_half_thickness, plate_half_thickness,
       0.08, 0.14, 0.25, 0.45, 0.7, 1.0};

    GeneratedMesh mesh;
    std::map<std::pair<long long, long long>, int> node_id;
    add_rectangular_patch(mesh, node_id, upstream_x, full_y);
    add_rectangular_patch(mesh, node_id, plate_x, lower_y);
    add_rectangular_patch(mesh, node_id, plate_x, upper_y);

    std::map<std::pair<int, int>, int> edge_to_side;
    auto register_edge = [&](const int element_index,
                             const int local_side,
                             const int a,
                             const int b) {
      const std::pair<int, int> key = std::minmax(a, b);
      auto [it, inserted] =
        edge_to_side.emplace(key, static_cast<int>(mesh.sides.size()));
      if (inserted)
        {
          SideData side;
          side.a = key.first;
          side.b = key.second;
          side.ea = element_index;
          mesh.sides.push_back(side);
        }
      else
        {
          SideData &side = mesh.sides[it->second];
          if (side.eb >= 0)
            throw std::runtime_error("nonmanifold generated Blasius side");
          side.eb = element_index;
        }
      mesh.elements[element_index].side[local_side] = it->second;
    };

    for (int e = 0; e < static_cast<int>(mesh.elements.size()); ++e)
      {
        const auto node = mesh.elements[e].node;
        register_edge(e, 0, node[1], node[2]);
        register_edge(e, 1, node[2], node[0]);
        register_edge(e, 2, node[0], node[1]);
      }

    for (int e = 0; e < static_cast<int>(mesh.elements.size()); ++e)
      for (int s = 0; s < 3; ++s)
        {
          const SideData &side = mesh.sides[mesh.elements[e].side[s]];
          mesh.elements[e].neighbor[s] =
            side.ea == e ? side.eb : side.ea;
        }

    mesh.min_edge = std::numeric_limits<double>::max();
    mesh.max_edge = 0.0;
    for (SideData &side : mesh.sides)
      {
        if (side.eb < 0)
          {
            const double midpoint_x =
              0.5 * (mesh.nodes[side.a].x + mesh.nodes[side.b].x);
            const double midpoint_y =
              0.5 * (mesh.nodes[side.a].y + mesh.nodes[side.b].y);
            side.boundary_mark = blasius_boundary_mark(midpoint_x,
                                                       midpoint_y);
            if (side.boundary_mark == 0)
              throw std::runtime_error(
                "unclassified generated Blasius boundary side");
          }
        else
          side.boundary_mark = 0;
        const double length =
          edge_length(mesh.nodes[side.a], mesh.nodes[side.b]);
        mesh.min_edge = std::min(mesh.min_edge, length);
        mesh.max_edge = std::max(mesh.max_edge, length);
      }

    return mesh;
  }

  void write_generated_mesh(const GeneratedMesh &mesh,
                            const std::string &mesh_base)
  {
    std::filesystem::create_directories(
      std::filesystem::path(mesh_base).parent_path());

    {
      std::ofstream out(mesh_base + ".n");
      out << std::setw(6) << mesh.nodes.size() << std::setw(8)
          << mesh.elements.size() << std::setw(8) << mesh.sides.size()
          << " **(Nnd, Nee, Nsd)**\n";
      out << std::scientific << std::setprecision(16);
      for (std::size_t i = 0; i < mesh.nodes.size(); ++i)
        out << std::setw(6) << i << " " << std::setw(24) << mesh.nodes[i].x
            << " " << std::setw(24) << mesh.nodes[i].y << " "
            << std::setw(6) << mesh.nodes[i].boundary_mark << "\n";
    }

    {
      std::ofstream out(mesh_base + ".e");
      out << std::setw(6) << mesh.elements.size() << std::setw(8)
          << mesh.nodes.size() << std::setw(8) << mesh.sides.size()
          << " **(Nee, Nnd, Nsd)**\n";
      for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        {
          const auto &element = mesh.elements[i];
          out << std::setw(6) << i;
          for (const int node : element.node)
            out << std::setw(8) << node;
          for (const int neighbor : element.neighbor)
            out << std::setw(8) << neighbor;
          for (const int side : element.side)
            out << std::setw(8) << side;
          out << "\n";
        }
    }

    {
      std::ofstream out(mesh_base + ".s");
      out << mesh.sides.size() << "\n";
      for (std::size_t i = 0; i < mesh.sides.size(); ++i)
        {
          const auto &side = mesh.sides[i];
          out << std::setw(6) << i << std::setw(8) << side.a
              << std::setw(8) << side.b << std::setw(8) << side.ea
              << std::setw(8) << side.eb << std::setw(8)
              << side.boundary_mark << " \n";
        }
    }
  }

  bool on_plate(const double *p)
  {
    const double half_thickness = active_plate_half_thickness;
    if (half_thickness <= boundary_tolerance)
      return p[0] > -boundary_tolerance &&
             p[0] < 5.0 + boundary_tolerance &&
             near(p[1], 0.0);

    return p[0] > -boundary_tolerance &&
           p[0] < 5.0 + boundary_tolerance &&
           std::abs(p[1]) < half_thickness + boundary_tolerance;
  }

  double velocity_boundary_value(const int component, const double *p)
  {
    if (on_plate(p))
      return 0.0;
    if (component == 0)
      return free_stream_speed;
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

  double zero_boundary_value(const double *)
  {
    return 0.0;
  }

  void apply_homogeneous_blasius_dirichlet_boundary(
    StiffMatrix<dim, double> &matrix,
    FEMFunction<double, dim> &solution,
    Vector<double> &rhs,
    const FEMSpace<double, dim> &space)
  {
    BoundaryFunction<double, dim> inflow(
      BoundaryConditionInfo::DIRICHLET, 1, &zero_boundary_value);
    BoundaryFunction<double, dim> lower_wall(
      BoundaryConditionInfo::DIRICHLET, 2, &zero_boundary_value);
    BoundaryFunction<double, dim> upper_wall(
      BoundaryConditionInfo::DIRICHLET, 4, &zero_boundary_value);
    BoundaryFunction<double, dim> plate(
      BoundaryConditionInfo::DIRICHLET, 6, &zero_boundary_value);
    BoundaryConditionAdmin<double, dim> boundary_admin(space);
    boundary_admin.add(inflow);
    boundary_admin.add(lower_wall);
    boundary_admin.add(upper_wall);
    boundary_admin.add(plate);
    boundary_admin.apply(matrix, solution, rhs);
  }

  bool is_velocity_dirichlet_dof(const FEMSpace<double, dim> &velocity_space,
                                 const int dof)
  {
    if (velocity_space.dofBoundaryMark(dof) == 0)
      return false;

    const auto &point = velocity_space.dofInfo(dof).interp_point;
    const double x = point[0];
    const double y = point[1];

    const double p[dim] = {x, y};
    if (near(x, 5.0) && !near(y, -1.0) && !near(y, 1.0) &&
        !on_plate(p))
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

  void assemble_picard_system(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    const FEMFunction<double, dim> &previous_x,
    const FEMFunction<double, dim> &previous_y,
    const double viscosity,
    const bool newton_linearization,
    RowSparseMatrix &matrix,
    RowSparseMatrix &pressure_mass_matrix,
    std::vector<double> &rhs,
    const int pressure_pin)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    matrix.reinit(n_total);
    pressure_mass_matrix.reinit(n_pressure_dofs);
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
        const auto previous_x_grads =
          previous_x.gradient(q_point, velocity_element);
        const auto previous_y_grads =
          previous_y.gradient(q_point, velocity_element);

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
                      {
                        matrix.add(velocity_index(component,
                                                  vi,
                                                  n_velocity_dofs),
                                   velocity_index(component,
                                                  vj,
                                                  n_velocity_dofs),
                                   velocity_block);
                        if (newton_linearization)
                          {
                            const auto &previous_gradient =
                              component == 0 ? previous_x_grads[q] :
                                               previous_y_grads[q];
                            for (int trial_component = 0;
                                 trial_component < dim;
                                 ++trial_component)
                              matrix.add(
                                velocity_index(component,
                                               vi,
                                               n_velocity_dofs),
                                velocity_index(trial_component,
                                               vj,
                                               n_velocity_dofs),
                                jxw * velocity_values[i][q] *
                                  velocity_values[j][q] *
                                  previous_gradient[trial_component]);
                          }
                      }
                  }

                for (int component = 0; component < dim; ++component)
                  {
                    double effective_rhs = f[component];
                    if (newton_linearization)
                      {
                        const auto &previous_gradient =
                          component == 0 ? previous_x_grads[q] :
                                           previous_y_grads[q];
                        effective_rhs +=
                          wind[0] * previous_gradient[0] +
                          wind[1] * previous_gradient[1];
                      }
                    rhs[velocity_index(component, vi, n_velocity_dofs)] +=
                      jxw * effective_rhs * velocity_values[i][q];
                  }

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj = pressure_element_dofs[j];
                    for (int component = 0; component < dim; ++component)
                      matrix.add(velocity_index(component,
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
                      matrix.add(pressure_index(pi_dof, n_velocity_dofs),
                                 velocity_index(component,
                                                vj,
                                                n_velocity_dofs),
                                 -jxw * pressure_values[i][q] *
                                   velocity_grads[j][q][component]);
                  }

                for (int j = 0; j < n_local_pressure_dofs; ++j)
                  {
                    const int pj_dof = pressure_element_dofs[j];
                    const double mass_entry =
                      jxw * pressure_values[i][q] * pressure_values[j][q];
                    pressure_mass_matrix.add(pi_dof, pj_dof, mass_entry);
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
                                       velocity_index(component,
                                                      i,
                                                      n_velocity_dofs),
                                       velocity_boundary_value(component, p));
        }

    apply_dirichlet_constraint(matrix,
                               rhs,
                               pressure_index(pressure_pin, n_velocity_dofs),
                               0.0);
    apply_homogeneous_constraint(pressure_mass_matrix, pressure_pin);
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

  std::vector<double>
  pack_solution(const FEMFunction<double, dim> &u_x,
                const FEMFunction<double, dim> &u_y,
                const FEMFunction<double, dim> &p_h)
  {
    const int n_velocity_dofs = u_x.femSpace().n_dof();
    const int n_pressure_dofs = p_h.femSpace().n_dof();
    std::vector<double> solution(2 * n_velocity_dofs + n_pressure_dofs, 0.0);
    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        solution[velocity_index(0, i, n_velocity_dofs)] = u_x(i);
        solution[velocity_index(1, i, n_velocity_dofs)] = u_y(i);
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      solution[pressure_index(i, n_velocity_dofs)] = p_h(i);
    return solution;
  }

  double nonlinear_residual_norm(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    const FEMFunction<double, dim> &u_x,
    const FEMFunction<double, dim> &u_y,
    const FEMFunction<double, dim> &p_h,
    const double viscosity,
    const int pressure_pin)
  {
    const int n_total =
      2 * velocity_space.n_dof() + pressure_space.n_dof();
    RowSparseMatrix residual_matrix(n_total);
    RowSparseMatrix pressure_mass_matrix(pressure_space.n_dof());
    std::vector<double> residual_rhs;
    assemble_picard_system(velocity_space,
                           pressure_space,
                           u_x,
                           u_y,
                           viscosity,
                           false,
                           residual_matrix,
                           pressure_mass_matrix,
                           residual_rhs,
                           pressure_pin);

    std::vector<double> residual = residual_matrix.vmult(
      pack_solution(u_x, u_y, p_h));
    add_scaled(residual, -1.0, residual_rhs);
    return norm(residual) / std::max(norm(residual_rhs), 1.0e-30);
  }

  double apply_newton_line_search(
    const FEMSpace<double, dim> &velocity_space,
    const FEMSpace<double, dim> &pressure_space,
    FEMFunction<double, dim> &u_x,
    FEMFunction<double, dim> &u_y,
    FEMFunction<double, dim> &p_h,
    const std::vector<double> &old_solution,
    const std::vector<double> &newton_solution,
    std::vector<double> &accepted_solution,
    const double viscosity,
    const int pressure_pin,
    const double old_residual)
  {
    constexpr double armijo_slope = 1.0e-4;
    constexpr int max_backtracks = 14;

    std::vector<double> best_solution = old_solution;
    double best_residual = old_residual;
    double step = 1.0;
    for (int backtrack = 0; backtrack <= max_backtracks; ++backtrack)
      {
        std::vector<double> candidate(old_solution.size(), 0.0);
        for (std::size_t i = 0; i < candidate.size(); ++i)
          candidate[i] = old_solution[i] +
                         step * (newton_solution[i] - old_solution[i]);
        assign_solution(candidate, u_x, u_y, p_h);
        const double candidate_residual = nonlinear_residual_norm(
          velocity_space,
          pressure_space,
          u_x,
          u_y,
          p_h,
          viscosity,
          pressure_pin);

        if (candidate_residual < best_residual)
          {
            best_residual = candidate_residual;
            best_solution = candidate;
          }
        if (candidate_residual <=
            (1.0 - armijo_slope * step) * old_residual)
          {
            accepted_solution = std::move(candidate);
            return candidate_residual;
          }
        step *= 0.5;
      }

    if (best_residual < old_residual)
      {
        accepted_solution = std::move(best_solution);
        assign_solution(accepted_solution, u_x, u_y, p_h);
        return best_residual;
      }

    assign_solution(old_solution, u_x, u_y, p_h);
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
            << "Newton line search failed for nu=" << viscosity
            << "; nonlinear relative residual=" << old_residual;
    throw std::runtime_error(message.str());
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

  void relax_velocity(FEMFunction<double, dim> &u_x,
                      FEMFunction<double, dim> &u_y,
                      const FEMFunction<double, dim> &old_x,
                      const FEMFunction<double, dim> &old_y,
                      std::vector<double> &linear_solution,
                      const double relaxation)
  {
    const int n_velocity_dofs = u_x.femSpace().n_dof();
    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        u_x(i) = old_x(i) + relaxation * (u_x(i) - old_x(i));
        u_y(i) = old_y(i) + relaxation * (u_y(i) - old_y(i));
        linear_solution[velocity_index(0, i, n_velocity_dofs)] = u_x(i);
        linear_solution[velocity_index(1, i, n_velocity_dofs)] = u_y(i);
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

  double boundary_layer_thickness(const FEMFunction<double, dim> &u_x,
                                  const double x,
                                  const bool upper_half)
  {
    constexpr double threshold = 0.99 * free_stream_speed;
    constexpr double strip_half_width = 0.08;

    std::vector<std::pair<double, double>> samples;
    const auto &space = u_x.femSpace();
    const double plate_offset = active_plate_half_thickness;
    for (std::size_t i = 0; i < u_x.size(); ++i)
      {
        const auto &point = space.dofInfo(static_cast<int>(i)).interp_point;
        if (std::abs(point[0] - x) > strip_half_width)
          continue;
        if (upper_half && point[1] < plate_offset - boundary_tolerance)
          continue;
        if (!upper_half && point[1] > -plate_offset + boundary_tolerance)
          continue;
        const double distance = upper_half ?
                                  point[1] - plate_offset :
                                  -plate_offset - point[1];
        samples.emplace_back(std::max(0.0, distance), u_x(i));
      }

    std::sort(samples.begin(), samples.end());
    if (samples.size() < 2)
      return -1.0;

    if (samples.front().second >= threshold)
      return samples.front().first;

    for (std::size_t i = 1; i < samples.size(); ++i)
      {
        const double previous_distance = samples[i - 1].first;
        const double current_distance = samples[i].first;
        const double previous_value = samples[i - 1].second;
        const double current_value = samples[i].second;

        if (current_value < threshold)
          continue;
        if (current_distance <= previous_distance ||
            std::abs(current_value - previous_value) < 1.0e-14)
          return current_distance;

        const double theta =
          (threshold - previous_value) / (current_value - previous_value);
        return previous_distance + theta * (current_distance -
                                            previous_distance);
      }

    return -1.0;
  }

  double blasius_thickness_scale(const double viscosity, const double x)
  {
    if (viscosity <= 0.0 || x <= 0.0 || free_stream_speed <= 0.0)
      return -1.0;
    return std::sqrt(viscosity * x / free_stream_speed);
  }

  double normalized_layer_thickness(const double thickness,
                                    const double viscosity,
                                    const double x)
  {
    const double scale = blasius_thickness_scale(viscosity, x);
    if (thickness < 0.0 || scale <= 0.0)
      return -1.0;
    return thickness / scale;
  }

  double average_valid_thickness(const double upper, const double lower)
  {
    if (upper >= 0.0 && lower >= 0.0)
      return 0.5 * (upper + lower);
    if (upper >= 0.0)
      return upper;
    if (lower >= 0.0)
      return lower;
    return -1.0;
  }

  std::string viscosity_label(const double viscosity)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(10) << viscosity;
    std::string label = stream.str();
    while (!label.empty() && label.back() == '0')
      label.pop_back();
    if (!label.empty() && label.back() == '.')
      label.pop_back();
    std::replace(label.begin(), label.end(), '.', 'p');
    return label;
  }

  double dof_coordinate_checksum(const FEMSpace<double, dim> &space)
  {
    long double checksum = 0.0;
    for (int i = 0; i < space.n_dof(); ++i)
      {
        const auto &point = space.dofInfo(i).interp_point;
        checksum += static_cast<long double>(i + 1) *
                    (0.125L + static_cast<long double>(point[0]) +
                     7.0L * static_cast<long double>(point[1]));
      }
    return static_cast<double>(checksum);
  }

  void write_blasius_state(const std::string &filename,
                           const double viscosity,
                           const FEMFunction<double, dim> &u_x,
                           const FEMFunction<double, dim> &u_y,
                           const FEMFunction<double, dim> &p_h)
  {
    const std::string temporary_filename = filename + ".tmp";
    std::ofstream output(temporary_filename);
    if (!output)
      throw std::runtime_error("cannot open Blasius state file for writing: " +
                               temporary_filename);

    output << "FEM_BLASIUS_STATE 1\n"
           << u_x.femSpace().n_dof() << ' ' << p_h.femSpace().n_dof() << '\n'
           << std::scientific << std::setprecision(17)
           << dof_coordinate_checksum(u_x.femSpace()) << ' '
           << dof_coordinate_checksum(p_h.femSpace()) << '\n'
           << viscosity << '\n';
    for (std::size_t i = 0; i < u_x.size(); ++i)
      output << u_x(i) << ' ' << u_y(i) << '\n';
    for (std::size_t i = 0; i < p_h.size(); ++i)
      output << p_h(i) << '\n';
    output.close();
    if (!output)
      throw std::runtime_error("failed while writing Blasius state file: " +
                               temporary_filename);

    std::filesystem::rename(temporary_filename, filename);
  }

  double read_blasius_state(const std::string &filename,
                            FEMFunction<double, dim> &old_x,
                            FEMFunction<double, dim> &old_y,
                            FEMFunction<double, dim> &u_x,
                            FEMFunction<double, dim> &u_y,
                            FEMFunction<double, dim> &p_h)
  {
    std::ifstream input(filename);
    if (!input)
      throw std::runtime_error("cannot open Blasius restart state: " + filename);

    std::string magic;
    int version = 0;
    int n_velocity_dofs = 0;
    int n_pressure_dofs = 0;
    double velocity_checksum = 0.0;
    double pressure_checksum = 0.0;
    double viscosity = 0.0;
    input >> magic >> version >> n_velocity_dofs >> n_pressure_dofs >>
      velocity_checksum >> pressure_checksum >> viscosity;
    if (!input || magic != "FEM_BLASIUS_STATE" || version != 1)
      throw std::runtime_error("invalid Blasius restart state header: " +
                               filename);
    if (n_velocity_dofs != old_x.femSpace().n_dof() ||
        n_pressure_dofs != p_h.femSpace().n_dof())
      throw std::runtime_error("Blasius restart state has incompatible dof counts: " +
                               filename);

    const auto checksum_matches = [](const double stored,
                                     const double current) {
      return std::abs(stored - current) <=
             1.0e-12 * std::max({1.0, std::abs(stored), std::abs(current)});
    };
    if (!checksum_matches(velocity_checksum,
                          dof_coordinate_checksum(old_x.femSpace())) ||
        !checksum_matches(pressure_checksum,
                          dof_coordinate_checksum(p_h.femSpace())))
      throw std::runtime_error("Blasius restart state belongs to a different mesh: " +
                               filename);
    if (!(viscosity > 0.0) || !std::isfinite(viscosity))
      throw std::runtime_error("invalid viscosity in Blasius restart state: " +
                               filename);

    for (std::size_t i = 0; i < old_x.size(); ++i)
      {
        double x_value = 0.0;
        double y_value = 0.0;
        input >> x_value >> y_value;
        if (!input || !std::isfinite(x_value) || !std::isfinite(y_value))
          throw std::runtime_error("invalid velocity in Blasius restart state: " +
                                   filename);
        old_x(i) = x_value;
        old_y(i) = y_value;
        u_x(i) = x_value;
        u_y(i) = y_value;
      }
    for (std::size_t i = 0; i < p_h.size(); ++i)
      {
        double value = 0.0;
        input >> value;
        if (!input || !std::isfinite(value))
          throw std::runtime_error("invalid pressure in Blasius restart state: " +
                                   filename);
        p_h(i) = value;
      }

    std::string trailing_data;
    if (input >> trailing_data)
      throw std::runtime_error("unexpected trailing data in Blasius restart state: " +
                               filename);
    return viscosity;
  }

  BlasiusResult solve_viscosity_step(
    const EasyMesh &mesh,
    FEMSpace<double, dim> &velocity_space,
    FEMSpace<double, dim> &pressure_space,
    const double viscosity,
    const std::string &output_prefix,
    FEMFunction<double, dim> &old_x,
    FEMFunction<double, dim> &old_y,
    FEMFunction<double, dim> &u_x,
    FEMFunction<double, dim> &u_y,
    FEMFunction<double, dim> &p_h,
    EigenUmfpackSolverCache &umfpack_solver)
  {
    const int n_velocity_dofs = velocity_space.n_dof();
    const int n_pressure_dofs = pressure_space.n_dof();
    const int n_total = 2 * n_velocity_dofs + n_pressure_dofs;

    RowSparseMatrix matrix(n_total);
    RowSparseMatrix pressure_mass_matrix(n_pressure_dofs);
    std::vector<double> rhs(n_total, 0.0);
    std::vector<double> linear_solution(n_total, 0.0);
    unsigned int last_system_nonzeros = 0;
    int total_gmres_iterations = 0;
    double final_gmres_relative_residual = 1.0;
    double final_relative_update = 1.0;
    double final_nonlinear_residual = 1.0;
    int picard_iteration = 0;

    const int pressure_pin =
      nearest_pressure_dof(pressure_space, 5.0, 0.5);
    const double gmres_tolerance = blasius_gmres_tolerance();
    const int gmres_max_iterations = blasius_gmres_max_iterations();
    const BlasiusLinearSolver linear_solver = blasius_linear_solver();
    const BlasiusNonlinearSolver nonlinear_solver =
      blasius_nonlinear_solver();
    const double picard_relaxation = blasius_picard_relaxation();
    const int max_picard_iterations = blasius_picard_max_iterations();

    for (int i = 0; i < n_velocity_dofs; ++i)
      {
        linear_solution[velocity_index(0, i, n_velocity_dofs)] = old_x(i);
        linear_solution[velocity_index(1, i, n_velocity_dofs)] = old_y(i);
      }
    for (int i = 0; i < n_pressure_dofs; ++i)
      linear_solution[pressure_index(i, n_velocity_dofs)] = p_h(i);

    std::unique_ptr<StiffMatrix<dim, double>> velocity_stiffness;
    std::unique_ptr<AMGPreconditioner> velocity_amg;
    if (linear_solver == BlasiusLinearSolver::schur_gmres)
      {
        velocity_stiffness =
          std::make_unique<StiffMatrix<dim, double>>(velocity_space);
        velocity_stiffness->algebricAccuracy() = 4;
        velocity_stiffness->build();
        FEMFunction<double, dim> velocity_boundary_solution(velocity_space);
        Vector<double> velocity_zero_rhs(n_velocity_dofs);
        velocity_zero_rhs = 0.0;
        apply_homogeneous_blasius_dirichlet_boundary(
          *velocity_stiffness,
          velocity_boundary_solution,
          velocity_zero_rhs,
          velocity_space);
        velocity_amg =
          std::make_unique<AMGPreconditioner>(*velocity_stiffness, 2, 2);
        velocity_amg->smoothStep() = 2;
      }

    for (; picard_iteration < max_picard_iterations; ++picard_iteration)
      {
        assemble_picard_system(velocity_space,
                               pressure_space,
                               old_x,
                               old_y,
                               viscosity,
                               nonlinear_solver ==
                                 BlasiusNonlinearSolver::newton,
                               matrix,
                               pressure_mass_matrix,
                               rhs,
                               pressure_pin);

        if (static_cast<int>(rhs.size()) != n_total)
          throw std::runtime_error("unexpected linear system size");

        last_system_nonzeros = matrix.nonzeros();
        const std::vector<double> old_solution =
          pack_solution(old_x, old_y, p_h);
        const double old_nonlinear_residual =
          nonlinear_solver == BlasiusNonlinearSolver::newton ?
            nonlinear_residual_norm(velocity_space,
                                    pressure_space,
                                    old_x,
                                    old_y,
                                    p_h,
                                    viscosity,
                                    pressure_pin) :
            0.0;
        SolveInfo gmres_info;
        if (linear_solver == BlasiusLinearSolver::eigen_sparse_lu)
          gmres_info = solve_eigen_sparse_lu(matrix, linear_solution, rhs);
        else if (linear_solver == BlasiusLinearSolver::eigen_umfpack)
          gmres_info = umfpack_solver.solve(matrix, linear_solution, rhs);
        else
          {
            SchurMassPreconditioner preconditioner(
              matrix,
              pressure_mass_matrix,
              *velocity_amg,
              viscosity,
              n_velocity_dofs,
              n_pressure_dofs,
              pressure_pin,
              "Navier-Stokes Blasius");
            gmres_info = solve_right_preconditioned_gmres(
              matrix,
              linear_solution,
              rhs,
              preconditioner,
              gmres_tolerance,
              120,
              gmres_max_iterations);
          }
        if (!gmres_info.converged)
          {
            std::ostringstream message;
            message << std::scientific << std::setprecision(6)
                    << blasius_linear_solver_name(linear_solver)
                    << " failed for nu="
                    << viscosity << " in nonlinear step "
                    << picard_iteration + 1 << " after "
                    << gmres_info.iterations
                    << " iterations; final relative residual="
                    << gmres_info.relative_residual;
            throw std::runtime_error(message.str());
          }

        total_gmres_iterations += gmres_info.iterations;
        final_gmres_relative_residual = gmres_info.relative_residual;
        if (nonlinear_solver == BlasiusNonlinearSolver::newton)
          {
            const std::vector<double> newton_solution = linear_solution;
            final_nonlinear_residual = apply_newton_line_search(
              velocity_space,
              pressure_space,
              u_x,
              u_y,
              p_h,
              old_solution,
              newton_solution,
              linear_solution,
              viscosity,
              pressure_pin,
              old_nonlinear_residual);
          }
        else
          assign_solution(linear_solution, u_x, u_y, p_h);
        if (nonlinear_solver == BlasiusNonlinearSolver::picard &&
            picard_relaxation < 1.0)
          relax_velocity(u_x,
                         u_y,
                         old_x,
                         old_y,
                         linear_solution,
                         picard_relaxation);

        final_relative_update =
          relative_velocity_update(old_x, old_y, u_x, u_y);
        ++picard_iteration;

        copy_velocity(old_x, old_y, u_x, u_y);

        if (final_relative_update < picard_tolerance &&
            (nonlinear_solver == BlasiusNonlinearSolver::picard ||
             final_nonlinear_residual < picard_tolerance))
          break;
      }

    final_nonlinear_residual = nonlinear_residual_norm(velocity_space,
                                                       pressure_space,
                                                       u_x,
                                                       u_y,
                                                       p_h,
                                                       viscosity,
                                                       pressure_pin);
    const bool nonlinear_converged =
      final_relative_update < picard_tolerance &&
      (nonlinear_solver == BlasiusNonlinearSolver::picard ||
       final_nonlinear_residual < picard_tolerance);
    if (!nonlinear_converged)
      {
        std::ostringstream message;
        message << std::scientific << std::setprecision(6)
                << blasius_nonlinear_solver_name(nonlinear_solver)
                << " iteration did not converge for nu=" << viscosity
                << "; final relative update=" << final_relative_update
                << "; final nonlinear relative residual="
                << final_nonlinear_residual;
        throw std::runtime_error(message.str());
      }

    const std::string step_prefix =
      output_prefix + "_nu" + viscosity_label(viscosity);
    u_x.writeOpenDXData(step_prefix + "_ux.dx");
    u_y.writeOpenDXData(step_prefix + "_uy.dx");
    p_h.writeOpenDXData(step_prefix + "_p.dx");
    write_blasius_state(step_prefix + ".state", viscosity, u_x, u_y, p_h);

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

    const double inflow_rate = vertical_line_flow_rate(u_x, -1.0, -1.0, 1.0);
    const double outflow_rate = vertical_line_flow_rate(u_x, 5.0, -1.0, 1.0);

    return {viscosity,
            plate_length * free_stream_speed / viscosity,
            static_cast<int>(mesh.n_geometry(dim)),
            velocity_space.n_dof(),
            pressure_space.n_dof(),
            last_system_nonzeros,
            picard_iteration,
            total_gmres_iterations,
            final_relative_update,
            final_gmres_relative_residual,
            final_nonlinear_residual,
            std::sqrt(ux_l2 * ux_l2 + uy_l2 * uy_l2),
            std::sqrt(ux_h1 * ux_h1 + uy_h1 * uy_h1),
            max_speed(u_x, u_y),
            inflow_rate,
            outflow_rate,
            outflow_rate - inflow_rate,
            pressure_min,
            pressure_max,
            pressure_integral(p_h) / active_blasius_domain_area,
            divergence_l2_norm(u_x, u_y),
            boundary_layer_thickness(u_x, 1.0, true),
            boundary_layer_thickness(u_x, 1.0, false),
            boundary_layer_thickness(u_x, 4.0, true),
            boundary_layer_thickness(u_x, 4.0, false)};
  }

  std::vector<BlasiusResult> solve_blasius_continuation(
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

    if (const char *restart_state = std::getenv("BLASIUS_RESTART_STATE"))
      {
        const double restart_viscosity = read_blasius_state(restart_state,
                                                            old_x,
                                                            old_y,
                                                            u_x,
                                                            u_y,
                                                            p_h);
        std::cout << std::scientific << std::setprecision(6)
                  << "loaded Blasius restart state " << restart_state
                  << " at nu=" << restart_viscosity << '\n';
      }

    std::vector<BlasiusResult> results;
    results.reserve(viscosities.size());
    EigenUmfpackSolverCache umfpack_solver;

    for (const double viscosity : viscosities)
      results.push_back(solve_viscosity_step(mesh,
                                             velocity_space,
                                             pressure_space,
                                             viscosity,
                                             output_prefix,
                                             old_x,
                                             old_y,
                                             u_x,
                                             u_y,
                                             p_h,
                                             umfpack_solver));

    return results;
  }

  void print_thickness_summary(const std::vector<BlasiusResult> &results,
                               const std::string &title)
  {
    std::cout << "\n" << title << "\n";
    std::cout << "nu Re Nonlinear delta_x1_avg delta_x4_avg "
                 "delta_x1/sqrt(nu) delta_x4/sqrt(4nu)\n";
    for (const auto &result : results)
      {
        const double delta_x1 =
          average_valid_thickness(result.upper_layer_x1,
                                  result.lower_layer_x1);
        const double delta_x4 =
          average_valid_thickness(result.upper_layer_x4,
                                  result.lower_layer_x4);
        std::cout << result.viscosity << ' '
                  << result.reynolds_number << ' '
                  << result.nonlinear_iterations << ' '
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

  std::vector<double> parse_viscosities(const std::string &argument)
  {
    if (argument.empty())
      return {1.0, 0.2, 0.1, 0.05, 0.02, 0.01};

    std::vector<double> viscosities;
    std::stringstream   stream(argument);
    std::string         token;
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
               std::string(AFEPACK_DEFAULT_MESH_DIR) + "/blasius_plate_h0p25";
  const std::string output_prefix =
    argc > 2 ? argv[2] :
               std::string(AFEPACK_DEFAULT_OUTPUT_DIR) +
                 "/navier_stokes_blasius_afepack";
  const auto viscosities =
    parse_viscosities(argc > 3 ? argv[3] : std::string());

  try
    {
      const auto results =
        solve_blasius_continuation(mesh_file, output_prefix, viscosities);

      std::cout << std::scientific << std::setprecision(6);
      std::cout << "AFEPack Navier-Stokes Blasius thin-plate proxy "
                   "continuation\n";
      for (const auto &result : results)
        {
          std::cout << "nu: " << result.viscosity
                    << " book Re=5/nu: " << result.reynolds_number << '\n';
          std::cout << "linear solver: "
                    << blasius_linear_solver_name(blasius_linear_solver())
                    << '\n';
          std::cout << "nonlinear solver: "
                    << blasius_nonlinear_solver_name(
                         blasius_nonlinear_solver())
                    << '\n';
          std::cout << "cells: " << result.cells << '\n';
          std::cout << "velocity dofs: " << result.velocity_dofs << '\n';
          std::cout << "pressure dofs: " << result.pressure_dofs << '\n';
          std::cout << "system nonzeros: " << result.system_nonzeros << '\n';
          std::cout << "nonlinear iterations: "
                    << result.nonlinear_iterations << '\n';
          std::cout << "total Schur-preconditioned GMRES iterations: "
                    << result.total_gmres_iterations << '\n';
          std::cout << "final GMRES relative residual: "
                    << result.final_gmres_relative_residual << '\n';
          std::cout << "final nonlinear relative residual: "
                    << result.final_nonlinear_relative_residual << '\n';
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
          std::cout << "ux=0.99 layer thickness at x=1 upper/lower: "
                    << result.upper_layer_x1 << " "
                    << result.lower_layer_x1 << '\n';
          std::cout << "delta/sqrt(nu*x/U) at x=1 upper/lower: "
                    << normalized_layer_thickness(result.upper_layer_x1,
                                                  result.viscosity,
                                                  1.0)
                    << " "
                    << normalized_layer_thickness(result.lower_layer_x1,
                                                  result.viscosity,
                                                  1.0)
                    << '\n';
          std::cout << "ux=0.99 layer thickness at x=4 upper/lower: "
                    << result.upper_layer_x4 << " "
                    << result.lower_layer_x4 << '\n';
          std::cout << "delta/sqrt(nu*x/U) at x=4 upper/lower: "
                    << normalized_layer_thickness(result.upper_layer_x4,
                                                  result.viscosity,
                                                  4.0)
                    << " "
                    << normalized_layer_thickness(result.lower_layer_x4,
                                                  result.viscosity,
                                                  4.0)
                    << '\n';
        }

      print_thickness_summary(results, "Blasius thickness-law proxy summary");

      if (argc == 1)
        {
          const std::string stretched_mesh_file =
            std::string(AFEPACK_DEFAULT_MESH_DIR) +
            "/generated_blasius_plate_stretched";
          const GeneratedMesh stretched_mesh = generate_stretched_blasius_mesh();
          write_generated_mesh(stretched_mesh, stretched_mesh_file);

          const std::vector<double> stretched_viscosities = {0.1, 0.02};
          const auto stretched_results =
            solve_blasius_continuation(stretched_mesh_file,
                                       output_prefix + "_stretched",
                                       stretched_viscosities);

          std::cout << "\nGenerated stretched Blasius proxy mesh\n";
          std::cout << "cells: " << stretched_mesh.elements.size()
                    << " nodes: " << stretched_mesh.nodes.size()
                    << " min_edge: " << stretched_mesh.min_edge
                    << " max_edge: " << stretched_mesh.max_edge << '\n';
          print_thickness_summary(
            stretched_results,
            "Blasius stretched-grid thickness-law proxy summary");
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
