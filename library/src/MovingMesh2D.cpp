/**
 * @file   MovingMesh2D.cpp
 * @author Robert Lie
 * @date   Thu Feb 22 11:51:37 2007
 * 
 * @brief  
 * 
 * 
 */

#include "MovingMesh2D.h"

AFEPACK_OPEN_NAMESPACE

#define DIM 2

typedef SparseMatrix<double> sparse_matrix_t;
typedef Vector<double> vector_t;

/// 素数表
int MovingMesh2D::primes[] = {
  02, 03, 05, 07, 11, 
  13, 17, 19, 23, 29, 
  31, 37, 41, 43, 47,
  53, 59, 61, 67, 71,
  73, 79, 83, 89, 91
};

MovingMesh2D::MovingMesh2D() :
  tol(1.0e-02), solve_step(5), 
  max_step(100), n_move_step(1)
{}



MovingMesh2D::~MovingMesh2D()
{}



void MovingMesh2D::readDomain(const std::string& filename)
{
  int i, j, k, l, l0, l1;
  double h;

  readData(filename);
  logical_node.resize(n_geometry(0));
  move_direction.resize(n_geometry(0));
  logical_move_direction.resize(n_geometry(0));
  monitor().resize(n_geometry(DIM));

  std::ifstream is;
  is.open((filename + ".d").c_str());
  if (!is) {
    std::cerr << "Open the easymesh input file "
              << filename + ".d failure, aborting ... " 
              << std::endl;
    abort();
  }
  readDummy(is); is >> domain.n_vertex;
  domain.physical_domain_vertex.resize(domain.n_vertex);
  domain.logical_domain_vertex.resize(domain.n_vertex);
  for (i = 0;i < domain.n_vertex;i ++) {
    readDummy(is); is >> domain.physical_domain_vertex[i].index;
    readDummy(is); is >> domain.physical_domain_vertex[i][0];
    readDummy(is); is >> domain.physical_domain_vertex[i][1];
    readDummy(is); is >> h;
    readDummy(is); is >> domain.physical_domain_vertex[i].boundary_mark;
    domain.physical_domain_vertex[i].boundary_mark = 1; /// 缺省值设为 1
  }
  readDummy(is); is >> domain.n_edge;
  domain.edge.resize(domain.n_edge);
  for (i = 0;i < domain.n_edge;i ++) {
    readDummy(is); is >> domain.edge[i].index;
    readDummy(is); is >> domain.edge[i].vertex[0];
    readDummy(is); is >> domain.edge[i].vertex[1];
    readDummy(is); is >> domain.edge[i].boundary_mark;
  }
  is.close();
	
  is.open((filename + ".log").c_str());
  if (!is) {
    std::cerr << "Open the logical domain description file "
              << filename + ".log failure.\n" 
              << "+-----------------------------------------------------+\n"
              << "| Warning: The vertex of the physical domain is used. |\n"
              << "+-----------------------------------------------------+\n"
              << std::endl;
    for (i = 0;i < domain.n_vertex;i ++) {
      domain.logical_domain_vertex[i][0] = domain.physical_domain_vertex[i][0];
      domain.logical_domain_vertex[i][1] = domain.physical_domain_vertex[i][1];
    }
  }
  else {
    for (i = 0;i < domain.n_vertex;i ++) {
      is >> domain.logical_domain_vertex[i][0];
      is >> domain.logical_domain_vertex[i][1];
    }
    is.close();
  }

  parseBoundary();

  /// 准备稀疏矩阵的稀疏模板
  u_int n_node = n_geometry(0);
  std::vector<u_int> n_coupling_node(n_node, 1);
  for (i = 0;i < n_geometry(1);i ++) {
    n_coupling_node[geometry(1,i).vertex(0)] ++;
    n_coupling_node[geometry(1,i).vertex(1)] ++;
  }
  spM.reinit(n_node, n_node, n_coupling_node);
  for (i = 0;i < n_geometry(1);i ++) { /**< 对每个边做循环 */
    l0 = geometry(1,i).vertex(0);
    l1 = geometry(1,i).vertex(1);
    spM.add(l0, l0); spM.add(l0, l1);
    spM.add(l1, l0); spM.add(l1, l1);
  }
  spM.compress();

  getLogicalMesh();
}



void MovingMesh2D::moveMesh()
{
  int i, j, k, m, step = 0;
  while (1) {
    getMoveDirection();

    double error = 0;
    for (i = 0;i < n_geometry(2);i ++) {
      const Point<2>& x0 = point(geometry(2,i).vertex(0));
      const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
      const Point<2>& x2 = point(geometry(2,i).vertex(2));
      double d[3], l[3] = {
        (x2[0] - x1[0])*(x2[0] - x1[0]) + (x2[1] - x1[1])*(x2[1] - x1[1]),
        (x0[0] - x2[0])*(x0[0] - x2[0]) + (x0[1] - x2[1])*(x0[1] - x2[1]),
        (x1[0] - x0[0])*(x1[0] - x0[0]) + (x1[1] - x0[1])*(x1[1] - x0[1])
      };
      double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
      double h = 0.5*area/sqrt(*std::max_element(&l[0], &l[3]));
      for (k = 0;k < 3;++ k) {
	m = geometry(2,i).vertex(k);
	for (d[k] = 0.0, j = 0;j < 2;j ++) {
	  d[k] += move_direction[m][j]*move_direction[m][j];
	}
      }
      h = sqrt(*std::max_element(&d[0], &d[3]))/h;
      if (error < h) error = h;
    }

    getMoveStepLength();
    for (i = 0;i < n_move_step;i ++) {
      updateSolution();
      updateMesh();
    }
    std::cerr << "step " << ++ step
              << ": energy reduction = " << energy_reduction
              << ", mesh moving error = " << error 
              << "." << std::endl;
    if (error < tolerence()) break;
    if (maxStep() > 0 && step > maxStep()) break; /// 超过预设步数最大值，则停止
  };
}



void MovingMesh2D::outputPhysicalMesh(const std::string& filename)
{
  writeData(filename);
}



void MovingMesh2D::outputLogicalMesh(const std::string& filename)
{
  Mesh<DIM,DIM> mesh = Mesh<DIM,DIM>(*this);
  for (int i = 0;i < n_geometry(0);i ++) {
    mesh.point(i) = logical_node[i];
  }
  mesh.writeData(filename);
}



void MovingMesh2D::getMoveDirection()
{
  int i, j, k, l, l0, l1;
  getMonitor();
  
  /// 取出控制函数中的最小值，以便做正则化
  double min_monitor = *std::min_element(monitor().begin(), monitor().end());
  M.reinit(spM);
  for (i = 0;i < n_geometry(DIM);i ++) {
    const Point<DIM>& x0 = point(geometry(DIM,i).vertex(0));
    const Point<DIM>& x1 = point(geometry(DIM,i).vertex(1)); 
    const Point<DIM>& x2 = point(geometry(DIM,i).vertex(2));
    double omega[2][3] = {
      {x2[0] - x1[0], x0[0] - x2[0], x1[0] - x0[0]},
      {x2[1] - x1[1], x0[1] - x2[1], x1[1] - x0[1]}
    };
    double area = ((x1[0] - x0[0])*(x2[1] - x0[1]) - 
                   (x2[0] - x0[0])*(x1[1] - x0[1]));
    double G = (monitor(i)/min_monitor)/area;
    for (j = 0;j < 3;j ++) {
      l0 = geometry(DIM,i).vertex(j);
      for (k = 0;k < 3;k ++) {
        l1 = geometry(DIM,i).vertex(k);
        double d = G*(omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k]);
        M.add(l0, l1, d);
      }
    }
  }
  std::vector<vector_t> r(DIM, vector_t(n_geometry(0)));
  std::vector<vector_t> x(DIM, vector_t(n_geometry(0)));
  for (j = 0;j < n_geometry(0);j ++) {
    for (k = 0;k < DIM;++ k) {
      x[k](j) = logical_node[j][k];
    }
  }
  solver.reinit(M, boundary_mark, domain);
  solver.solve(x, r, solveStep());

  for (j = 0;j < n_geometry(0);j ++) {
    for (k = 0;k < DIM;++ k) {
      logical_move_direction[j][k] = logical_node[j][k] - x[k](j);
      move_direction[j][k] = 0;
    }
  }

  /// 计算能量的下降幅度，用来判断是否可以退出
  double energy_before, energy_after;
  for (energy_after = 0.0, k = 0;k < DIM;++ k) {
    M.vmult(r[k], x[k]);
    energy_after += innerProduct(r[k], x[k]);
  }
  for (j = 0;j < n_geometry(0);j ++) {
    for (k = 0;k < DIM;++ k) {
      x[k](j) = logical_node[j][k];
    }
  }
  for (energy_before = 0.0, k = 0;k < DIM;++ k) {
    M.vmult(r[k], x[k]);
    energy_before += innerProduct(r[k], x[k]);
  }
  energy_reduction = 2*fabs(energy_before - energy_after)/
                           (energy_before + energy_after);

  std::vector<double> mass_lumping(n_geometry(0), 0.);
  for (i = 0;i < n_geometry(DIM);i ++) {
    const int& v0 = geometry(DIM,i).vertex(0);
    const int& v1 = geometry(DIM,i).vertex(1);
    const int& v2 = geometry(DIM,i).vertex(2);
    const Point<DIM>& x0 = point(v0); const Point<DIM>& xi0 = logical_node[v0];
    const Point<DIM>& x1 = point(v1); const Point<DIM>& xi1 = logical_node[v1];
    const Point<DIM>& x2 = point(v2); const Point<DIM>& xi2 = logical_node[v2];
    double jacobi[2][2];
    double area = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1]) - 
                   (xi1[1] - xi0[1])*(xi2[0] - xi0[0]));
    jacobi[0][0] = ((x1[0] - x0[0])*(xi2[1] - xi0[1]) -
                    (x2[0] - x0[0])*(xi1[1] - xi0[1]));
    jacobi[1][0] = ((x1[1] - x0[1])*(xi2[1] - xi0[1]) -
                    (x2[1] - x0[1])*(xi1[1] - xi0[1]));
    jacobi[0][1] = ((x2[0] - x0[0])*(xi1[0] - xi0[0]) -
                    (x1[0] - x0[0])*(xi2[0] - xi0[0]));
    jacobi[1][1] = ((x2[1] - x0[1])*(xi1[0] - xi0[0]) -
                    (x1[1] - x0[1])*(xi2[0] - xi0[0]));
    for (j = 0;j < 3;j ++) {
      k = geometry(DIM,i).vertex(j);
      move_direction[k][0] += (jacobi[0][0]*logical_move_direction[k][0] + 
                               jacobi[0][1]*logical_move_direction[k][1]);
      move_direction[k][1] += (jacobi[1][0]*logical_move_direction[k][0] + 
                               jacobi[1][1]*logical_move_direction[k][1]);
      mass_lumping[k] += area;
    }
  }
  for (i = 0;i < n_geometry(0);i ++) {
    move_direction[i][0] /= mass_lumping[i];
    move_direction[i][1] /= mass_lumping[i];
    if (boundary_mark[i] == 1) continue;
    int n_constraint = 0;
    for (j = 0;j < domain.n_edge;++ j) {
      if (boundary_mark[i]%domain.edge[j].boundary_mark != 0) continue;
      double d = (move_direction[i][0]*domain.edge[j].normal[0] +
                  move_direction[i][1]*domain.edge[j].normal[1]);
      move_direction[i][0] -= d*domain.edge[j].normal[0];
      move_direction[i][1] -= d*domain.edge[j].normal[1];
      n_constraint += 1;
    }
    if (n_constraint > 1) {
      move_direction[i][0] = 0.0;
      move_direction[i][1] = 0.0;
    }
  }
}



void MovingMesh2D::getMoveStepLength()
{
  int i, j;
  double a, b, c;
  move_step_length = 1.0;
  for (i = 0;i < n_geometry(DIM);i ++) {
    const int& v0 = geometry(DIM,i).vertex(0);
    const int& v1 = geometry(DIM,i).vertex(1);
    const int& v2 = geometry(DIM,i).vertex(2);
    const Point<DIM>& x0 = point(v0);
    const Point<DIM>& x1 = point(v1); 
    const Point<DIM>& x2 = point(v2);
    a = ((move_direction[v1][0] - move_direction[v0][0])*
         (move_direction[v2][1] - move_direction[v0][1]) - 
         (move_direction[v1][1] - move_direction[v0][1])*
         (move_direction[v2][0] - move_direction[v0][0]));
    b = move_direction[v0][0]*(x1[1] - x2[1])
      - move_direction[v0][1]*(x1[0] - x2[0])
      + move_direction[v1][0]*(x2[1] - x0[1])
      - move_direction[v1][1]*(x2[0] - x0[0])
      + move_direction[v2][0]*(x0[1] - x1[1])
      - move_direction[v2][1]*(x0[0] - x1[0]);
    c = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
    if (fabs(a)/(fabs(b) + fabs(c)) < 1.0e-04) {
      if (fabs(b) < 1.0e-04*fabs(c))
        ; //do nothing
      else if (c/b > 0)
        ; //do nothing
      else
        move_step_length = std::min(move_step_length, -c/b);
    }
    else if (b*b - 4*a*c < 0)
      ; //do nothing
    else {
      if (a < 0) {
        a = -a;
        b = -b;
        c = -c;
      }
      double d = (-b - sqrt(b*b - 4*a*c))/(2*a);
      if (d < 0) {
        d = (-b + sqrt(b*b - 4*a*c))/(2*a);
        if (d > 0)
          move_step_length = std::min(move_step_length, d);
      }
      else
        move_step_length = std::min(move_step_length, d);
    }
  }
  move_step_length *= 0.5;
}



void MovingMesh2D::readDummy(std::ifstream& is)
{
  char c;
  for (is.get(c);true;is.get(c)) {
    if (c == '#') {
      while (true) {
        is.get(c);
        if (c == '#') break;
      }
    }
    else if (('0' <= c && c <= '9')||(c == '.')||(c =='-'))
      break;
  }
  is.putback(c);		
}



void MovingMesh2D::getMonitor()
{
  std::fill(monitor().begin(), monitor().end(), 1);
}



void MovingMesh2D::smoothMonitor(int s)
{
  int i, j, k;
  std::vector<double> area(n_geometry(DIM));
  std::vector<double> mass_lumping(n_geometry(0), 0);
  std::vector<double> monitor1(n_geometry(0));
  for (i = 0;i < n_geometry(DIM);i ++) {
    const Point<DIM>& x0 = point(geometry(DIM,i).vertex(0));
    const Point<DIM>& x1 = point(geometry(DIM,i).vertex(1)); 
    const Point<DIM>& x2 = point(geometry(DIM,i).vertex(2));
    area[i] = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
    for (j = 0;j < 3;j ++)
      mass_lumping[geometry(DIM,i).vertex(j)] += area[i];
  }
  for (i = 0;i < s;i ++) {
    std::fill(monitor1.begin(), monitor1.end(), 0);
    for (j = 0;j < n_geometry(DIM);j ++) {
      for (k = 0;k < 3;k ++) {
        monitor1[geometry(DIM,j).vertex(k)] += monitor(j)*area[j];
      }
    }
    for (j = 0;j < n_geometry(0);j ++)
      monitor1[j] /= 3*mass_lumping[j];
    std::fill(monitor().begin(), monitor().end(), 0);
    for (j = 0;j < n_geometry(DIM);j ++) {
      for (k = 0;k < 3;k ++) {
        monitor(j) += monitor1[geometry(DIM,j).vertex(k)];
      }
    }
  }
}



void MovingMesh2D::updateMesh()
{
  for (int i = 0;i < n_geometry(0);i ++) {
    point(i)[0] += move_step_length * move_direction[i][0];
    point(i)[1] += move_step_length * move_direction[i][1];
  }
}



void MovingMesh2D::getLogicalMesh()
{
  std::cout << "Computing logical mesh ..." << std::endl;
  int i, j, k, l, l0, l1;
  /**
   * 计算逻辑网格的步骤分为两步：
   *
   *   1. 我们在边界上给定线性映射，解狄氏边值问题的调和映照，作为一个
   *      初始值；
   *   2. 以刚才的初始值，解边界控制问题，得到更合适的调和映照；
   * 
   */

  /**
   * 现在先做第一步
   * 
   */
  M.reinit(spM);
  for (i = 0;i < n_geometry(DIM);i ++) {
    const Point<DIM>& x0 = point(geometry(DIM,i).vertex(0));
    const Point<DIM>& x1 = point(geometry(DIM,i).vertex(1)); 
    const Point<DIM>& x2 = point(geometry(DIM,i).vertex(2));
    double omega[2][3] = {
      {x2[0] - x1[0], x0[0] - x2[0], x1[0] - x0[0]},
      {x2[1] - x1[1], x0[1] - x2[1], x1[1] - x0[1]}
    };
    double area = ((x1[0] - x0[0])*(x2[1] - x0[1]) - 
                   (x2[0] - x0[0])*(x1[1] - x0[1]));
    for (j = 0;j < 3;j ++) {
      l0 = geometry(DIM,i).vertex(j);
      for (k = 0;k < 3;k ++) {
        l1 = geometry(DIM,i).vertex(k);
        double d = (omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k])/area;
        M.add(l0, l1, d);
      }
    }
  }
  std::vector<vector_t> r(DIM, vector_t(n_geometry(0)));
  std::vector<vector_t> x(DIM, vector_t(n_geometry(0)));
  const std::size_t * row_start = spM.get_rowstart_indices();
  const u_int * column_index = spM.get_column_numbers();
  for (j = 0;j < n_geometry(0);j ++) {
    if (boundary_mark[j] == 1) continue; /// 跳过内部节点
    const GeometryBM& node = geometry(0,j);
    const Point<DIM>& pnt = point(node.vertex(0));
    u_int edge_idx = 0;
    for (;edge_idx < domain.n_edge;++ edge_idx) {
      if (boundary_mark[j]%domain.edge[edge_idx].boundary_mark == 0) {
        break; /// 检索此点位于哪条边上
      }
    }
    const Edge& edge = domain.edge[edge_idx]; /// 就是这条边

    /// 计算物理区域上，此节点在边上所在的位置的比例
    const Vertex& vtx0 = domain.physical_domain_vertex[edge.vertex[0]];
    const Vertex& vtx1 = domain.physical_domain_vertex[edge.vertex[1]];
    double lambda = (pnt - vtx0).length()/(vtx1 - vtx0).length();

    /// 计算逻辑区域上，相应的逻辑节点的坐标，并对线性系统进行修正
    const Vertex& lvtx0 = domain.logical_domain_vertex[edge.vertex[0]];
    const Vertex& lvtx1 = domain.logical_domain_vertex[edge.vertex[1]];
    for (k = 0;k < DIM;++ k) {
      x[k](j) = (1 - lambda)*lvtx0[k] + lambda*lvtx1[k];
      r[k](j) = M.diag_element(j)*x[k](j);
    }
    for (k = row_start[j] + 1;k < row_start[j + 1];++ k) {
      M.global_entry(k) = 0.0;
      l = column_index[k];
      const u_int * p = std::find(&column_index[row_start[l] + 1],
                                  &column_index[row_start[l + 1]], j);
      if (p != &column_index[row_start[l + 1]]) {
        l0 = p - &column_index[row_start[0]];
        for (l1 = 0;l1 < DIM;++ l1) {
          r[l1](l) -= M.global_entry(l0)*x[l1](j);
        }
        M.global_entry(l0) = 0.;
      }
    }
  }
  AMGSolver amg_solver;
  amg_solver.lazyReinit(M);
  for (k = 0;k < DIM;++ k) {
    amg_solver.solve(x[k], r[k], 1.0e-08, 50);
  }

  /**
   * 现在做第二步
   * 
   */
  /// Todo: ...
#if 0
  M.reinit(spM);
  for (i = 0;i < n_geometry(DIM);i ++) {
    const Point<DIM>& x0 = point(geometry(DIM,i).vertex(0));
    const Point<DIM>& x1 = point(geometry(DIM,i).vertex(1)); 
    const Point<DIM>& x2 = point(geometry(DIM,i).vertex(2));
    double omega[2][3] = {
      {x2[0] - x1[0], x0[0] - x2[0], x1[0] - x0[0]},
      {x2[1] - x1[1], x0[1] - x2[1], x1[1] - x0[1]}
    };
    double area = ((x1[0] - x0[0])*(x2[1] - x0[1]) - 
                   (x2[0] - x0[0])*(x1[1] - x0[1]));
    for (j = 0;j < 3;j ++) {
      l0 = geometry(DIM,i).vertex(j);
      for (k = 0;k < 3;k ++) {
        l1 = geometry(DIM,i).vertex(k);
        double d = (omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k])/area;
        M.add(l0, l1, d);
      }
    }
  }
  for (k = 0;k < DIM;++ k) r[k] = 0.0;
  solver.reinit(M, boundary_mark, domain);
  solver.solve(x, r, solveStep());
#endif 

  for (j = 0;j < n_geometry(0);j ++) {
    for (k = 0;k < DIM;++ k) {
      logical_node[j][k] = x[k](j);
    }
  }
}



void MovingMesh2D::parseBoundary()
{
  std::cout << "Parsing boundary nodes ..." << std::endl;

  u_int n_node = n_geometry(0);
  u_int n_edge = n_geometry(1);

  /// 将第 j 条边上的材料标识设置为第 j 个素数
  std::vector<int> edge_idx(n_edge); /// 网格的每条边位于区域的边的序号
  for (u_int i = 0;i < n_edge;++ i) { /// 对所有的边做循环
    const GeometryBM& edge = geometry(1, i);
    const bound_t& bm = edge.boundaryMark();
    if (bm == 0) { /// 寻找位于边界上的边，而这条边在内部
      edge_idx[i] = -1;
      continue; 
    }
    for (u_int j = 0;j < domain.n_edge;++ j) { /// 否则搜索
      if (bm == domain.edge[j].boundary_mark) {
        edge_idx[i] = j;
        break;
      }
    }
  }
    
  /// 对 domain 对象中的信息进行更新和进一步的处理
  for (u_int i = 0;i < domain.n_edge;++ i) {
    /// 将边界的材料标识设为第 i 个素数
    domain.edge[i].boundary_mark = MovingMesh2D::primes[i]; 

    /// 设置每个边上的单位法向向量
    const u_int& v0 = domain.edge[i].vertex[0];
    const u_int& v1 = domain.edge[i].vertex[1];

    /// 物理区域上的法向量
    const Vertex& p0 = domain.physical_domain_vertex[v0];
    const Vertex& p1 = domain.physical_domain_vertex[v1];
    double dp = (p0 - p1).length();
    domain.edge[i].normal[0] = (p1[1] - p0[1])/dp;
    domain.edge[i].normal[1] = (p0[0] - p1[0])/dp;

    /// 逻辑区域上的法向量
    const Vertex& lp0 = domain.logical_domain_vertex[v0];
    const Vertex& lp1 = domain.logical_domain_vertex[v1];
    double dlp = (lp0 - lp1).length();
    domain.edge[i].logical_normal[0] = (lp1[1] - lp0[1])/dlp;
    domain.edge[i].logical_normal[1] = (lp0[0] - lp1[0])/dlp;

    /// 将顶点的材料标识设为其所属于的两个边的材料标识的乘积
    domain.physical_domain_vertex[v0].boundary_mark *= domain.edge[i].boundary_mark;
    domain.physical_domain_vertex[v1].boundary_mark *= domain.edge[i].boundary_mark;
  }

  /// 统计每条区域边界上的节点的个数，并为每个节点设置正确的材料标识
  boundary_mark.resize(n_node, 1);
  for (u_int i = 0;i < n_edge;++ i) { /// 对所有的边做循环
    const GeometryBM& edge = geometry(1,i);
    const int& idx = edge_idx[i];
    if (idx == -1) continue; /// 这条边在内部，我们跳过去
    const int& bm = domain.edge[idx].boundary_mark;
    for (u_int j = 0;j < 2;++ j) { /// 考虑其两个端点
      const u_int& v = edge.vertex(j);
      if (boundary_mark[v]%bm != 0) {
        boundary_mark[v] *= bm; /// 那么就先设置为本条边的材料标识
      }
    }
  }
}



std::vector<double> 
MovingMesh2D::moveDirection(const Point<DIM>& p, 
                            const int& n) const
{
  const int& v0 = geometry(DIM,n).vertex(0);
  const int& v1 = geometry(DIM,n).vertex(1);
  const int& v2 = geometry(DIM,n).vertex(2);
  const Point<DIM>& x0 = point(v0);
  const Point<DIM>& x1 = point(v1); 
  const Point<DIM>& x2 = point(v2);
  const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
  double lambda[3];
  lambda[0] = ((x1[0] - p[0])*(x2[1] - p[1]) - 
               (x1[1] - p[1])*(x2[0] - p[0]))/area;
  lambda[1] = ((x2[0] - p[0])*(x0[1] - p[1]) - 
               (x2[1] - p[1])*(x0[0] - p[0]))/area;
  lambda[2] = ((x0[0] - p[0])*(x1[1] - p[1]) - 
               (x0[1] - p[1])*(x1[0] - p[0]))/area;
  std::vector<double> v(2);
  v[0] = lambda[0]*m00 + lambda[1]*m10 + lambda[2]*m20;
  v[1] = lambda[0]*m01 + lambda[1]*m11 + lambda[2]*m21;
  return v;
}




std::vector<std::vector<double> > 
MovingMesh2D::moveDirection(const std::vector<Point<DIM> >& p, 
                            const int& n) const
{
  const int& v0 = geometry(DIM,n).vertex(0);
  const int& v1 = geometry(DIM,n).vertex(1);
  const int& v2 = geometry(DIM,n).vertex(2);
  const Point<DIM>& x0 = point(v0);
  const Point<DIM>& x1 = point(v1); 
  const Point<DIM>& x2 = point(v2);
  const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
  std::vector<std::vector<double> > v(p.size(), std::vector<double>(2));
  double lambda[3];
  for (int j = 0;j < p.size();j ++) {
    lambda[0] = ((x1[0] - p[j][0])*(x2[1] - p[j][1]) - 
                 (x1[1] - p[j][1])*(x2[0] - p[j][0]))/area;
    lambda[1] = ((x2[0] - p[j][0])*(x0[1] - p[j][1]) - 
                 (x2[1] - p[j][1])*(x0[0] - p[j][0]))/area;
    lambda[2] = ((x0[0] - p[j][0])*(x1[1] - p[j][1]) - 
                 (x0[1] - p[j][1])*(x1[0] - p[j][0]))/area;
    v[j][0] = lambda[0]*m00 + lambda[1]*m10 + lambda[2]*m20;
    v[j][1] = lambda[0]*m01 + lambda[1]*m11 + lambda[2]*m21;
  }
  return v;
}



double MovingMesh2D::moveDirectionDivergence(const int& n) const
{
  const int& v0 = geometry(DIM,n).vertex(0);
  const int& v1 = geometry(DIM,n).vertex(1);
  const int& v2 = geometry(DIM,n).vertex(2);
  const Point<DIM>& x0 = point(v0);
  const Point<DIM>& x1 = point(v1); 
  const Point<DIM>& x2 = point(v2);
  const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
  return ((m10 - m00)*(x2[1] - x0[1]) - (x1[1] - x0[1])*(m20 - m00) + 
          (x1[0] - x0[0])*(m21 - m01) - (m11 - m01)*(x2[0] - x0[0]))/area;
}



void MovingMesh2D::outputMoveDirection(const std::string& file)
{
  std::ofstream os(file.c_str());
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  for (u_int i = 0;i < n_geometry(0);++ i) {
    os << point(i)[0] << " " << point(i)[1] << " "
       << logical_node[i][0] << " " << logical_node[i][1] << " "
       << move_direction[i][0] << " " << move_direction[i][1] << " "
       << logical_move_direction[i][0] << " " << logical_move_direction[i][1] << "\n";
  }
  os.close();
}

////////////////////////////////////////////////////////////////////////////////

void MovingMesh2D::Solver::Project(const sparse_matrix_t& M, 
                                   const std::vector<int>& bm,
                                   sparse_matrix_t *& P,
                                   sparse_matrix_t *& PMPT,
                                   sparse_matrix_t *& Pt,
                                   std::vector<int> *& pbm)
{
  std::vector<u_int> counter(M.m());
  std::vector<u_int> core_point(M.m());
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const std::size_t * M_rowstart = spM.get_rowstart_indices();
  const u_int * M_colnums = spM.get_column_numbers();

  /**
   * 选择核心点，我们优先将材料标识较高的点选择为核心点。
   *
   */
  u_int n_core_point = 0;
  for (u_int i = 0;i < M.m();i ++) {
    if (counter[i] > 0) continue; /**< 这一行已经被处理过了 */

    std::list<u_int> idx_set;
    idx_set.push_back(i);
    while (!idx_set.empty()) {
      u_int& ii = idx_set.back();
      if (counter[ii] > 0) {
        idx_set.pop_back();
        continue;
      }

      u_int j = M_rowstart[ii];
      for (;j < M_rowstart[ii + 1];j ++) {
        const u_int& k = M_colnums[j];
        if (bm[ii]%bm[k] != 0 &&
            bm[k]%bm[ii] == 0 &&
            counter[k] == 0) break;
      }
      if (j < M_rowstart[ii + 1]) {
        idx_set.push_back(M_colnums[j]);
        continue;
      }

      for (j = M_rowstart[ii];j < M_rowstart[ii + 1];j ++) {
        const u_int& k = M_colnums[j];
        if (bm[ii]%bm[k] == 0) counter[k] ++;
      }
      core_point[n_core_point ++] = ii;
      idx_set.pop_back();
    }
  }

  SparsityPattern& spP = *new SparsityPattern(n_core_point, M.n(), spM.max_entries_per_row());
  SparsityPattern& spPt = *new SparsityPattern(M.n(), n_core_point, spM.max_entries_per_row());
  pbm = new std::vector<int>(n_core_point); /**< 为投影过的自由度的边界标号分配空间 */
  for (u_int i = 0;i < n_core_point;i ++) {
    u_int& k = core_point[i];
    (*pbm)[i] = bm[k]; /**< 设定投影过的自由度的边界标号 */
    for (u_int j = M_rowstart[k];j < M_rowstart[k+1];j ++) {
      const u_int& l = M_colnums[j];
      if (bm[k]%bm[l] != 0) continue; /// 不做聚类的情况
      spP.add(i, l); spPt.add(l, i);
    }
  }
  spP.compress();
  spPt.compress();

  P = new sparse_matrix_t(spP);
  Pt = new sparse_matrix_t(spPt);
  for (u_int i = 0;i < n_core_point;i ++) {
    u_int& k = core_point[i];
    for (u_int j = M_rowstart[k];j < M_rowstart[k+1];j ++) {
      const u_int& l = M_colnums[j];
      if (bm[k]%bm[l] != 0) continue; /// 不做聚类的情况
      P->add(i, l, 1.0/counter[l]);
      Pt->add(l, i, 1.0/counter[l]);
    }
  }

  PMPT = getPMPT(*P, M, *Pt);
}



sparse_matrix_t * 
MovingMesh2D::Solver::getPMPT(const sparse_matrix_t& P, 
                              const sparse_matrix_t& M,
                              const sparse_matrix_t& Pt) const
{
  const SparsityPattern& spP = P.get_sparsity_pattern();
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const SparsityPattern& spPt = Pt.get_sparsity_pattern();
  const std::size_t * P_rowstart = spP.get_rowstart_indices();
  const u_int * P_colnums = spP.get_column_numbers();
  const std::size_t * M_rowstart = spM.get_rowstart_indices();
  const u_int * M_colnums = spM.get_column_numbers();
  const std::size_t * Pt_rowstart = spPt.get_rowstart_indices();
  const u_int * Pt_colnums = spPt.get_column_numbers();

  std::vector<u_int> row_length(P.m(), 0);
  std::vector<bool> flag(P.m(), true);
  std::vector<u_int> index(P.m());
  std::vector<std::vector<u_int> > col_index(P.m());
  for (u_int i = 0;i < P.m();i ++) {
    row_length[i] = 1; /**< add the diagonal entry at first */
    flag[i] = false;
    index[0] = i;
    for (u_int j = P_rowstart[i];j < P_rowstart[i+1];j ++) {
      const u_int& a = P_colnums[j];
      for (u_int k = M_rowstart[a];k < M_rowstart[a+1];k ++) {
        const u_int& b = M_colnums[k];
        for (u_int l = Pt_rowstart[b];l < Pt_rowstart[b+1];l ++) {
          const u_int&  c = Pt_colnums[l];
          if (flag[c]) {
            index[row_length[i] ++] = c;
            flag[c] = false;
          }
        }
      }
    }
    col_index[i].resize(row_length[i]);
    for (u_int j = 0;j < row_length[i];j ++) {
      col_index[i][j] = index[j];
      flag[index[j]] = true;
    }
  }

  SparsityPattern& spA = *(new SparsityPattern(P.m(), row_length));
  for (u_int i = 0;i < P.m();i ++) {
    for (u_int j = 0;j < row_length[i];j ++) {
      spA.add(i, col_index[i][j]);
    }
  }
  spA.compress();

  sparse_matrix_t *A = new sparse_matrix_t(spA);
  lazyPMPT(P, M, Pt, *A);
  return A;
}



void MovingMesh2D::Solver::lazyPMPT(const sparse_matrix_t& P, 
                                    const sparse_matrix_t& M,
                                    const sparse_matrix_t& Pt,
                                    sparse_matrix_t& A) const
{
  const SparsityPattern& spP = P.get_sparsity_pattern();
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const SparsityPattern& spPt = Pt.get_sparsity_pattern();
  const SparsityPattern& spA = A.get_sparsity_pattern();
  const std::size_t * P_rowstart = spP.get_rowstart_indices();
  const u_int * P_colnums = spP.get_column_numbers();
  const std::size_t * M_rowstart = spM.get_rowstart_indices();
  const u_int * M_colnums = spM.get_column_numbers();
  const std::size_t * Pt_rowstart = spPt.get_rowstart_indices();
  const u_int * Pt_colnums = spPt.get_column_numbers();
  const std::size_t * A_rowstart = spA.get_rowstart_indices();
  const u_int * A_colnums = spA.get_column_numbers();

  u_int i, j, k, l;
  std::vector<double> row_entry(P.m(), 0.0);
  for (i = 0;i < P.m();i ++) {
    for (j = P_rowstart[i];j < P_rowstart[i+1];j ++) {
      const u_int& a = P_colnums[j];
      for (k = M_rowstart[a];k < M_rowstart[a+1];k ++) {
        const u_int& b = M_colnums[k];
        for (l = Pt_rowstart[b];l < Pt_rowstart[b+1];l ++) {
          const u_int&  c = Pt_colnums[l];
          row_entry[c] += P.global_entry(j)*M.global_entry(k)*Pt.global_entry(l);
        }
      }
    }
    for (j = A_rowstart[i];j < A_rowstart[i+1];j ++) {
      const u_int& a = A_colnums[j];
      A.global_entry(j) = row_entry[a];
      row_entry[a] = 0.0;
    }
  }
}



void MovingMesh2D::Solver::GaussSidel(const sparse_matrix_t& M, 
                                      std::vector<vector_t>& x, 
                                      const std::vector<vector_t>& r,
                                      const std::vector<int>& bm,
                                      const u_int& s) const
{
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const std::size_t * rowstart = spM.get_rowstart_indices();
  const u_int * colnums = spM.get_column_numbers();
  for (unsigned i = 0;i < s;i ++) {
    for (u_int j = 0;j < M.m();j ++) {
      double rj[DIM] = {r[0](j), r[1](j)};
      for (u_int k = rowstart[j] + 1;k < rowstart[j + 1];k ++) {
        const double& a = M.global_entry(k); /**< M 的元素 */
        rj[0] -= a*x[0](colnums[k]);
        rj[1] -= a*x[1](colnums[k]);
      }
      const double& a = M.global_entry(rowstart[j]); /**< 对角元 */
      if (bm[j] == 1) { /**< 内部的点就不投影了 */
        x[0](j) = rj[0]/a; x[1](j) = rj[1]/a;
        continue;
      }
      rj[0] = rj[0]/a - x[0](j); /**< 首先计算增量 */
      rj[1] = rj[1]/a - x[1](j);
      u_int n_constraint = 0;
      for (u_int k = 0;k < domain->n_edge;k ++) { /**< 对边界上的点逐一进行投影操作 */
        const Edge& edge = domain->edge[k];
        if (bm[j]%edge.boundary_mark != 0) continue; /**< 不在这个表面上 */
        const double * n = edge.logical_normal;
        double ip = rj[0]*n[0] + rj[1]*n[1];
        rj[0] -= ip*n[0]; rj[1] -= ip*n[1];
        n_constraint += 1;
      }
      if (n_constraint == 1) { /// 当只有一个约束的时候才更新，多个约
                               /// 束意味着这个点不允许移动。
        x[0](j) += rj[0]; x[1](j) += rj[1];
      }
    }
  }
}


void MovingMesh2D::Solver::Init(sparse_matrix_t& M,
                                const std::vector<int>& bm,
                                const Domain& d)
{
  std::cerr << "Solver initializing in lazy mode ..." << std::flush;
  u_int k, m;
  projected_matrix.push_back(&M);
  boundary_mark.push_back(&bm);
  domain = &d;
  k = m = M.m();
  n_project = 0;
  while(1) {
#ifdef DEBUG
    std::cerr << "\tlevel " << n_project
              << " with order " << k
              << std::endl;
#endif // DEGUG
    if (k < min_order) break;

    sparse_matrix_t *PMPT, *P, *Pt;
    std::vector<int> *pbm;
    Project(*projected_matrix[n_project], *boundary_mark[n_project],
            P, PMPT, Pt, pbm);
    m = PMPT->m();
    if (k < 2*m) {
      const SparsityPattern& spP = P->get_sparsity_pattern();
      delete P; delete (&spP);
      const SparsityPattern& spPt = Pt->get_sparsity_pattern();
      delete Pt; delete (&spPt);
      const SparsityPattern& spM = PMPT->get_sparsity_pattern();
      delete PMPT; delete (&spM);
      delete pbm;
      break;
    }
    k = m;
    project_matrix.push_back(P);
    project_matrix_r.push_back(Pt);
    projected_matrix.push_back(PMPT);
    boundary_mark.push_back(pbm);
    n_project ++;
  };
  is_initialized = true;
  std::cerr << " OK! grid levels: " << n_project << std::endl;
}



/**
 * to solve a linear system with the same sparsity pattern as already initialized
 * matrix
 * 
 */
void MovingMesh2D::Solver::reinit(sparse_matrix_t& M,
                                  const std::vector<int>& bm,
                                  const Domain& d)
{
  if (!is_initialized) {
    Init(M, bm, d);
    return;
  }

  std::cerr << "Solver reinitializing in lazy mode ..." << std::flush;
  projected_matrix[0] = &M;
  for (u_int j = 0;j < n_project;j ++) {
    lazyPMPT(*project_matrix[j], 
             *projected_matrix[j], 
             *project_matrix_r[j],
             *projected_matrix[j + 1]);
  };
  std::cerr << " OK! grid levels: " << n_project << std::endl;
}



void MovingMesh2D::Solver::solve(std::vector<vector_t>& x, 
                                 const std::vector<vector_t>& r,
                                 u_int steps) const
{
  int i;
  u_int j, k, n_iter = 0;
  std::vector<std::vector<vector_t> *> projected_r(n_project + 1);
  std::vector<std::vector<vector_t> *> projected_x(n_project + 1);
  std::vector<vector_t> r1(r);
  projected_x[0] = &x;
  projected_r[0] = &r1;

  for (j = 1;j <= n_project;j ++) {
    projected_r[j] = new std::vector<vector_t>
      (DIM, vector_t(projected_matrix[j]->m()));
    projected_x[j] = new std::vector<vector_t>
      (DIM, vector_t(projected_matrix[j]->m()));
  }
  double init_residual = 0, residual;
  vector_t v0(projected_matrix[0]->m());
  for (j = 0;j < DIM;j ++) {
    projected_matrix[0]->vmult(v0, x[j]);
    v0.add(-1, r[j]);
    for (k = 0;k < x[j].size();k ++) { /**< 清除边界能量 */
      v0(k) *= ((*boundary_mark[0])[k] == 0);
    }
    init_residual += v0.l2_norm();
  }
  std::cerr << "Solver begin with initial residual " 
            << init_residual << " ..." << std::endl;
  while(1) {
    for (i = 0;i < n_project;i ++) {
      GaussSidel(*projected_matrix[i], *projected_x[i], *projected_r[i], 
                 *boundary_mark[i], smooth_step);
      v0.reinit(projected_matrix[i]->m());
      for (j = 0;j < DIM;j ++) {
        projected_matrix[i]->vmult(v0, (*projected_x[i])[j]);
        v0.sadd(-1, (*projected_r[i])[j]);
        project_matrix[i]->vmult((*projected_r[i+1])[j], v0);
        (*projected_x[i+1])[j] = 0;
      }
    }
		
    GaussSidel(*projected_matrix[i], *projected_x[i], *projected_r[i], 
               *boundary_mark[i], smooth_step);

    for (i = n_project - 1;i >= 0;i --) {
      v0.reinit(projected_matrix[i]->m());
      for (j = 0;j < DIM;j ++) {
        project_matrix[i]->Tvmult(v0, (*projected_x[i+1])[j]);
        (*projected_x[i])[j] += v0;
      }
      GaussSidel(*projected_matrix[i], *projected_x[i], *projected_r[i], 
                 *boundary_mark[i], smooth_step);
    }
	
    residual = 0;
    v0.reinit(projected_matrix[0]->m());
    for (j = 0;j < DIM;j ++) {
      projected_matrix[0]->vmult(v0, x[j]);
      v0.add(-1, r[j]);
      for (k = 0;k < x[j].size();k ++) { /**< 清除边界能量 */
        v0(k) *= ((*boundary_mark[0])[k] == 0);
      }
      residual += v0.l2_norm();
    }
#ifdef DEBUG
    std::cerr << "\tresidual = " << residual << std::endl;
#endif // DEBUG
    if ((++ n_iter) > steps) break;
  };
  for (j = 1;j <= n_project;j ++) {
    delete projected_r[j];
    delete projected_x[j];
  }
  std::cerr << "Solver exit with residual as " 
            << residual << std::endl;
}



void MovingMesh2D::Solver::clear()
{	
  if (!is_initialized) return;
  for (u_int i = 0;i < n_project;i ++) {
    const SparsityPattern& spP = project_matrix[i]->get_sparsity_pattern();
    delete project_matrix[i]; delete (&spP);
    const SparsityPattern& spM = projected_matrix[i+1]->get_sparsity_pattern();
    delete projected_matrix[i+1]; delete (&spM);
    const SparsityPattern& spPt = project_matrix_r[i]->get_sparsity_pattern();
    delete project_matrix_r[i]; delete (&spPt);
    delete boundary_mark[i + 1];
  }
  project_matrix.clear();
  projected_matrix.clear();
  project_matrix_r.clear();
  boundary_mark.clear();
  is_initialized = false;
}



MovingMesh2D::Solver::Solver() : 
  is_initialized(false),
  n_project(0),
  smooth_step(3),
  min_order(50)
{}



MovingMesh2D::Solver::~Solver()
{
  clear ();
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
