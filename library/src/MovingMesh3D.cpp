/**
 * @file   MovingMesh3D.cpp
 * @author Yana Di, Ruo Li
 * @date   Wed Dec 21 20:16:26 2005
 * 
 * @brief  the implementation of class MovingMesh3D
 * 
 * 
 */

#include "MovingMesh3D.h"

AFEPACK_OPEN_NAMESPACE

#define DIM 3

typedef SparseMatrix<double> sparse_matrix_t;
typedef Vector<double> vector_t;

/// 素数表
int MovingMesh3D::primes[] = {
  02, 03, 05, 07, 11, 
  13, 17, 19, 23, 29, 
  31, 37, 41, 43, 47,
  53, 59, 61, 67, 71,
  73, 79, 83, 89, 91
};

MovingMesh3D::MovingMesh3D() :
  tol(0.1), solve_step(5)
{}



MovingMesh3D::~MovingMesh3D()
{}



void MovingMesh3D::readDomain(const std::string& filename)
{
  u_int i, j, k, l0, l1;

  readData(filename);
  logical_node.resize(n_geometry(0));
  move_direction.resize(n_geometry(0));
  logical_move_direction.resize(n_geometry(0));
  monitor().resize(n_geometry(DIM));

  std::ifstream is;
  is.open((filename + ".d").c_str());
  if (!is) { 
    std::cout << "Open domain description file "
	      << filename + ".d failure, aborting... "
              << std::endl;
    abort();
  }
  readDummy(is); is >> domain.n_vertex;
  domain.physical_domain_vertex.resize(domain.n_vertex);
  domain.logical_domain_vertex.resize(domain.n_vertex);
  for (i = 0;i < domain.n_vertex;i ++) {
    readDummy(is); is >> domain.physical_domain_vertex[i].index;
    for (j = 0;j < 3;j ++) {
      readDummy(is); is >> domain.physical_domain_vertex[i][j];
    }
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
    domain.edge[i].boundary_mark = 1; /// 缺省值设为 1
  }
  readDummy(is); is >> domain.n_surf;
  domain.surface.resize(domain.n_surf);
  for (i = 0;i < domain.n_surf;i ++) {
    readDummy(is); is >> domain.surface[i].index;
    readDummy(is); is >> domain.surface[i].n_edge;
    domain.surface[i].edge.resize(domain.surface[i].n_edge);
    for (j = 0; j < domain.surface[i].n_edge; j ++) {
      readDummy(is); is >> domain.surface[i].edge[j];
    }
    readDummy(is); is >> domain.surface[i].boundary_mark;
    for (j = 0;j < DIM; j ++) { /**< 读入物理区域的表面法向 */
      readDummy(is); is >> domain.surface[i].normal[j];
    }
    double * n = domain.surface[i].normal; /**< 对表面法向向量做单位化 */
    double ln = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    n[0] /= ln; n[1] /= ln; n[2] /= ln;

    for (j = 0;j < DIM; j ++) { /**< 读入逻辑区域的表面法向 */
      readDummy(is); is >> domain.surface[i].logic_normal[j];
    }
    n = domain.surface[i].logic_normal; /**< 对表面法向向量做单位化 */
    ln = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    n[0] /= ln; n[1] /= ln; n[2] /= ln;
  }
  is.close();
	
  is.open((filename + ".log").c_str());
  if (!is) { 
    std::cerr << "Open the logical domain description file "
              << filename + ".log failure, aborting ... " 
              << std::endl;
  }
  for (i = 0;i < domain.n_vertex;i ++) {
    for (j = 0;j < DIM;j ++) {
      is >> domain.logical_domain_vertex[i][j];
    }
  }
  is.close();

  parseBoundary();

  std::vector<u_int> n_coupling_node(n_geometry(0), 1);
  for (i = 0;i < n_geometry(1);i ++) {
    n_coupling_node[geometry(1,i).vertex(0)] ++;
    n_coupling_node[geometry(1,i).vertex(1)] ++;
  }
  spM.reinit(n_geometry(0), n_geometry(0), n_coupling_node);
  for (i = 0;i < n_geometry(1);i ++) { /**< 对每个边做循环 */
    l0 = geometry(1,i).vertex(0);
    l1 = geometry(1,i).vertex(1);
    spM.add(l0, l0); spM.add(l0, l1);
    spM.add(l1, l0); spM.add(l1, l1);
  }
  spM.compress();

  getLogicalMesh();
}



void MovingMesh3D::moveMesh()
{
  u_int i, j;
  double epsilon = tolerence();
  double error = 2*epsilon;
  while (error > epsilon) {
    getMoveDirection();
    error = 0;
    for (i = 0;i < n_geometry(0);i ++) {
      double d = 0.;
      for (j = 0;j < DIM;j ++) {
	d += move_direction[i][j]*move_direction[i][j];
      }
      d = sqrt(d);
      if (error < d)
	error = d;
    }
#ifdef DEBUG
    std::cerr << "mesh moving error = " << error << std::endl;
#endif // DEBUG
    getMoveStepLength();
    for (i = 0;i < n_move_step;i ++) {
      updateSolution();
      updateMesh();
    }
  };
}



void MovingMesh3D::outputPhysicalMesh(const std::string& filename)
{
  writeData(filename);
}



void MovingMesh3D::outputLogicalMesh(const std::string& filename)
{
  Mesh<DIM,DIM> mesh = Mesh<DIM,DIM>(*this);
  for (u_int i = 0;i < n_geometry(0);i ++)
    mesh.point(i) = logical_node[i];
  mesh.writeData(filename);
}



void MovingMesh3D::getMoveDirection()
{
  u_int i, j, k, l, l0, l1;
  getMonitor();
  M.reinit(spM);
  for (i = 0;i < n_geometry(DIM);i ++) {
    const Point<DIM>& x0 = point(geometry(DIM,i).vertex(0));
    const Point<DIM>& x1 = point(geometry(DIM,i).vertex(1)); 
    const Point<DIM>& x2 = point(geometry(DIM,i).vertex(2));
    const Point<DIM>& x3 = point(geometry(DIM,i).vertex(3));
    double omega[3][4];
    double volume =   ((x1[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
                       (x1[1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
                       (x1[2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
                       (x1[0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
                       (x1[1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
                       (x1[2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0]));
    omega[0][0] = -((x2[1] - x1[1])*(x3[2] - x1[2]) - (x2[2] - x1[2])*(x3[1] - x1[1]));
    omega[0][1] =  ((x3[1] - x2[1])*(x0[2] - x2[2]) - (x3[2] - x2[2])*(x0[1] - x2[1]));
    omega[0][2] = -((x0[1] - x3[1])*(x1[2] - x3[2]) - (x0[2] - x3[2])*(x1[1] - x3[1]));
    omega[0][3] =  ((x1[1] - x0[1])*(x2[2] - x0[2]) - (x1[2] - x0[2])*(x2[1] - x0[1]));
    omega[1][0] = -((x2[2] - x1[2])*(x3[0] - x1[0]) - (x2[0] - x1[0])*(x3[2] - x1[2]));
    omega[1][1] =  ((x3[2] - x2[2])*(x0[0] - x2[0]) - (x3[0] - x2[0])*(x0[2] - x2[2]));
    omega[1][2] = -((x0[2] - x3[2])*(x1[0] - x3[0]) - (x0[0] - x3[0])*(x1[2] - x3[2]));
    omega[1][3] =  ((x1[2] - x0[2])*(x2[0] - x0[0]) - (x1[0] - x0[0])*(x2[2] - x0[2]));
    omega[2][0] = -((x2[0] - x1[0])*(x3[1] - x1[1]) - (x2[1] - x1[1])*(x3[0] - x1[0]));
    omega[2][1] =  ((x3[0] - x2[0])*(x0[1] - x2[1]) - (x3[1] - x2[1])*(x0[0] - x2[0]));
    omega[2][2] = -((x0[0] - x3[0])*(x1[1] - x3[1]) - (x0[1] - x3[1])*(x1[0] - x3[0]));
    omega[2][3] =  ((x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]));
    for (j = 0;j < 4;j ++) {
      l0 = geometry(DIM,i).vertex(j);
      for (k = 0;k < 4;k ++) {
	double d = monitor(i)*(omega[0][j]*omega[0][k] +
			       omega[1][j]*omega[1][k] + 
			       omega[2][j]*omega[2][k])/volume;
	l1 = geometry(DIM,i).vertex(k);
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

  std::vector<double> mass_lumping(n_geometry(0), 0.);
  for (i = 0;i < n_geometry(DIM);i ++) {
    const int& v0 = geometry(DIM,i).vertex(0);
    const int& v1 = geometry(DIM,i).vertex(1);
    const int& v2 = geometry(DIM,i).vertex(2);
    const int& v3 = geometry(DIM,i).vertex(3);
    const Point<DIM>& x0 = point(v0); const Point<DIM>& xi0 = logical_node[v0];
    const Point<DIM>& x1 = point(v1); const Point<DIM>& xi1 = logical_node[v1];
    const Point<DIM>& x2 = point(v2); const Point<DIM>& xi2 = logical_node[v2];
    const Point<DIM>& x3 = point(v3); const Point<DIM>& xi3 = logical_node[v3];
    double jacobi[3][3];
    double volume = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
                     (xi1[1] - xi0[1])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
                     (xi1[2] - xi0[2])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
                     (xi1[0] - xi0[0])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
                     (xi1[1] - xi0[1])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
                     (xi1[2] - xi0[2])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));
    jacobi[0][0] = ((x1[0] - x0[0])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
                    (xi1[1] - xi0[1])*(xi2[2] - xi0[2])*(x3[0] - x0[0]) +
                    (xi1[2] - xi0[2])*(x2[0] - x0[0])*(xi3[1] - xi0[1]) -
                    (x1[0] - x0[0])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
                    (xi1[1] - xi0[1])*(x2[0] - x0[0])*(xi3[2] - xi0[2]) -
                    (xi1[2] - xi0[2])*(xi2[1] - xi0[1])*(x3[0] - x0[0]));
    jacobi[1][0] = ((x1[1] - x0[1])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
                    (xi1[1] - xi0[1])*(xi2[2] - xi0[2])*(x3[1] - x0[1]) +
                    (xi1[2] - xi0[2])*(x2[1] - x0[1])*(xi3[1] - xi0[1]) -
                    (x1[1] - x0[1])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
                    (xi1[1] - xi0[1])*(x2[1] - x0[1])*(xi3[2] - xi0[2]) -
                    (xi1[2] - xi0[2])*(xi2[1] - xi0[1])*(x3[1] - x0[1]));
    jacobi[2][0] = ((x1[2] - x0[2])*(xi2[1] - xi0[1])*(xi3[2] - xi0[2]) +
                    (xi1[1] - xi0[1])*(xi2[2] - xi0[2])*(x3[2] - x0[2]) +
                    (xi1[2] - xi0[2])*(x2[2] - x0[2])*(xi3[1] - xi0[1]) -
                    (x1[2] - x0[2])*(xi2[2] - xi0[2])*(xi3[1] - xi0[1]) -
                    (xi1[1] - xi0[1])*(x2[2] - x0[2])*(xi3[2] - xi0[2]) -
                    (xi1[2] - xi0[2])*(xi2[1] - xi0[1])*(x3[2] - x0[2]));
    jacobi[0][1] = ((xi1[0] - xi0[0])*(x2[0] - x0[0])*(xi3[2] - xi0[2]) +
                    (x1[0] - x0[0])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
                    (xi1[2] - xi0[2])*(xi2[0] - xi0[0])*(x3[0] - x0[0]) -
                    (xi1[0] - xi0[0])*(xi2[2] - xi0[2])*(x3[0] - x0[0]) -
                    (x1[0] - x0[0])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
                    (xi1[2] - xi0[2])*(x2[0] - x0[0])*(xi3[0] - xi0[0]));
    jacobi[1][1] = ((xi1[0] - xi0[0])*(x2[1] - x0[1])*(xi3[2] - xi0[2]) +
                    (x1[1] - x0[1])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
                    (xi1[2] - xi0[2])*(xi2[0] - xi0[0])*(x3[1] - x0[1]) -
                    (xi1[0] - xi0[0])*(xi2[2] - xi0[2])*(x3[1] - x0[1]) -
                    (x1[1] - x0[1])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
                    (xi1[2] - xi0[2])*(x2[1] - x0[1])*(xi3[0] - xi0[0]));
    jacobi[2][1] = ((xi1[0] - xi0[0])*(x2[2] - x0[2])*(xi3[2] - xi0[2]) +
                    (x1[2] - x0[2])*(xi2[2] - xi0[2])*(xi3[0] - xi0[0]) +
                    (xi1[2] - xi0[2])*(xi2[0] - xi0[0])*(x3[2] - x0[2]) -
                    (xi1[0] - xi0[0])*(xi2[2] - xi0[2])*(x3[2] - x0[2]) -
                    (x1[2] - x0[2])*(xi2[0] - xi0[0])*(xi3[2] - xi0[2]) -
                    (xi1[2] - xi0[2])*(x2[2] - x0[2])*(xi3[0] - xi0[0]));
    jacobi[0][2] = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*(x3[0] - x0[0]) +
                    (xi1[1] - xi0[1])*(x2[0] - x0[0])*(xi3[0] - xi0[0]) +
                    (x1[0] - x0[0])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
                    (xi1[0] - xi0[0])*(x2[0] - x0[0])*(xi3[1] - xi0[1]) -
                    (xi1[1] - xi0[1])*(xi2[0] - xi0[0])*(x3[0] - x0[0]) -
                    (x1[0] - x0[0])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));
    jacobi[1][2] = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*(x3[1] - x0[1]) +
                    (xi1[1] - xi0[1])*(x2[1] - x0[1])*(xi3[0] - xi0[0]) +
                    (x1[1] - x0[1])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
                    (xi1[0] - xi0[0])*(x2[1] - x0[1])*(xi3[1] - xi0[1]) -
                    (xi1[1] - xi0[1])*(xi2[0] - xi0[0])*(x3[1] - x0[1]) -
                    (x1[1] - x0[1])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));
    jacobi[2][2] = ((xi1[0] - xi0[0])*(xi2[1] - xi0[1])*(x3[2] - x0[2]) +
                    (xi1[1] - xi0[1])*(x2[2] - x0[2])*(xi3[0] - xi0[0]) +
                    (x1[2] - x0[2])*(xi2[0] - xi0[0])*(xi3[1] - xi0[1]) -
                    (xi1[0] - xi0[0])*(x2[2] - x0[2])*(xi3[1] - xi0[1]) -
                    (xi1[1] - xi0[1])*(xi2[0] - xi0[0])*(x3[2] - x0[2]) -
                    (x1[2] - x0[2])*(xi2[1] - xi0[1])*(xi3[0] - xi0[0]));
    for (j = 0;j < 4;j ++) {
      k = geometry(DIM,i).vertex(j);
      for (l0 = 0;l0 < DIM;++ l0) {
        for (l1 = 0;l1 < DIM;++ l1) {
          move_direction[k][l0] += jacobi[l0][l1]*logical_move_direction[k][l1];
        }
      }
      mass_lumping[k] += volume;
    }
  }
  for (i = 0;i < n_geometry(0);i ++) {
    for (k = 0;k < DIM;++ k) {
      move_direction[i][k] /= mass_lumping[i];
    }
  }

  /// 下面这一段是有问题的，Todo: 已经修正
  for (j = 0;j < n_geometry(0);j ++) {
    if (boundary_mark[j] == 1) continue;

    u_int n_constraint = 0;
    std::vector<u_int> constraint_idx(domain.n_surf);
    for (k = 0;k < domain.n_surf;k ++) { /// 计算约束的个数
      const Surface& surface = domain.surface[k];
      if (boundary_mark[j]%surface.boundary_mark != 0) continue; /**< 不在这个表面上 */
      constraint_idx[n_constraint ++] = k;
    }
    if (n_constraint == 1) { /// 在表面上
      const Surface& surface = domain.surface[constraint_idx[0]];
      const double * n = surface.normal;
      double ip = (move_direction[j][0]*n[0] + 
                   move_direction[j][1]*n[1] + 
                   move_direction[j][2]*n[2]);
      move_direction[j][0] -= ip*n[0]; 
      move_direction[j][1] -= ip*n[1]; 
      move_direction[j][2] -= ip*n[2]; 
    }
    else if (n_constraint == 2) { /// 在棱上
      const Surface& surface0 = domain.surface[constraint_idx[0]];
      const double * n0 = surface0.normal;
      const Surface& surface1 = domain.surface[constraint_idx[1]];
      const double * n1 = surface1.normal;
      double n[3] = {n0[1]*n1[2] - n0[2]*n1[1],
                     n0[2]*n1[0] - n0[0]*n1[2],
                     n0[0]*n1[1] - n0[1]*n1[0]};
      double ip = (move_direction[j][0]*n[0] + 
                   move_direction[j][1]*n[1] + 
                   move_direction[j][2]*n[2])/(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
      move_direction[j][0] = ip*n[0]; 
      move_direction[j][1] = ip*n[1]; 
      move_direction[j][2] = ip*n[2]; 
    }
    else { /// 在顶点上
      move_direction[j][0] = 0.0; 
      move_direction[j][1] = 0.0; 
      move_direction[j][2] = 0.0;
    }
  }
}



void MovingMesh3D::getMoveStepLength()
{
  u_int i, j, k;
  double a[4], pi = 4*atan(1.0);
  n_move_step = 1;
  move_step_length = 1;
  for (i = 0; i < n_geometry(DIM); i ++) {
    const GeometryBM& geo = geometry(DIM,i);
    double x[DIM][DIM], mx[DIM][DIM];
    for (j = 0;j < DIM;j ++) {
      for (k = 0;k < DIM;k ++) {
	x[j][k] = point(geo.vertex(j+1))[k] - point(geo.vertex(0))[k];
	mx[j][k] = move_direction[geo.vertex(j+1)][k] - move_direction[geo.vertex(0)][k];
      }
    }

#define DETERMINANT3(r0, r1, r2) \
   (r0[0]*r1[1]*r2[2] + r0[1]*r1[2]*r2[0] + r0[2]*r1[0]*r2[1] \
  - r0[0]*r1[2]*r2[1] - r0[2]*r1[1]*r2[0] - r0[1]*r1[0]*r2[2])
    a[3] = DETERMINANT3(mx[0], mx[1], mx[2]);
    a[2] = (DETERMINANT3( x[0], mx[1], mx[2]) +
            DETERMINANT3(mx[0],  x[1], mx[2]) +
            DETERMINANT3(mx[0], mx[1],  x[2]));
    a[1] = (DETERMINANT3(mx[0],  x[1],  x[2]) +
            DETERMINANT3( x[0], mx[1],  x[2]) +
            DETERMINANT3( x[0],  x[1], mx[2]));
    a[0] = DETERMINANT3( x[0],  x[1],  x[2]);
#undef DETERMINANT3
    a[0] /= a[3]; a[1] /= a[3]; a[2] /= a[3]; a[3] = 1.0;

    double Q = (3*a[1] - a[2]*a[2])/9;
    double R = (9*a[1]*a[2] - 27*a[0] - 2*a[2]*a[2]*a[2])/54;
    double D = R*R + Q*Q*Q;
    if (D > 0) { /**< 此时只有一个实根 */
      double S = cbrt(R + sqrt(D)), T = cbrt(R - sqrt(D));
      double root = -a[2]/3 + (S + T);
      if (root > 0 && move_step_length > root) {
	move_step_length = root;
      }
    }
    else { /**< 三个实根的情形 */
      double theta = acos(R/sqrt(-Q*Q*Q));
      double root[3] = {2*sqrt(-Q)*cos(theta/3) - a[2]/3,
			2*sqrt(-Q)*cos((theta + 2*pi)/3) - a[2]/3,
			2*sqrt(-Q)*cos((theta - 2*pi)/3) - a[2]/3};
      for (j = 0;j < 3;j ++) {
	if (root[j] > 0 && move_step_length > root[j]) {
	  move_step_length = root[j];
	}
      }
    }
  }
  move_step_length *= 0.5;
  std::cout << "move step length = " << move_step_length << std::endl;
}



void MovingMesh3D::readDummy(std::ifstream& is)
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



void MovingMesh3D::getMonitor()
{
  std::fill(monitor().begin(), monitor().end(), 1);
}



void MovingMesh3D::smoothMonitor(u_int s)
{
  u_int i, j, k;
  std::vector<float> volume(n_geometry(DIM));
  std::vector<float> mass_lumping(n_geometry(0), 0);
  std::vector<float> monitor1(n_geometry(0));
  for (i = 0;i < n_geometry(DIM);i ++) {
    const Point<DIM>& x0 = point(geometry(DIM,i).vertex(0));
    const Point<DIM>& x1 = point(geometry(DIM,i).vertex(1)); 
    const Point<DIM>& x2 = point(geometry(DIM,i).vertex(2));
    const Point<DIM>& x3 = point(geometry(DIM,i).vertex(3));
    volume[i] = ((x1[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
                 (x1[1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
                 (x1[2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
                 (x1[0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
                 (x1[1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
                 (x1[2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0]));
    for (j = 0;j < 4;j ++)
      mass_lumping[geometry(DIM,i).vertex(j)] += volume[i];
  }
  for (i = 0;i < s;i ++) {
    std::fill(monitor1.begin(), monitor1.end(), 0);
    for (j = 0;j < n_geometry(DIM);j ++) {
      for (k = 0;k < 4;k ++) {
	monitor1[geometry(DIM,j).vertex(k)] += monitor(j)*volume[j];
      }
    }
    for (j = 0;j < n_geometry(0);j ++)
      monitor1[j] /= 4*mass_lumping[j];
    std::fill(monitor().begin(), monitor().end(), 0);
    for (j = 0;j < n_geometry(DIM);j ++) {
      for (k = 0;k < 4;k ++) {
	monitor(j) += monitor1[geometry(DIM,j).vertex(k)];
      }
    }
  }
}



void MovingMesh3D::updateMesh()
{
  for (u_int i = 0;i < n_geometry(0);i ++) {
    for (u_int k = 0;k < DIM;++ k) {
      point(i)[k] += move_step_length * move_direction[i][k];
    }
  }
}



void MovingMesh3D::getLogicalMesh()
{
  std::cout << "Computing logical mesh ..." << std::endl;
  u_int i, j, k, l0, l1;
  bmark_t bm;
  for (i = 0;i < n_geometry(0);i ++) {
    logical_node[i][0] = point(i)[0];
    logical_node[i][1] = point(i)[1];
    logical_node[i][2] = point(i)[2];
  }
  //   for (i = 0;i < n_boundary_node;i ++) {
  //     j = boundary_node_index[i];
  //     bm = geometry(0,j).boundaryMark();
  //     for (k = 0;k < domain.n_edge;k ++) {
  //       if (domain.edge[k].boundary_mark == bm) {
  // 	break;
  //       }
  //     }
  //     double p[2][2];
  //     for (l0 = 0;l0 < domain.n_vertex;l0 ++) {
  //       if(domain.physical_domain_vertex[l0].index == domain.edge[k].vertex[0]) {
  // 	break;
  //       }
  //     }
  //     p[0][0] = domain.physical_domain_vertex[l0][0];
  //     p[0][1] = domain.physical_domain_vertex[l0][1];
  //     for (l1 = 0;l1 < domain.n_vertex;l1 ++) {
  //       if(domain.physical_domain_vertex[l1].index == domain.edge[k].vertex[1]) {
  // 	break;
  //       }
  //     }
  //     p[1][0] = domain.physical_domain_vertex[l1][0];
  //     p[1][1] = domain.physical_domain_vertex[l1][1];
  //     if (fabs(domain.physical_domain_vertex[l0][0] - domain.physical_domain_vertex[l1][0]) > 
  // 	fabs(domain.physical_domain_vertex[l0][1] - domain.physical_domain_vertex[l1][1]))
  //       k = 0;
  //     else
  //       k = 1;
  //     double lambda = (point(j)[k] - domain.physical_domain_vertex[l0][k]) 
  //       / (domain.physical_domain_vertex[l1][k] - domain.physical_domain_vertex[l0][k]);
  //     logical_node[j][0] = (1-lambda)*domain.logical_domain_vertex[l0][0] + lambda*domain.logical_domain_vertex[l1][0];
  //     logical_node[j][1] = (1-lambda)*domain.logical_domain_vertex[l0][1] + lambda*domain.logical_domain_vertex[l1][1];
  //   }
  //   sparse_matrix_t M(spM);
  //   sparse_matrix_t N(spN);
  //   for (i = 0;i < n_geometry(2);i ++) {
  //     const Point<2>& x0 = point(geometry(2,i).vertex(0));
  //     const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
  //     const Point<2>& x2 = point(geometry(2,i).vertex(2));
  //     double omega[2][3];
  //     double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
  //     omega[0][0] = x2[0] - x1[0];
  //     omega[0][1] = x0[0] - x2[0];
  //     omega[0][2] = x1[0] - x0[0];
  //     omega[1][0] = x2[1] - x1[1];
  //     omega[1][1] = x0[1] - x2[1];
  //     omega[1][2] = x1[1] - x0[1];
  //     for (j = 0;j < 3;j ++) {
  //       l0 = geometry(2,i).vertex(j);
  //       if (geometry(0,l0).boundaryMark() == 0) {
  // 	for (k = 0;k < 3;k ++) {
  // 	  double d = (omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k])/area;
  // 	  l1 = geometry(2,i).vertex(k);
  // 	  if (geometry(0,l1).boundaryMark())
  // 	    N.add(index[l0], index[l1], d);
  // 	  else
  // 	    M.add(index[l0], index[l1], d);
  // 	}
  //       }
  //     }
  //   }
  //   vector_t b(n_boundary_node);
  //   vector_t r(n_interior_node);
  //   vector_t x(n_interior_node);
  //   for (i = 0;i < n_boundary_node;i ++) {
  //     b(i) = -logical_node[boundary_node_index[i]][0];
  //   }
  //   N.vmult(r, b);
  //   AMGSolver solver(M);
  //   solver.solve(x, r);
  //   for (i = 0;i < n_interior_node;i ++) {
  //     logical_node[interior_node_index[i]][0] = x(i);
  //   }

  //   for (i = 0;i < n_boundary_node;i ++) {
  //     b(i) = -logical_node[boundary_node_index[i]][1];
  //   }
  //   N.vmult(r, b);
  //   solver.solve(x, r);
  //   for (i = 0;i < n_interior_node;i ++) {
  //     logical_node[interior_node_index[i]][1] = x(i);
  //   }
}



void MovingMesh3D::parseBoundary()
{
  std::cout << "Parsing boundary nodes and faces ..." << std::endl;

  u_int n_node = n_geometry(0);
  u_int n_surf = n_geometry(DIM-1);

  /// 将第 j 个面上的材料标识设置为第 j 个素数
  std::vector<int> surf_idx(n_surf); /// 网格的每个面位于区域的面的序号
  for (u_int i = 0;i < n_surf;++ i) { /// 对所有的面做循环
    const GeometryBM& surf = geometry(DIM-1, i);
    const bound_t& bm = surf.boundaryMark();
    if (bm == 0) { /// 寻找位于边界上的面，而这个面在内部
      surf_idx[i] = -1;
      continue; 
    }
    for (u_int j = 0;j < domain.n_surf;++ j) { /// 否则搜索
      if (bm == domain.surface[j].boundary_mark) {
        surf_idx[i] = j;
        break;
      }
    }
  }

  /// 将边界的材料标识设为第 i 个素数
  for (u_int i = 0;i < domain.n_surf;++ i) {
    domain.surface[i].boundary_mark = MovingMesh3D::primes[i]; 
  }

  /// 统计每个区域的面上的节点的个数，并为每个节点设置正确的材料标识
  boundary_mark.resize(n_node, 1);
  for (u_int i = 0;i < n_surf;++ i) { /// 对所有的边做循环
    const GeometryBM& surf = geometry(DIM-1,i);
    const int& idx = surf_idx[i];
    if (idx == -1) continue; /// 这条边在内部，我们跳过去
    const int& bm = domain.surface[idx].boundary_mark;
    for (u_int j = 0;j < 3;++ j) { /// 考虑其三个顶点
      const u_int& v = surf.vertex(j);
      if (boundary_mark[v]%bm != 0) { /// 这个端点应该还没有考虑过
        boundary_mark[v] *= bm; /// 在材料标识上乘以本边的标识
      }
    }
  }        
}



std::vector<double> 
MovingMesh3D::moveDirection(const Point<DIM>& p, 
                            const int& n) const
{
  const int& v0 = geometry(DIM,n).vertex(0);
  const int& v1 = geometry(DIM,n).vertex(1);
  const int& v2 = geometry(DIM,n).vertex(2);
  const int& v3 = geometry(DIM,n).vertex(3);

  const Point<DIM>& x0 = point(v0);
  const Point<DIM>& x1 = point(v1); 
  const Point<DIM>& x2 = point(v2);
  const Point<DIM>& x3 = point(v3);

  const double& m00 = move_direction[v0][0]; 
  const double& m01 = move_direction[v0][1]; 
  const double& m02 = move_direction[v0][2];

  const double& m10 = move_direction[v1][0]; 
  const double& m11 = move_direction[v1][1]; 
  const double& m12 = move_direction[v1][2];

  const double& m20 = move_direction[v2][0]; 
  const double& m21 = move_direction[v2][1]; 
  const double& m22 = move_direction[v2][2];

  const double& m30 = move_direction[v3][0]; 
  const double& m31 = move_direction[v3][1]; 
  const double& m32 = move_direction[v3][2];

  double volume = ((x1[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
                   (x1[1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
                   (x1[2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
                   (x1[0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
                   (x1[1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
                   (x1[2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0]));
  double lambda[4];
  lambda[0] = ((x1[0] - p[0])*(x2[1] - p[1])*(x3[2] - p[2]) +
               (x1[1] - p[1])*(x2[2] - p[2])*(x3[0] - p[0]) +
               (x1[2] - p[2])*(x2[0] - p[0])*(x3[1] - p[1]) -
               (x1[0] - p[0])*(x2[2] - p[2])*(x3[1] - p[1]) -
               (x1[1] - p[1])*(x2[0] - p[0])*(x3[2] - p[2]) -
               (x1[2] - p[2])*(x2[1] - p[1])*(x3[0] - p[0]))/volume;
  lambda[1] = ((p[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
               (p[1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
               (p[2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
               (p[0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
               (p[1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
               (p[2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0]))/volume;
  lambda[2] = ((x1[0] - x0[0])*(p[1] - x0[1])*(x3[2] - x0[2]) +
               (x1[1] - x0[1])*(p[2] - x0[2])*(x3[0] - x0[0]) +
               (x1[2] - x0[2])*(p[0] - x0[0])*(x3[1] - x0[1]) -
               (x1[0] - x0[0])*(p[2] - x0[2])*(x3[1] - x0[1]) -
               (x1[1] - x0[1])*(p[0] - x0[0])*(x3[2] - x0[2]) -
               (x1[2] - x0[2])*(p[1] - x0[1])*(x3[0] - x0[0]))/volume;
  lambda[3] = ((x1[0] - x0[0])*(x2[1] - x0[1])*(p[2] - x0[2]) +
               (x1[1] - x0[1])*(x2[2] - x0[2])*(p[0] - x0[0]) +
               (x1[2] - x0[2])*(x2[0] - x0[0])*(p[1] - x0[1]) -
               (x1[0] - x0[0])*(x2[2] - x0[2])*(p[1] - x0[1]) -
               (x1[1] - x0[1])*(x2[0] - x0[0])*(p[2] - x0[2]) -
               (x1[2] - x0[2])*(x2[1] - x0[1])*(p[0] - x0[0]))/volume;
  std::vector<double> v(DIM);
  v[0] = lambda[0]*m00 + lambda[1]*m10 + lambda[2]*m20 + lambda[3]*m30;
  v[1] = lambda[0]*m01 + lambda[1]*m11 + lambda[2]*m21 + lambda[3]*m31;
  v[2] = lambda[0]*m02 + lambda[1]*m12 + lambda[2]*m22 + lambda[3]*m32;
  return v;
}




std::vector<std::vector<double> > 
MovingMesh3D::moveDirection(const std::vector<Point<DIM> >& p, 
                            const int& n) const
{
  const int& v0 = geometry(DIM,n).vertex(0);
  const int& v1 = geometry(DIM,n).vertex(1);
  const int& v2 = geometry(DIM,n).vertex(2);
  const int& v3 = geometry(DIM,n).vertex(3);

  const Point<DIM>& x0 = point(v0);
  const Point<DIM>& x1 = point(v1); 
  const Point<DIM>& x2 = point(v2);
  const Point<DIM>& x3 = point(v3);

  const double& m00 = move_direction[v0][0]; 
  const double& m01 = move_direction[v0][1]; 
  const double& m02 = move_direction[v0][2];

  const double& m10 = move_direction[v1][0]; 
  const double& m11 = move_direction[v1][1]; 
  const double& m12 = move_direction[v1][2];

  const double& m20 = move_direction[v2][0]; 
  const double& m21 = move_direction[v2][1]; 
  const double& m22 = move_direction[v2][2];

  const double& m30 = move_direction[v3][0]; 
  const double& m31 = move_direction[v3][1]; 
  const double& m32 = move_direction[v3][2];

  double volume = ((x1[0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
                   (x1[1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
                   (x1[2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
                   (x1[0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
                   (x1[1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
                   (x1[2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0]));
  double lambda[4];
  std::vector<std::vector<double> > v(p.size(), std::vector<double>(DIM));
  for (u_int j = 0; j < p.size(); j ++) {
    lambda[0] = ((x1[0] - p[j][0])*(x2[1] - p[j][1])*(x3[2] - p[j][2]) +
                 (x1[1] - p[j][1])*(x2[2] - p[j][2])*(x3[0] - p[j][0]) +
                 (x1[2] - p[j][2])*(x2[0] - p[j][0])*(x3[1] - p[j][1]) -
                 (x1[0] - p[j][0])*(x2[2] - p[j][2])*(x3[1] - p[j][1]) -
                 (x1[1] - p[j][1])*(x2[0] - p[j][0])*(x3[2] - p[j][2]) -
                 (x1[2] - p[j][2])*(x2[1] - p[j][1])*(x3[0] - p[j][0]))/volume;
    lambda[1] = ((p[j][0] - x0[0])*(x2[1] - x0[1])*(x3[2] - x0[2]) +
                 (p[j][1] - x0[1])*(x2[2] - x0[2])*(x3[0] - x0[0]) +
                 (p[j][2] - x0[2])*(x2[0] - x0[0])*(x3[1] - x0[1]) -
                 (p[j][0] - x0[0])*(x2[2] - x0[2])*(x3[1] - x0[1]) -
                 (p[j][1] - x0[1])*(x2[0] - x0[0])*(x3[2] - x0[2]) -
                 (p[j][2] - x0[2])*(x2[1] - x0[1])*(x3[0] - x0[0]))/volume;
    lambda[2] = ((x1[0] - x0[0])*(p[j][1] - x0[1])*(x3[2] - x0[2]) +
                 (x1[1] - x0[1])*(p[j][2] - x0[2])*(x3[0] - x0[0]) +
                 (x1[2] - x0[2])*(p[j][0] - x0[0])*(x3[1] - x0[1]) -
                 (x1[0] - x0[0])*(p[j][2] - x0[2])*(x3[1] - x0[1]) -
                 (x1[1] - x0[1])*(p[j][0] - x0[0])*(x3[2] - x0[2]) -
                 (x1[2] - x0[2])*(p[j][1] - x0[1])*(x3[0] - x0[0]))/volume;
    lambda[3] = ((x1[0] - x0[0])*(x2[1] - x0[1])*(p[j][2] - x0[2]) +
                 (x1[1] - x0[1])*(x2[2] - x0[2])*(p[j][0] - x0[0]) +
                 (x1[2] - x0[2])*(x2[0] - x0[0])*(p[j][1] - x0[1]) -
                 (x1[0] - x0[0])*(x2[2] - x0[2])*(p[j][1] - x0[1]) -
                 (x1[1] - x0[1])*(x2[0] - x0[0])*(p[j][2] - x0[2]) -
                 (x1[2] - x0[2])*(x2[1] - x0[1])*(p[j][0] - x0[0]))/volume;
    v[j][0] = lambda[0]*m00 + lambda[1]*m10 + lambda[2]*m20 + lambda[3]*m30;
    v[j][1] = lambda[0]*m01 + lambda[1]*m11 + lambda[2]*m21 + lambda[3]*m31;
    v[j][2] = lambda[0]*m02 + lambda[1]*m12 + lambda[2]*m22 + lambda[3]*m32;
  }
  return v;
}



double MovingMesh3D::moveDirectionDivergence(const u_int& n) const
{
  //   const int& v0 = geometry(2,n).vertex(0);
  //   const int& v1 = geometry(2,n).vertex(1);
  //   const int& v2 = geometry(2,n).vertex(2);
  //   const Point<2>& x0 = point(v0);
  //   const Point<2>& x1 = point(v1); 
  //   const Point<2>& x2 = point(v2);
  //   const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  //   const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  //   const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  //   return ((m10 - m00)*(x2[1] - x0[1]) - (x1[1] - x0[1])*(m20 - m00)
  // 	  + (x1[0] - x0[0])*(m21 - m01) - (m11 - m01)*(x2[0] - x0[0]))
  //     / ((x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]));
  return 0.0;
}


////////////////////////////////////////////////////////////////////////////////

void MovingMesh3D::Solver::Project(const sparse_matrix_t& M, 
                                   const std::vector<bound_t>& bm,
                                   sparse_matrix_t *& P,
                                   sparse_matrix_t *& PMPT,
                                   sparse_matrix_t *& Pt,
                                   std::vector<bound_t> *& pbm)
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
        /// 只有高材料标识的点能对低材料标识的点进行聚类
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
  pbm = new std::vector<bound_t>(n_core_point); /**< 为投影过的自由度的边界标号分配空间 */
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
      if (bm[k]%bm[l] != 0) continue; /// 边界点不向内部点做聚类
      P->add(i, l, 1.0/counter[l]);
      Pt->add(l, i, 1.0/counter[l]);
    }
  }

  PMPT = getPMPT(*P, M, *Pt);
}



sparse_matrix_t * MovingMesh3D::Solver::getPMPT(const sparse_matrix_t& P, 
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



void MovingMesh3D::Solver::lazyPMPT(const sparse_matrix_t& P, 
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



void MovingMesh3D::Solver::GaussSidel(const sparse_matrix_t& M, 
				      std::vector<vector_t>& x, 
				      const std::vector<vector_t>& r,
				      const std::vector<bound_t>& bm,
				      const u_int& s) const
{
  const SparsityPattern& spM = M.get_sparsity_pattern();
  const std::size_t * rowstart = spM.get_rowstart_indices();
  const u_int * colnums = spM.get_column_numbers();
  for (unsigned i = 0;i < s;i ++) {
    for (u_int j = 0;j < M.m();j ++) {
      double rj[3] = {r[0](j), r[1](j), r[2](j)};
      for (u_int k = rowstart[j] + 1;k < rowstart[j+1];k ++) {
	const double& a = M.global_entry(k); /**< M 的元素 */
	rj[0] -= a*x[0](colnums[k]);
	rj[1] -= a*x[1](colnums[k]);
	rj[2] -= a*x[2](colnums[k]);
      }
      const double& a = M.global_entry(rowstart[j]); /**< 对角元 */
      if (bm[j] == 1) { /**< 内部的点就不投影了 */
	x[0](j) = rj[0]/a; x[1](j) = rj[1]/a; x[2](j) = rj[2]/a;
	continue;
      }
      rj[0] = rj[0]/a - x[0](j); /**< 首先计算增量 */
      rj[1] = rj[1]/a - x[1](j);
      rj[2] = rj[2]/a - x[2](j);

      /// 这个投影过程有问题的，Todo: 已经修正
      u_int n_constraint = 0;
      std::vector<u_int> constraint_idx(domain->n_surf);
      for (u_int k = 0;k < domain->n_surf;k ++) { /// 计算约束的个数
	const Surface& surface = domain->surface[k];
	if (bm[j]%surface.boundary_mark != 0) continue; /**< 不在这个表面上 */
        constraint_idx[n_constraint ++] = k;
      }
      if (n_constraint == 1) { /// 在表面上
        const Surface& surface = domain->surface[constraint_idx[0]];
        const double * n = surface.logic_normal;
        double ip = rj[0]*n[0] + rj[1]*n[1] + rj[2]*n[2];
        rj[0] -= ip*n[0]; rj[1] -= ip*n[1]; rj[2] -= ip*n[2]; 
      }
      else if (n_constraint == 2) { /// 在棱上
        const Surface& surface0 = domain->surface[constraint_idx[0]];
        const double * n0 = surface0.logic_normal;
        const Surface& surface1 = domain->surface[constraint_idx[1]];
        const double * n1 = surface1.logic_normal;
        double n[3] = {n0[1]*n1[2] - n0[2]*n1[1],
                       n0[2]*n1[0] - n0[0]*n1[2],
                       n0[0]*n1[1] - n0[1]*n1[0]};
        double ip = (rj[0]*n[0] + rj[1]*n[1] + rj[2]*n[2])/
                    (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        rj[0] = ip*n[0]; rj[1] = ip*n[1]; rj[2] = ip*n[2]; 
      }
      else { /// 在顶点上
        rj[0] = 0.0; rj[1] = 0.0; rj[2] = 0.0;
      }
      x[0](j) += rj[0]; x[1](j) += rj[1]; x[2](j) += rj[2]; 
    }
  }
}


void MovingMesh3D::Solver::Init(sparse_matrix_t& M,
				const std::vector<bound_t>& bm,
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
    if (k < min_order) break;

    sparse_matrix_t *PMPT, *P, *Pt;
    std::vector<bound_t> *pbm;
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
#ifdef DEBUG
    std::cerr << "\tlevel " << n_project - 1
	      << " with order " << m
	      << std::endl;
#endif // DEGUG
  };
  is_initialized = true;
  std::cerr << " OK! grid levels: " << n_project << std::endl;
}



/**
 * to solve a linear system with the same sparsity pattern as already initialized
 * matrix
 * 
 */
void MovingMesh3D::Solver::reinit(sparse_matrix_t& M,
				  const std::vector<bound_t>& bm,
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



void MovingMesh3D::Solver::solve(std::vector<vector_t>& x, 
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
    projected_r[j] = new std::vector<vector_t>(DIM,
                                               vector_t(projected_matrix[j]->m()));
    projected_x[j] = new std::vector<vector_t>(DIM,
                                               vector_t(projected_matrix[j]->m()));
  }
  double init_residual = 0, residual;
  vector_t v0(projected_matrix[0]->m());
  for (j = 0;j < 3;j ++) {
    projected_matrix[0]->vmult(v0, x[j]);
    v0.add(-1, r[j]);
    for (k = 0;k < x[j].size();k ++) { /**< 清除边界能量 */
      v0(k) *= ((*boundary_mark[0])[k] == 1);
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
      for (j = 0;j < 3;j ++) {
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
      for (j = 0;j < 3;j ++) {
	project_matrix[i]->Tvmult(v0, (*projected_x[i+1])[j]);
	(*projected_x[i])[j] += v0;
      }
      GaussSidel(*projected_matrix[i], *projected_x[i], *projected_r[i], 
		 *boundary_mark[i], smooth_step);
    }
	
    residual = 0;
    v0.reinit(projected_matrix[0]->m());
    for (j = 0;j < 3;j ++) {
      projected_matrix[0]->vmult(v0, x[j]);
      v0.add(-1, r[j]);
      for (k = 0;k < x[j].size();k ++) { /**< 清除边界能量 */
	v0(k) *= ((*boundary_mark[0])[k] == 1);
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



void MovingMesh3D::Solver::clear()
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



MovingMesh3D::Solver::Solver() : 
  is_initialized(false),
  n_project(0),
  smooth_step(3)
{}



MovingMesh3D::Solver::~Solver()
{
  clear ();
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
