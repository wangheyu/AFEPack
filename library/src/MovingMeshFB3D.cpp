/**
 * @file   MovingMeshFB3D.cpp
 * @author Yana Di, Ruo Li
 * @date   Mon May 28 22:28:00 2007
 * 
 * @brief  
 * 
 * 
 */

#include "MovingMeshFB3D.h"

AFEPACK_OPEN_NAMESPACE

#define DIM 3

typedef SparseMatrix<double> sparse_matrix_t;
typedef Vector<double> vector_t;

MovingMeshFB3D::MovingMeshFB3D()
{}



MovingMeshFB3D::~MovingMeshFB3D()
{}



void MovingMeshFB3D::readDomain(const std::string& filename)
{
  u_int i, j, k, l0, l1;

  readData(filename);
  index.resize(n_geometry(0));
  n_interior_node = 0;
  n_boundary_node = 0;
  for (i = 0;i < n_geometry(0);i ++) {
    Assert(geometry(0,i).vertex(0) == i, ExcInternalError());
    if (geometry(0, i).boundaryMark() != 0)
      index[i] = n_boundary_node ++;
    else
      index[i] = n_interior_node ++;
  }
  interior_node_index.resize(n_interior_node);
  boundary_node_index.resize(n_boundary_node);
  for (i = 0,j = 0, k = 0;i < n_geometry(0);i ++) {
    if (geometry(0, i).boundaryMark() != 0)
      boundary_node_index[j ++] = i;
    else
      interior_node_index[k ++] = i;
  }

  logical_node.resize(n_geometry(0));
  move_direction.resize(n_geometry(0));
  logical_move_direction.resize(n_geometry(0));
  monitor().resize(n_geometry(DIM));

  std::vector<u_int> n_coupling_node(n_geometry(0), 1);
  for (i = 0;i < n_geometry(1);i ++) {
    n_coupling_node[geometry(1,i).vertex(0)] ++;
    n_coupling_node[geometry(1,i).vertex(1)] ++;
  }
  spM.reinit(n_interior_node, n_interior_node, n_coupling_node);
  spN.reinit(n_interior_node, n_boundary_node, n_coupling_node);
  for (i = 0;i < n_geometry(3);i ++) {
    for (j = 0;j < 3;j ++) {
      l0 = geometry(2,i).vertex(j);
      if (geometry(0,l0).boundaryMark() == 0) {
	for (k = 0;k < 3;k ++) {
	  l1 = geometry(2,i).vertex(k);
	  if (geometry(0,l1).boundaryMark())
	    spN.add(index[l0], index[l1]);
	  else
	    spM.add(index[l0], index[l1]);
	}
      }
    }
  }
  spM.compress();
  spN.compress();

  getLogicalMesh();
}



void MovingMeshFB3D::moveMesh()
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



void MovingMeshFB3D::outputPhysicalMesh(const std::string& filename)
{
  writeData(filename);
}



void MovingMeshFB3D::outputLogicalMesh(const std::string& filename)
{
  Mesh<DIM,DIM> mesh = Mesh<DIM,DIM>(*this);
  for (u_int i = 0;i < n_geometry(0);i ++)
    mesh.point(i) = logical_node[i];
  mesh.writeData(filename);
}



void MovingMeshFB3D::getMoveDirection()
{
  u_int i, j, k, l, l0, l1;
  getMonitor();
  M.reinit(spM);
  N.reinit(spN);
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
      if (boundaryMark(0,l0) == 0) {
        for (k = 0;k < 4;k ++) {
          double d = monitor(i)*(omega[0][j]*omega[0][k] +
                                 omega[1][j]*omega[1][k] + 
                                 omega[2][j]*omega[2][k])/volume;
          l1 = geometry(DIM,i).vertex(k);
          if (boundaryMark(0,l1)) 
	    N.add(index[l0], index[l1], d);
	  else
	    M.add(index[l0], index[l1], d);
        }
      }
    }
  }
  solver.lazyReinit(M);

  Vector<double> b(n_boundary_node);
  Vector<double> r(n_interior_node);
  Vector<double> x(n_interior_node);
  for (k = 0;k < DIM;++ k) {
    for (i = 0;i < n_boundary_node;i ++) {
      b(i) = -logical_node[boundary_node_index[i]][k];
      logical_move_direction[j][k] = 0;
    }
    N.vmult(r, b);
    for (i = 0;i < n_interior_node;i ++)
      x(i) = logical_node[interior_node_index[i]][k];
    solver.solve(x, r);
    for (i = 0;i < n_interior_node;i ++) {
      logical_move_direction[interior_node_index[i]][k] = logical_node[interior_node_index[i]][k] - x(i);
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

  for (i = 0;i < n_boundary_node;i ++) {
    j = boundary_node_index[i];
    for (k = 0;k < DIM;++ k) {
      move_direction[j][k] = 0.0;
    }
  }
}



void MovingMeshFB3D::getMoveStepLength()
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



void MovingMeshFB3D::readDummy(std::ifstream& is)
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



void MovingMeshFB3D::getMonitor()
{
  std::fill(monitor().begin(), monitor().end(), 1);
}



void MovingMeshFB3D::smoothMonitor(u_int s)
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



void MovingMeshFB3D::updateMesh()
{
  for (u_int i = 0;i < n_geometry(0);i ++) {
    for (u_int k = 0;k < DIM;++ k) {
      point(i)[k] += move_step_length * move_direction[i][k];
    }
  }
}



void MovingMeshFB3D::getLogicalMesh()
{
  std::cout << "Computing logical mesh ..." << std::endl;
  u_int i, j, k, l0, l1;
  bmark_t bm;
  for (i = 0;i < n_geometry(0);i ++) {
    logical_node[i][0] = point(i)[0];
    logical_node[i][1] = point(i)[1];
    logical_node[i][2] = point(i)[2];
  }
}



std::vector<double> 
MovingMeshFB3D::moveDirection(const Point<DIM>& p, 
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
MovingMeshFB3D::moveDirection(const std::vector<Point<DIM> >& p, 
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



double MovingMeshFB3D::moveDirectionDivergence(const u_int& n) const
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

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
