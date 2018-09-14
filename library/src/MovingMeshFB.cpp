/**
 * @file   MovingMeshFB.cpp
 * @author Robert Lie
 * @date   Mon May 28 22:14:49 2007
 * 
 * @brief  
 * 
 * 
 */

#include "MovingMeshFB.h"

AFEPACK_OPEN_NAMESPACE

MovingMeshFB::MovingMeshFB()
{
}



MovingMeshFB::~MovingMeshFB()
{
}



void MovingMeshFB::readDomain(const std::string& filename)
{
  int i, j, k, l, l0, l1;
  double d;

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
  monitor().resize(n_geometry(2));

  std::vector<int> n_coupling_node(n_geometry(0), 0);
  for (i = 0;i < n_geometry(1);i ++) {
    n_coupling_node[geometry(1,i).vertex(0)] ++;
    n_coupling_node[geometry(1,i).vertex(1)] ++;
  }
  int n_max_coupling_node = *std::max_element(n_coupling_node.begin(), n_coupling_node.end());
  n_max_coupling_node ++;
  spM.reinit(n_interior_node, n_interior_node, n_max_coupling_node);
  spN.reinit(n_interior_node, n_boundary_node, n_max_coupling_node);
  for (i = 0;i < n_geometry(2);i ++) {
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



void MovingMeshFB::moveMesh()
{
  int i, j, k, m;
  double epsilon = 0.2;
  double error = 1.;
  double area, h, l[3], d[3];
  while (error > epsilon) {
    getMoveDirection();
    error = 0;
    for (i = 0;i < n_geometry(2);i ++) {
      const Point<2>& x0 = point(geometry(2,i).vertex(0));
      const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
      const Point<2>& x2 = point(geometry(2,i).vertex(2));
      l[0] = (x2[0] - x1[0])*(x2[0] - x1[0]) + (x2[1] - x1[1])*(x2[1] - x1[1]);
      l[1] = (x0[0] - x2[0])*(x0[0] - x2[0]) + (x0[1] - x2[1])*(x0[1] - x2[1]);
      l[2] = (x1[0] - x0[0])*(x1[0] - x0[0]) + (x1[1] - x0[1])*(x1[1] - x0[1]);
      area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
      h = 0.5*area/sqrt(*std::max_element(&l[0], &l[3]));
      for (k = 0;k < 3;++ k) {
	m = geometry(2,i).vertex(k);
	for (d[k] = 0.0, j = 0;j < 2;j ++) {
	  d[k] += move_direction[m][j]*move_direction[m][j];
	}
      }
      h = sqrt(*std::max_element(&d[0], &d[3]))/h;
      if (error < h) error = h;
    }
    std::cerr << "mesh moving error = " << error << std::endl;
    getMoveStepLength();
    for (i = 0;i < n_move_step;i ++) {
      updateSolution();
      updateMesh();
    };
  };
}



void MovingMeshFB::outputPhysicalMesh(const std::string& filename)
{
  writeData(filename);
}



void MovingMeshFB::outputLogicalMesh(const std::string& filename)
{
  Mesh<2,2> mesh = Mesh<2,2>(*this);
  for (int i = 0;i < n_geometry(0);i ++)
    mesh.point(i) = logical_node[i];
  mesh.writeData(filename);
}



void MovingMeshFB::getMoveDirection()
{
  int i, j, k, l, l0, l1;
  getMonitor();
  M.reinit(spM);
  N.reinit(spN);
  for (i = 0;i < n_geometry(2);i ++) {
    const Point<2>& x0 = point(geometry(2,i).vertex(0));
    const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
    const Point<2>& x2 = point(geometry(2,i).vertex(2));
    double omega[2][3];
    double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
    omega[0][0] = x2[0] - x1[0]; omega[0][1] = x0[0] - x2[0];
    omega[0][2] = x1[0] - x0[0]; omega[1][0] = x2[1] - x1[1];
    omega[1][1] = x0[1] - x2[1]; omega[1][2] = x1[1] - x0[1];
    for (j = 0;j < 3;j ++) {
      l0 = geometry(2,i).vertex(j);
      if (boundaryMark(0,l0) == 0) {
	for (k = 0;k < 3;k ++) {
	  double d = monitor(i)*(omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k])/area;
	  l1 = geometry(2,i).vertex(k);
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
  for (k = 0;k < 2;++ k) {
    for (i = 0;i < n_boundary_node;i ++) {
      j = boundary_node_index[i];
      b(i) = -logical_node[j][k];
      logical_move_direction[j][k] = 0.0;
    }
    N.vmult(r, b);
    for (i = 0;i < n_interior_node;i ++)
      x(i) = logical_node[interior_node_index[i]][k];
    solver.solve(x, r);
    for (i = 0;i < n_interior_node;i ++) {
      j = interior_node_index[i];
      logical_move_direction[j][k] = logical_node[j][k] - x(i);
      move_direction[j] = 0.0;
    }
  }

  std::vector<double> mass_lumping(n_geometry(0), 0.);
  for (i = 0;i < n_geometry(2);i ++) {
    const int& v0 = geometry(2,i).vertex(0);
    const int& v1 = geometry(2,i).vertex(1);
    const int& v2 = geometry(2,i).vertex(2);
    const Point<2>& x0 = point(v0);
    const Point<2>& x1 = point(v1); 
    const Point<2>& x2 = point(v2);
    double jacobi[2][2];
    double area = ((logical_node[v1][0] - logical_node[v0][0])*(logical_node[v2][1] - logical_node[v0][1]) -
                   (logical_node[v1][1] - logical_node[v0][1])*(logical_node[v2][0] - logical_node[v0][0]));
    jacobi[0][0] = ((x1[0] - x0[0])*(logical_node[v2][1] - logical_node[v0][1]) -
                    (logical_node[v1][1] - logical_node[v0][1])*(x2[0] - x0[0]));
    jacobi[1][0] = ((x1[1] - x0[1])*(logical_node[v2][1] - logical_node[v0][1]) -
                    (logical_node[v1][1] - logical_node[v0][1])*(x2[1] - x0[1]));
    jacobi[0][1] = ((logical_node[v1][0] - logical_node[v0][0])*(x2[0] - x0[0]) -
                    (x1[0] - x0[0])*(logical_node[v2][0] - logical_node[v0][0]));
    jacobi[1][1] = ((logical_node[v1][0] - logical_node[v0][0])*(x2[1] - x0[1]) -
                    (x1[1] - x0[1]) * (logical_node[v2][0] - logical_node[v0][0]));
    for (j = 0;j < 3;j ++) {
      k = geometry(2,i).vertex(j);
      move_direction[k][0] += jacobi[0][0]*logical_move_direction[k][0] + jacobi[0][1]*logical_move_direction[k][1];
      move_direction[k][1] += jacobi[1][0]*logical_move_direction[k][0] + jacobi[1][1]*logical_move_direction[k][1];
      mass_lumping[k] += area;
    }
  }
  for (i = 0;i < n_geometry(0);i ++) {
    move_direction[i][0] /= mass_lumping[i];
    move_direction[i][1] /= mass_lumping[i];
  }

  for (i = 0;i < n_boundary_node;i ++) {
    j = boundary_node_index[i];
    move_direction[j][0] = 0.0;
    move_direction[j][1] = 0.0;
  }
}



void MovingMeshFB::getMoveStepLength()
{
  int i, j;
  double a, b, c;
  n_move_step = 1;
  move_step_length = 1.;
  for (i = 0;i < n_geometry(2);i ++) {
    const Point<2>& x0 = point(geometry(2,i).vertex(0));
    const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
    const Point<2>& x2 = point(geometry(2,i).vertex(2));
    a = (move_direction[geometry(2,i).vertex(1)][0] - move_direction[geometry(2,i).vertex(0)][0])
      * (move_direction[geometry(2,i).vertex(2)][1] - move_direction[geometry(2,i).vertex(0)][1])
      - (move_direction[geometry(2,i).vertex(1)][1] - move_direction[geometry(2,i).vertex(0)][1])
      * (move_direction[geometry(2,i).vertex(2)][0] - move_direction[geometry(2,i).vertex(0)][0]);
    b = move_direction[geometry(2,i).vertex(0)][0] * (x1[1] - x2[1])
      - move_direction[geometry(2,i).vertex(0)][1] * (x1[0] - x2[0])
      + move_direction[geometry(2,i).vertex(1)][0] * (x2[1] - x0[1])
      - move_direction[geometry(2,i).vertex(1)][1] * (x2[0] - x0[0])
      + move_direction[geometry(2,i).vertex(2)][0] * (x0[1] - x1[1])
      - move_direction[geometry(2,i).vertex(2)][1] * (x0[0] - x1[0]);
    c = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
    if (fabs(a)/(fabs(b) + fabs(c)) < 1.0e-4) {
      if (fabs(b) < 1.0e-4*fabs(c))
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



void MovingMeshFB::readDummy(std::ifstream& is)
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



void MovingMeshFB::getMonitor()
{
  std::fill(monitor().begin(), monitor().end(), 1);
}



void MovingMeshFB::smoothMonitor(int s)
{
  int i, j, k;
  std::vector<float> area(n_geometry(2));
  std::vector<float> mass_lumping(n_geometry(0), 0);
  std::vector<float> monitor1(n_geometry(0));
  for (i = 0;i < n_geometry(2);i ++) {
    const Point<2>& x0 = point(geometry(2,i).vertex(0));
    const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
    const Point<2>& x2 = point(geometry(2,i).vertex(2));
    area[i] = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
    for (j = 0;j < 3;j ++) {
      mass_lumping[geometry(2,i).vertex(j)] += area[i];
    }
  }
  for (i = 0;i < s;i ++) {
    std::fill(monitor1.begin(), monitor1.end(), 0);
    for (j = 0;j < n_geometry(2);j ++) {
      for (k = 0;k < 3;k ++) {
	monitor1[geometry(2,j).vertex(k)] += monitor(j)*area[j];
      }
    }
    for (j = 0;j < n_geometry(0);j ++)
      monitor1[j] /= 3*mass_lumping[j];
    std::fill(monitor().begin(), monitor().end(), 0);
    for (j = 0;j < n_geometry(2);j ++) {
      for (k = 0;k < 3;k ++) {
	monitor(j) += monitor1[geometry(2,j).vertex(k)];
      }
    }
  }
}



void MovingMeshFB::updateMesh()
{
  for (int i = 0;i < n_geometry(0);i ++) {
    point(i)[0] += move_step_length * move_direction[i][0];
    point(i)[1] += move_step_length * move_direction[i][1];
  }
}



void MovingMeshFB::getLogicalMesh()
{
  std::cout << "Computing logical mesh ..." << std::endl;
  int i, j, k, l0, l1;

  for (i = 0;i < n_geometry(0);i ++) {
    logical_node[i][0] = point(i)[0];
    logical_node[i][1] = point(i)[1];
  }
}



std::vector<double> MovingMeshFB::moveDirection(const Point<2>& p, const int& n) const
{
  const int& v0 = geometry(2,n).vertex(0);
  const int& v1 = geometry(2,n).vertex(1);
  const int& v2 = geometry(2,n).vertex(2);
  const Point<2>& x0 = point(v0);
  const Point<2>& x1 = point(v1); 
  const Point<2>& x2 = point(v2);
  const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
  double lambda[3];
  lambda[0] = ((x1[0] - p[0])*(x2[1] - p[1]) - (x1[1] - p[1])*(x2[0] - p[0]))/area;
  lambda[1] = ((x2[0] - p[0])*(x0[1] - p[1]) - (x2[1] - p[1])*(x0[0] - p[0]))/area;
  lambda[2] = ((x0[0] - p[0])*(x1[1] - p[1]) - (x0[1] - p[1])*(x1[0] - p[0]))/area;
  std::vector<double> v(2);
  v[0] = lambda[0]*m00 + lambda[1]*m10 + lambda[2]*m20;
  v[1] = lambda[0]*m01 + lambda[1]*m11 + lambda[2]*m21;
  return v;
}




std::vector<std::vector<double> > MovingMeshFB::moveDirection(const std::vector<Point<2> >& p, const int& n) const
{
  const int& v0 = geometry(2,n).vertex(0);
  const int& v1 = geometry(2,n).vertex(1);
  const int& v2 = geometry(2,n).vertex(2);
  const Point<2>& x0 = point(v0);
  const Point<2>& x1 = point(v1); 
  const Point<2>& x2 = point(v2);
  const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]);
  std::vector<std::vector<double> > v(p.size(), std::vector<double>(2));
  double lambda[3];
  for (int j = 0;j < p.size();j ++) {
    lambda[0] = ((x1[0] - p[j][0])*(x2[1] - p[j][1]) - (x1[1] - p[j][1])*(x2[0] - p[j][0]))/area;
    lambda[1] = ((x2[0] - p[j][0])*(x0[1] - p[j][1]) - (x2[1] - p[j][1])*(x0[0] - p[j][0]))/area;
    lambda[2] = ((x0[0] - p[j][0])*(x1[1] - p[j][1]) - (x0[1] - p[j][1])*(x1[0] - p[j][0]))/area;
    v[j][0] = lambda[0]*m00 + lambda[1]*m10 + lambda[2]*m20;
    v[j][1] = lambda[0]*m01 + lambda[1]*m11 + lambda[2]*m21;
  }
  return v;
}



double MovingMeshFB::moveDirectionDivergence(const int& n) const
{
  const int& v0 = geometry(2,n).vertex(0);
  const int& v1 = geometry(2,n).vertex(1);
  const int& v2 = geometry(2,n).vertex(2);
  const Point<2>& x0 = point(v0);
  const Point<2>& x1 = point(v1); 
  const Point<2>& x2 = point(v2);
  const double& m00 = move_direction[v0][0]; const double& m01 = move_direction[v0][1];
  const double& m10 = move_direction[v1][0]; const double& m11 = move_direction[v1][1];
  const double& m20 = move_direction[v2][0]; const double& m21 = move_direction[v2][1];
  return ((m10 - m00)*(x2[1] - x0[1]) - (x1[1] - x0[1])*(m20 - m00)
	  + (x1[0] - x0[0])*(m21 - m01) - (m11 - m01)*(x2[0] - x0[0]))
    / ((x1[0] - x0[0])*(x2[1] - x0[1]) - (x1[1] - x0[1])*(x2[0] - x0[0]));
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

 

