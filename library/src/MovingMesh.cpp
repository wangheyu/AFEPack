/**
 * @file   MovingMesh.cpp
 * @author Robert Lie
 * @date   Wed Jun 21 08:17:34 2003
 * 
 * @brief  
 * 
 * 
 */
#include "MovingMesh.h"

AFEPACK_OPEN_NAMESPACE

MovingMesh::MovingMesh()
{
}



MovingMesh::~MovingMesh()
{
}



void MovingMesh::readDomain(const std::string& filename)
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
    readDummy(is); is >> d;
    readDummy(is); is >> domain.physical_domain_vertex[i].boundary_mark;
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

  std::vector<int> n_coupling_node(n_geometry(0), 0);
  for (i = 0;i < n_geometry(1);i ++) {
    n_coupling_node[geometry(1,i).vertex(0)] ++;
    n_coupling_node[geometry(1,i).vertex(1)] ++;
  }
  int n_max_coupling_node = *std::max_element(n_coupling_node.begin(), n_coupling_node.end());
  n_max_coupling_node ++;
  spM.reinit(n_interior_node, n_interior_node, n_max_coupling_node);
  spN.reinit(n_interior_node, n_boundary_node, n_max_coupling_node);
  spMb.reinit(2*n_boundary_node + n_boundary_constraint, 
              2*n_boundary_node + n_boundary_constraint, 
              n_max_coupling_node + 2);
  for (i = 0;i < n_geometry(2);i ++) {
    for (j = 0;j < 3;j ++) {
      l0 = geometry(2,i).vertex(j);
      if (geometry(0,l0).boundaryMark()) {
	for (k = 0;k < 3;k ++) {
	  l1 = geometry(2,i).vertex(k);
	  if (geometry(0,l1).boundaryMark()) {
	    spMb.add(index[l0], index[l1]);
	    spMb.add(index[l0] + n_boundary_node, index[l1] + n_boundary_node);
	  }
	}
      }
      else {
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

  for (i = 0, k = 2*n_boundary_node;i < domain.n_edge;i ++) {
    l = mb_node[i].size();
    for (j = 0;j < l;j ++) {
      l0 = index[mb_node[i][j]];
      spMb.add(k, l0);
      spMb.add(l0, k);
      spMb.add(k, l0 + n_boundary_node);
      spMb.add(l0 + n_boundary_node, k);
      k ++;
    }
  }
  spMb.compress();

  getLogicalMesh();
}



void MovingMesh::moveMesh()
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



void MovingMesh::outputPhysicalMesh(const std::string& filename)
{
  writeData(filename);
}



void MovingMesh::outputLogicalMesh(const std::string& filename)
{
  Mesh<2,2> mesh = Mesh<2,2>(*this);
  for (int i = 0;i < n_geometry(0);i ++)
    mesh.point(i) = logical_node[i];
  mesh.writeData(filename);
}



void MovingMesh::getMoveDirection()
{
  int i, j, k, l, l0, l1;
  getMonitor();
  SparseMatrix<double> M(spM);
  SparseMatrix<double> N(spN);
  SparseMatrix<double> Mb(spMb);
  for (i = 0;i < n_geometry(2);i ++) {
    const Point<2>& x0 = point(geometry(2,i).vertex(0));
    const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
    const Point<2>& x2 = point(geometry(2,i).vertex(2));
    double omega[2][3];
    double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
    omega[0][0] = x2[0] - x1[0];
    omega[0][1] = x0[0] - x2[0];
    omega[0][2] = x1[0] - x0[0];
    omega[1][0] = x2[1] - x1[1];
    omega[1][1] = x0[1] - x2[1];
    omega[1][2] = x1[1] - x0[1];
    for (j = 0;j < 3;j ++) {
      l0 = geometry(2,i).vertex(j);
      if (boundaryMark(0,l0)) {
	for (k = 0;k < 3;k ++) {
	  double d = monitor(i)*(omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k])/area;
	  l1 = geometry(2,i).vertex(k);
	  if (boundaryMark(0,l1)) {
	    Mb.add(index[l0], index[l1], d);
	    Mb.add(index[l0] + n_boundary_node, index[l1] + n_boundary_node, d);
	  }
	}
      }
      else {
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
  Vector<double> b(n_boundary_node);
  Vector<double> r(n_interior_node);
  Vector<double> x(n_interior_node);
  Vector<double> xb(2*n_boundary_node + n_boundary_constraint);
  Vector<double> rb(2*n_boundary_node + n_boundary_constraint);
  Vector<double> rb1(n_boundary_node);
  for (i = 0;i < n_boundary_node;i ++)
    b(i) = -logical_node[boundary_node_index[i]][0];
  N.vmult(r, b);
  for (i = 0;i < n_interior_node;i ++)
    x(i) = logical_node[interior_node_index[i]][0];
  AMGSolver solver(M);
  solver.solve(x, r);
  for (i = 0;i < n_interior_node;i ++)
    logical_move_direction[interior_node_index[i]][0] = logical_node[interior_node_index[i]][0] - x(i);
  N.Tvmult(rb1, x);
  for (i = 0;i < n_boundary_node;i ++) {
    rb(i) = -rb1(i);
    xb(i) = logical_node[boundary_node_index[i]][0];
  }

  for (i = 0;i < n_boundary_node;i ++)
    b(i) = -logical_node[boundary_node_index[i]][1];
  N.vmult(r, b);
  for (i = 0;i < n_interior_node;i ++)
    x(i) = logical_node[interior_node_index[i]][1];
  solver.solve(x, r);
  for (i = 0;i < n_interior_node;i ++)
    logical_move_direction[interior_node_index[i]][1] = logical_node[interior_node_index[i]][1] - x(i);
  N.Tvmult(rb1, x);
  for (i = 0;i < n_boundary_node;i ++) {
    rb(i + n_boundary_node) = -rb1(i);
    xb(i + n_boundary_node) = logical_node[boundary_node_index[i]][1];
  }

  for (i = 0, k = 2*n_boundary_node;i < domain.n_edge;i ++) {
    l = mb_node[i].size();
    double w[2];
    w[0] = logical_node[mb_node[i][l-1]][0] - logical_node[mb_node[i][0]][0];
    w[1] = logical_node[mb_node[i][l-1]][1] - logical_node[mb_node[i][0]][1];
    double d = logical_node[mb_node[i][0]][0]*w[1] - logical_node[mb_node[i][0]][1]*w[0];
    for (j = 0;j < l;j ++) {
      l0 = index[mb_node[i][j]];
      Mb.add(k, l0, w[1]);
      Mb.add(l0, k, w[1]);
      Mb.add(k, l0 + n_boundary_node, -w[0]);
      Mb.add(l0 + n_boundary_node, k, -w[0]);
      rb(k ++) = d;
    }
  }
  SparseILU<double> ilu(spMb);
  ilu.decompose(Mb);
  SolverControl solver_control(200, 1e-6);
  PrimitiveVectorMemory <> vector_memory;
  SolverGMRES<> gmres(solver_control, vector_memory);
  gmres.solve (Mb, xb, rb,
	       PreconditionUseMatrix <SparseILU<double>, Vector<double> >
	       (ilu, &SparseILU<double>::apply_decomposition<double>));
  for (i = 0;i < n_boundary_node;i ++) {
    j = boundary_node_index[i];
    logical_move_direction[j][0] = logical_node[j][0] - xb(i);
    logical_move_direction[j][1] = logical_node[j][1] - xb(i + n_boundary_node);
  }

  for (i = 0;i < n_geometry(0);i ++) {
    move_direction[i][0] = 0;
    move_direction[i][1] = 0;
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
    double area = (logical_node[v1][0] - logical_node[v0][0])*(logical_node[v2][1] - logical_node[v0][1])
      - (logical_node[v1][1] - logical_node[v0][1])*(logical_node[v2][0] - logical_node[v0][0]);
    jacobi[0][0] = (x1[0] - x0[0])*(logical_node[v2][1] - logical_node[v0][1])
      - (logical_node[v1][1] - logical_node[v0][1])*(x2[0] - x0[0]);
    jacobi[1][0] = (x1[1] - x0[1])*(logical_node[v2][1] - logical_node[v0][1])
      - (logical_node[v1][1] - logical_node[v0][1])*(x2[1] - x0[1]);
    jacobi[0][1] = (logical_node[v1][0] - logical_node[v0][0])*(x2[0] - x0[0])
      - (x1[0] - x0[0])*(logical_node[v2][0] - logical_node[v0][0]);
    jacobi[1][1] = (logical_node[v1][0] - logical_node[v0][0])*(x2[1] - x0[1])
      - (x1[1] - x0[1]) * (logical_node[v2][0] - logical_node[v0][0]);
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

  for (i = 0;i < domain.n_edge;i ++) {
    k = mb_node[i].size();
    double w[2], w2;
    w[0] = point(mb_node[i][0])[0] - point(mb_node[i][k-1])[0];
    w[1] = point(mb_node[i][0])[1] - point(mb_node[i][k-1])[1];
    w2 = w[0]*w[0] + w[1]*w[1];
    move_direction[mb_node[i][0]][0] = 0;
    move_direction[mb_node[i][0]][1] = 0;
    for (j = 1;j < k-1;j ++) {
      l = mb_node[i][j];
      double d = (move_direction[l][0]*w[0] + move_direction[l][1]*w[1])/w2;
      move_direction[l][0] = d*w[0];
      move_direction[l][1] = d*w[1];
    }
    move_direction[mb_node[i][k-1]][0] = 0;
    move_direction[mb_node[i][k-1]][1] = 0;
  }
}



void MovingMesh::getMoveStepLength()
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



void MovingMesh::readDummy(std::ifstream& is)
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



void MovingMesh::getMonitor()
{
  std::fill(monitor().begin(), monitor().end(), 1);
}



void MovingMesh::smoothMonitor(int s)
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
    for (j = 0;j < 3;j ++)
      mass_lumping[geometry(2,i).vertex(j)] += area[i];
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



void MovingMesh::updateMesh()
{
  for (int i = 0;i < n_geometry(0);i ++) {
    point(i)[0] += move_step_length * move_direction[i][0];
    point(i)[1] += move_step_length * move_direction[i][1];
  }
}



void MovingMesh::getLogicalMesh()
{
  std::cout << "Computing logical mesh ..." << std::endl;
  int i, j, k, l0, l1;

  for (i = 0;i < n_geometry(0);i ++) {
    logical_node[i][0] = point(i)[0];
    logical_node[i][1] = point(i)[1];
  }
  return;

  bmark_t bm;
  for (i = 0;i < n_boundary_node;i ++) {
    j = boundary_node_index[i];
    bm = geometry(0,j).boundaryMark();
    for (k = 0;k < domain.n_edge;k ++) {
      if (domain.edge[k].boundary_mark == bm) {
	break;
      }
    }
    double p[2][2];
    for (l0 = 0;l0 < domain.n_vertex;l0 ++) {
      if(domain.physical_domain_vertex[l0].index == domain.edge[k].vertex[0]) {
	break;
      }
    }
    p[0][0] = domain.physical_domain_vertex[l0][0];
    p[0][1] = domain.physical_domain_vertex[l0][1];
    for (l1 = 0;l1 < domain.n_vertex;l1 ++) {
      if(domain.physical_domain_vertex[l1].index == domain.edge[k].vertex[1]) {
	break;
      }
    }
    p[1][0] = domain.physical_domain_vertex[l1][0];
    p[1][1] = domain.physical_domain_vertex[l1][1];
    if (fabs(domain.physical_domain_vertex[l0][0] - domain.physical_domain_vertex[l1][0]) > 
	fabs(domain.physical_domain_vertex[l0][1] - domain.physical_domain_vertex[l1][1]))
      k = 0;
    else
      k = 1;
    double lambda = (point(j)[k] - domain.physical_domain_vertex[l0][k]) 
      / (domain.physical_domain_vertex[l1][k] - domain.physical_domain_vertex[l0][k]);
    logical_node[j][0] = (1-lambda)*domain.logical_domain_vertex[l0][0] + lambda*domain.logical_domain_vertex[l1][0];
    logical_node[j][1] = (1-lambda)*domain.logical_domain_vertex[l0][1] + lambda*domain.logical_domain_vertex[l1][1];
  }
  SparseMatrix<double> M(spM);
  SparseMatrix<double> N(spN);
  for (i = 0;i < n_geometry(2);i ++) {
    const Point<2>& x0 = point(geometry(2,i).vertex(0));
    const Point<2>& x1 = point(geometry(2,i).vertex(1)); 
    const Point<2>& x2 = point(geometry(2,i).vertex(2));
    double omega[2][3];
    double area = (x1[0] - x0[0])*(x2[1] - x0[1]) - (x2[0] - x0[0])*(x1[1] - x0[1]);
    omega[0][0] = x2[0] - x1[0];
    omega[0][1] = x0[0] - x2[0];
    omega[0][2] = x1[0] - x0[0];
    omega[1][0] = x2[1] - x1[1];
    omega[1][1] = x0[1] - x2[1];
    omega[1][2] = x1[1] - x0[1];
    for (j = 0;j < 3;j ++) {
      l0 = geometry(2,i).vertex(j);
      if (geometry(0,l0).boundaryMark() == 0) {
	for (k = 0;k < 3;k ++) {
	  double d = (omega[0][j]*omega[0][k] + omega[1][j]*omega[1][k])/area;
	  l1 = geometry(2,i).vertex(k);
	  if (geometry(0,l1).boundaryMark())
	    N.add(index[l0], index[l1], d);
	  else
	    M.add(index[l0], index[l1], d);
	}
      }
    }
  }
  Vector<double> b(n_boundary_node);
  Vector<double> r(n_interior_node);
  Vector<double> x(n_interior_node);
  for (i = 0;i < n_boundary_node;i ++) {
    b(i) = -logical_node[boundary_node_index[i]][0];
  }
  N.vmult(r, b);
  AMGSolver solver(M);
  solver.solve(x, r);
  for (i = 0;i < n_interior_node;i ++) {
    logical_node[interior_node_index[i]][0] = x(i);
  }

  for (i = 0;i < n_boundary_node;i ++) {
    b(i) = -logical_node[boundary_node_index[i]][1];
  }
  N.vmult(r, b);
  solver.solve(x, r);
  for (i = 0;i < n_interior_node;i ++) {
    logical_node[interior_node_index[i]][1] = x(i);
  }
}



void MovingMesh::parseBoundary()
{
  std::cout << "Parsing boundary nodes and faces ..." << std::endl;
  return;

  int i, j, k, m, n;
  bmark_t bm;
  mb_node.resize(domain.n_edge);
  mb_face.resize(domain.n_edge);
  std::vector<int> count(domain.n_edge, 0);
  face_index.resize(n_geometry(1));
  for (i = 0, j = 0, n_boundary_face = 0;i < n_geometry(1);i ++) {
    if (geometry(1,i).boundaryMark() != 0)
      face_index[i] = n_boundary_face ++;
    else
      face_index[i] = j ++;
  }
  boundary_face.resize(n_boundary_face);
  for (i = 0, k = 0;i < n_geometry(1);i ++)
    if (geometry(1,i).boundaryMark() != 0)
      boundary_face[k ++] = i;
  for (i = 0;i < n_boundary_face;i ++) {
    bm = geometry(1, boundary_face[i]).boundaryMark();
    for (j = 0;j < domain.n_edge;j ++) {
      if (bm == domain.edge[j].boundary_mark) {
	count[j] ++;
	break;
      }
    }
  }
  for (i = 0, n_boundary_constraint = 0;i < domain.n_edge;i ++) {
    mb_node[i].resize(count[i] + 1);
    mb_face[i].resize(count[i]);
    n_boundary_constraint += count[i] + 1;
  }
  std::vector<bool> flag(n_boundary_face, false);
  for (i = 0;i < domain.n_edge;i ++) {
    for (j = 0;j < n_boundary_face;j ++) {
      bm = geometry(1, boundary_face[j]).boundaryMark();
      if (bm == domain.edge[i].boundary_mark) {
	m = boundary_face[j];
	mb_node[i][0] = geometry(1,m).vertex(0);
	mb_node[i][1] = geometry(1,m).vertex(1);
	mb_face[i][0] = m;
	flag[j] = true;
	break;
      }
    }
    m = mb_node[i][1];
    k = 1;
    for (j = 0;j < n_boundary_face;) {
      if (flag[j]) {
	j ++;
	continue;
      }
      bm = geometry(1, boundary_face[j]).boundaryMark();
      if (bm != domain.edge[i].boundary_mark) {
	j ++;
	continue;
      }
      n = boundary_face[j];
      if (geometry(1, n).vertex(0) == m) {
	m = mb_node[i][k+1] = geometry(1, n).vertex(1);
	mb_face[i][k ++] = n;
	flag[j] = true;
	j = 0;
      }
      else if (geometry(1, n).vertex(1) == m) {
	m = mb_node[i][k+1] = geometry(1, n).vertex(0);
	mb_face[i][k ++] = n;
	flag[j] = true;
	j = 0;
      }
      else {
	j ++;
	continue;
      }
    }
    if (k == count[i]) continue;
    m = mb_node[i][0];
    std::reverse(&mb_node[i][0], &mb_node[i][k+1]);
    std::reverse(&mb_face[i][0], &mb_face[i][k]);
    for (j = 0;j < n_boundary_face;) {
      if (flag[j]) {
	j ++;
	continue;
      }
      bm = geometry(1, boundary_face[j]).boundaryMark();
      if (bm != domain.edge[i].boundary_mark) {
	j ++;
	continue;
      }
      n = boundary_face[j];
      if (geometry(1, n).vertex(0) == m) {
	m = mb_node[i][k+1] = geometry(1, n).vertex(1);
	mb_face[i][k ++] = n;
	flag[j] = true;
	j = 0;
      }
      else if (geometry(1, n).vertex(1) == m) {
	m = mb_node[i][k+1] = geometry(1, n).vertex(0);
	mb_face[i][k ++] = n;
	flag[j] = true;
	j = 0;
      }
      else {
	j ++;
	continue;
      }
    }
    Assert(k == count[i], ExcInternalError());
  }
}



std::vector<double> MovingMesh::moveDirection(const Point<2>& p, const int& n) const
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




std::vector<std::vector<double> > MovingMesh::moveDirection(const std::vector<Point<2> >& p, const int& n) const
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



double MovingMesh::moveDirectionDivergence(const int& n) const
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

 
