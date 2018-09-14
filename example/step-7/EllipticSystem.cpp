

// EllipticSystem.cpp

#include "EllipticSystem.h"

double _a0_00_(const double * p) {
	if (p[0]*p[0] > p[1]*p[1])
		return 10.0;
	else
		return 1.0;
};

double _a0_01_(const double * p) {
	return 0.0;
};

double _a0_10_(const double * p) {
        return 0.0;
};

double _a0_11_(const double * p) {
        if (p[0]*p[0] > p[1]*p[1])
                return 10.0;
        else
                return 1.0;
};

double _a1_00_(const double * p) {
        if (p[0]*p[0] +  p[1]*p[1] < 0.5)
                return 10.0;
        else
                return 1.0;
};

double _a1_01_(const double * p) {
        return 0.0;
};

double _a1_10_(const double * p) {
        return 0.0;
};

double _a1_11_(const double * p) {
        if (p[0]*p[0] + p[1]*p[1] < 0.5)
                return 10.0;
        else
                return 1.0;
};

double _c0_0_(const double * p) {
	return 1.0;
};

double _c0_1_(const double * p) {
        return 0.25;
};

double _c1_0_(const double * p) {
        return 0.25;
};

double _c1_1_(const double * p) {
        return 1.0;
};

double _u0_(const double * p) {
	return 0.0;
};

double _u1_(const double * p) {
        return 0.0;
};

double _f0_(const double * p) {
        return 1.0;
};

double _f1_(const double * p) {
        return 1.0;
};

void EllipticSystem::Matrix_A::getElementMatrix(
		const Element<double,DIM>& element0,
		const Element<double,DIM>& element1,
		const ActiveElementPairIterator<DIM>::State state)
{
	int j, k, l, j1, k1;
	int n_element_dof0 = elementDof0().size();
	int n_element_dof1 = elementDof1().size();
	double volume = element0.templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebricAccuracy());
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<double> > basis_value = element0.basis_function_value(q_point);
	std::vector<std::vector<std::vector<double> > > basis_gradient = element0.basis_function_gradient(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
		double Jxw = quad_info.weight(l)*jacobian[l]*volume;
		double a_value[DIM][DIM];
		for (j1 = 0;j1 < DIM;j1 ++)
			for (k1 = 0;k1 < DIM;k1 ++)
				a_value[j1][k1] = (*a[j1][k1])(q_point[l]);
		double c_value = (*c)(q_point[l]);
		for (j = 0;j < n_element_dof0;j ++) {
			for (k = 0;k < n_element_dof1;k ++) {
				for (j1 = 0;j1 < DIM;j1 ++) {
					for (k1 = 0;k1 < DIM;k1 ++) {
						elementMatrix(j, k) += Jxw * a_value[j1][k1] 
							* basis_gradient[j][l][j1] * basis_gradient[k][l][k1];
					}
				}
				elementMatrix(j, k) += Jxw * c_value * basis_value[j][l] * basis_value[k][l];
			}
		}
	}
};

// this matrix is a weighted L^2 inner product between two finite element spaces. this function is also a very
// standard example to calculate the matrix of a bilinear operator between two finite element spaces.
void EllipticSystem::Matrix_C::getElementMatrix(
		const Element<double,DIM>& element0,
		const Element<double,DIM>& element1,
		const ActiveElementPairIterator<DIM>::State state)
{
	int n_element_dof0 = elementDof0().size();
	int n_element_dof1 = elementDof1().size();
	if (state == ActiveElementPairIterator<DIM>::GREAT_THAN) {
		double volume = element1.templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = element1.findQuadratureInfo(algebricAccuracy());
		std::vector<double> jacobian = element1.local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<AFEPack::Point<DIM> > q_point = element1.local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > basis_value0 = element0.basis_function_value(q_point);
		std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			double c_value = (*c)(q_point[l]);
			for (int j = 0;j < n_element_dof0;j ++) {
				for (int k = 0;k < n_element_dof1;k ++) {
					elementMatrix(j,k) += Jxw * c_value * basis_value0[j][l] * basis_value1[k][l];
				}
			}
		}
	}
	else {
		double volume = element0.templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = element0.findQuadratureInfo(algebricAccuracy());
		std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<AFEPack::Point<DIM> > q_point = element0.local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > basis_value0 = element0.basis_function_value(q_point);
		std::vector<std::vector<double> > basis_value1 = element1.basis_function_value(q_point);
		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			double c_value = (*c)(q_point[l]);
			for (int j = 0;j < n_element_dof0;j ++) {
				for (int k = 0;k < n_element_dof1;k ++) {
					elementMatrix(j,k) += Jxw * c_value * basis_value0[j][l] * basis_value1[k][l];
				}
			}
		}
	}
};

EllipticSystem::EllipticSystem(const std::string& filename) :
	mesh_file(filename)
{
	init();
};

void EllipticSystem::run()
{
	do {
		buildFEMSpace();
		solve();
		saveData();
		getError();

		getIndicator();
		adaptMesh();
	} while (1);
};

void EllipticSystem::init()
{
	// read in the mesh data file to construct the hierarchy geometry tree. the mesh data is assumed to
	// be easymesh format
	h_geometry_tree.readEasyMesh(mesh_file);

	// construct the irregular meshes based on the hierarchy geometry tree
	for (int i = 0;i < DIM;i ++)
		irregular_mesh[i].reinit(h_geometry_tree);

	// construct the piecewise linear triangle element
	triangle_template_geometry.readData("triangle.tmp_geo");
	triangle_coord_transform.readData("triangle.crd_trs");
	triangle_1_template_dof.reinit(triangle_template_geometry);
	triangle_1_template_dof.readData("triangle.1.tmp_dof");
	triangle_1_basis_function.reinit(triangle_1_template_dof);
	triangle_1_basis_function.readData("triangle.1.bas_fun");

	// construct the piecewise linear twin triangle element. this element is special for this adaptive
	// package which is in fact the combination of two piecewise linear triangle element
	twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
	twin_triangle_coord_transform.readData("twin_triangle.crd_trs");
	twin_triangle_1_template_dof.reinit(twin_triangle_template_geometry);
	twin_triangle_1_template_dof.readData("twin_triangle.1.tmp_dof");
	twin_triangle_1_basis_function.reinit(twin_triangle_1_template_dof);
	twin_triangle_1_basis_function.readData("twin_triangle.1.bas_fun");

	// construct the two template elements using the information got by the last lines
	template_element.resize(2);
	template_element[0].reinit(triangle_template_geometry, 
		triangle_1_template_dof, 
		triangle_coord_transform, 
		triangle_1_basis_function);
	template_element[1].reinit(twin_triangle_template_geometry, 
		twin_triangle_1_template_dof, 
		twin_triangle_coord_transform, 
		twin_triangle_1_basis_function);
};

void EllipticSystem::buildFEMSpace()
{
	for (int n = 0;n < DIM;n ++) {
		irregular_mesh[n].semiregularize();
		irregular_mesh[n].regularize();
		RegularMesh<DIM>& regular_mesh = irregular_mesh[n].regularMesh();
		fem_space[n].reinit(regular_mesh, template_element);
		int n_element = regular_mesh.n_geometry(DIM);
		fem_space[n].element().resize(n_element);
		for (int i = 0;i < n_element;i ++) {
			const GeometryBM& geometry = regular_mesh.geometry(DIM, i);
			switch (geometry.n_vertex()) {
				case 3: // this is a triangle, we use the triangle template
					fem_space[n].element(i).reinit(fem_space[n],i,0);
					break;
				case 4: // this is a twin-triangle, we use the twin-triangle template
					fem_space[n].element(i).reinit(fem_space[n],i,1);
					break;
				default: // something must be wrong
					Assert(false, ExcNotImplemented());
			}
		}
		fem_space[n].buildElement();
		fem_space[n].buildDof();
		fem_space[n].buildDofBoundaryMark();
		solution[n].reinit(fem_space[n]);
	}
};

void EllipticSystem::adaptMesh()
{
	MeshAdaptor<DIM> mesh_adaptor;
	for (int i = 0;i < DIM;i ++) {
		mesh_adaptor.reinit(irregular_mesh[i]);
		mesh_adaptor.convergenceOrder() = 1.;
		mesh_adaptor.refineStep() = 1;
		mesh_adaptor.setIndicator(indicator[i]);
		mesh_adaptor.tolerence() = 1.0e-6;
		mesh_adaptor.adapt();
	}

};

// as can be seen easily, under this standard frame of an adaptive program, the solver for the original
// problem is the certral part.
void EllipticSystem::solve()
{
	// for the first equation
	// build the matrix A0
	Matrix_A A0(&_a0_00_, &_a0_01_, &_a0_10_, &_a0_11_, &_c0_0_, fem_space[0]);
	A0.algebricAccuracy() = 2;
	A0.build();

	// prepare the right hand side
	Vector<double> f0;
	Operator::L2Discretize(&_f0_, fem_space[0], f0, 3);

	// prepare the boundary condition
	BoundaryFunction<double,DIM> u0_boundary(BoundaryConditionInfo::DIRICHLET, 1, &_u0_);
	BoundaryConditionAdmin<double,DIM> u0_boundary_admin(fem_space[0]);
	u0_boundary_admin.add(u0_boundary);
	u0_boundary_admin.apply(A0, solution[0], f0);

	// prepare the linear solver for A0
	AMGSolver A0_solver(A0);

	// build the matrix C0
	Matrix_C C0(&_c0_1_, fem_space[0], fem_space[1]);
	C0.algebricAccuracy() = 2;
	C0.build();

	// for the second equation
	// build the matrix A1
	Matrix_A A1(&_a1_00_, &_a1_01_, &_a1_10_, &_a1_11_, &_c1_1_, fem_space[1]);
	A1.algebricAccuracy() = 2;
	A1.build();

	// prepare the right hand side
	Vector<double> f1;
	Operator::L2Discretize(&_f1_, fem_space[1], f1, 3);

	// prepare the boundary condition
	BoundaryFunction<double,DIM> u1_boundary(BoundaryConditionInfo::DIRICHLET, 1, &_u1_);
	BoundaryConditionAdmin<double,DIM> u1_boundary_admin(fem_space[1]);
	u1_boundary_admin.add(u1_boundary);
	u1_boundary_admin.apply(A1, solution[1], f1);

	// prepare the linear solver for A1
	AMGSolver A1_solver(A1);

	// build the matrix C1
	Matrix_C C1(&_c1_0_, fem_space[1], fem_space[0]);
	C1.algebricAccuracy() = 2;
	C1.build();

	// we adopt the block Gauss-Sidle iteration as the solver because $c^i$ are small
	int i;
	double error;
	double tolerence = 1.0e-6;
	do {
		// backup the old solution
		FEMFunction<double,DIM> old_solution[DIM];
		for (i = 0;i < DIM;i ++)
			old_solution[i] = solution[i];

		// update $u^0$
		Vector<double> C0xu1(solution[0]);
		C0.vmult(C0xu1, (Vector<double>)solution[1]);
		u0_boundary_admin.clearEntry(C0xu1);
		C0xu1.sadd(-1.0, f0);
		A0_solver.solve(solution[0], C0xu1);
	
		// update $u^1$
		Vector<double> C1xu0(solution[1]);
		C1.vmult(C1xu0, (Vector<double>)solution[0]);
		u1_boundary_admin.clearEntry(C1xu0);
		C1xu0.sadd(-1.0, f1);
		A1_solver.solve(solution[1], C1xu0);

		// calculate the error
		error = 0.0;
		for (i = 0;i < DIM;i ++) {
			old_solution[i].add(-1.0, solution[i]);
			error += pow(Functional::L2Norm(old_solution[i],3), 2.0);
		}
		error = sqrt(error);
		std::cout << "\r\terror = " << error << std::flush;
	} while (error > tolerence);
	std::cout << std::endl;
};

// we adopted only the leading term of the posteriori error estimator - saying the jump between the
// element boundary - to be calculated as the indicator. and the coefficient in the divergence operator
// is omitted, too, as this is only a very simple demonstration to show how to use the code.
void EllipticSystem::getIndicator()
{
	for (int n = 0;n < DIM;n ++) {
		RegularMesh<DIM>& regular_mesh = irregular_mesh[n].regularMesh();
		indicator[n].reinit(regular_mesh);
		int n_face = regular_mesh.n_geometry(DIM-1);
		std::vector<bool> flag(n_face, false);
		std::vector<double> jump(n_face);
		FEMSpace<double,DIM>::ElementIterator the_element = fem_space[n].beginElement();
		FEMSpace<double,DIM>::ElementIterator end_element = fem_space[n].endElement();
		for (;the_element != end_element;the_element ++) {
			const GeometryBM& geometry = the_element->geometry();
			for (int i = 0;i < geometry.n_boundary();i ++) {
				int j = geometry.boundary(i);
				const GeometryBM& side = regular_mesh.geometry(1, j);
				const AFEPack::Point<DIM>& p0 = regular_mesh.point(regular_mesh.geometry(0, side.boundary(0)).vertex(0));
				const AFEPack::Point<DIM>& p1 = regular_mesh.point(regular_mesh.geometry(0, side.boundary(1)).vertex(0));
				std::vector<double> gradient = solution[n].gradient(midpoint(p0, p1), *the_element);
				if (flag[j]) {
					jump[j] -= gradient[0]*(p0[1] - p1[1]) + gradient[1]*(p1[0] - p0[0]);
					flag[j] = false;
				}
				else {
					jump[j] = gradient[0]*(p0[1] - p1[1]) + gradient[1]*(p1[0] - p0[0]);
					flag[j] = true;
				}
			}
		}
		the_element = fem_space[n].beginElement();
		for (int m = 0;the_element != end_element;the_element ++, m ++) {
			const GeometryBM& geometry = the_element->geometry();
			indicator[n][m] = 0.0;
			for (int i = 0;i < geometry.n_boundary();i ++) {
				int j = geometry.boundary(i);
				if (flag[j]) continue; // this side is on the boundary of the domain
				indicator[n][m] += jump[j]*jump[j];
			}
		}
	}
};

void EllipticSystem::getError()
{
	double error;

	error = Functional::L2Error(solution[0], FunctionFunction<double>(&_u0_), 3);
	std::cout << "|| u_0 - u_0h ||_L^2 = " << error << std::endl;

	error = Functional::L2Error(solution[1], FunctionFunction<double>(&_u1_), 3);
	std::cout << "|| u_1 - u_1h ||_L^2 = " << error << std::endl;
};

void EllipticSystem::saveData()
{
	solution[0].writeOpenDXData("u_0h.dx");
	solution[1].writeOpenDXData("u_1h.dx");
//	solution[0].writeEasyMeshData("u_0h");
//	solution[1].writeEasyMeshData("u_1h");
};

// end of file

