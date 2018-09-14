///////////////////////////////////////////////////////////////////////////////
// burgers.cpp : by R.Lie
//

#include "burgers.h"


double ArgFunction::value(const double * p) const
{
	return 1./(1. + exp((p[0] + p[1] - t)/(2.*a)));
}

std::vector<double> ArgFunction::gradient(const double * p) const
{
	std::vector<double> v(2);
	return v;
};

///////////////////////////////////////////////////////////////////////////////

void burgers::Matrix::getElementMatrix(
		const Element<double,2>& element0,
		const Element<double,2>& element1,
		const ActiveElementPairIterator<2>::State state)
{
	int n_element_dof0 = elementDof0().size();
	int n_element_dof1 = elementDof1().size();
	double volume = element0.templateElement().volume();
	const QuadratureInfo<2>& quad_info = element0.findQuadratureInfo(algebricAccuracy());
	std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<AFEPack::Point<2> > q_point = element0.local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<double> > basis_value = element0.basis_function_value(q_point);
	std::vector<std::vector<std::vector<double> > > basis_gradient = element0.basis_function_gradient(q_point);
	for (int l = 0;l < n_quadrature_point;l ++) {
		double Jxw = quad_info.weight(l)*jacobian[l]*volume;
		for (int j = 0;j < n_element_dof0;j ++) {
			for (int k = 0;k < n_element_dof1;k ++) {
				elementMatrix(j,k) += Jxw*((1/dt)*basis_value[j][l]*basis_value[k][l]
						       + a*innerProduct(basis_gradient[j][l], basis_gradient[k][l]));
			}
		}
	}
};



burgers::burgers(const std::string& file) :
	mesh_file(file), a(0.005), t(0.25)
{};



burgers::~burgers()
{};



void burgers::initialize()
{
	readDomain(mesh_file);
	template_geometry.readData("triangle.tmp_geo");
	coord_transform.readData("triangle.crd_trs");
	template_dof.reinit(template_geometry);
	template_dof.readData("triangle.1.tmp_dof");
	basis_function.reinit(template_dof);
	basis_function.readData("triangle.1.bas_fun");
	template_element.resize(1);
	template_element[0].reinit(template_geometry,
			template_dof,
			coord_transform,
			basis_function);
	fem_space.reinit(*this, template_element);
	int n_element = n_geometry(2);
	fem_space.element().resize(n_element);
	for (int i = 0;i < n_element;i ++)
		fem_space.element(i).reinit(fem_space, i, 0);
	fem_space.buildElement();
	fem_space.buildDof();
	fem_space.buildDofBoundaryMark();
	u_h.reinit(fem_space);
	initialMesh();
};



void burgers::run()
{
	initialize();

	do {
		moveMesh();
		stepForward();
		outputSolution();
		std::cout << "t  = " << t << std::endl;
	} while (t < 0.4);
};



void burgers::initialMesh()
{
	std::cout << "Initialize mesh ... " << std::endl;
	double scale, scale_step = 0.2;
	scale = scale_step;
	do {
		initialValue();
		u_h.scale(scale);
		moveMesh();
		outputSolution();
		std::cout << "\r\tscale = " << scale << std::endl;
		scale += scale_step;
	} while (scale <= 1.0);
};



void burgers::initialValue()
{
	ArgFunction u(a, t);
	Operator::L2Project(u, u_h, Operator::LOCAL_LEAST_SQUARE, 3);
};



void burgers::boundaryValue()
{
};



void burgers::stepForward()
{
	int i, j, k, l;
	dt = 2.0e-3;
	FEMFunction<double,2> _u_h(u_h);
	for (i = 3;i > 0;i --) {
		Matrix matrix(fem_space, dt/i, a);
		matrix.algebricAccuracy() = 2;
		matrix.build();
		Vector<double> rhs(fem_space.n_dof());
		FEMSpace<double,2>::ElementIterator the_element = fem_space.beginElement();
		FEMSpace<double,2>::ElementIterator end_element = fem_space.endElement();
		for (;the_element != end_element;++ the_element) {
			double volume = the_element->templateElement().volume();
			const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(2);
			std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
			int n_quadrature_point = quad_info.n_quadraturePoint();
			std::vector<AFEPack::Point<2> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
			std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
			std::vector<double> _u_h_value = _u_h.value(q_point, *the_element);
			std::vector<double> u_h_value = u_h.value(q_point, *the_element);
			std::vector<std::vector<double> > u_h_gradient = u_h.gradient(q_point, *the_element);
			int n_element_dof = the_element->n_dof();
			const std::vector<int>& element_dof = the_element->dof();
			for (l = 0;l < n_quadrature_point;l ++) {
				double Jxw = quad_info.weight(l)*jacobian[l]*volume;
				for (j = 0;j < n_element_dof;j ++) {
					rhs(element_dof[j]) += Jxw*(u_h_value[l]*basis_value[j][l]/(dt/i)
							- u_h_value[l]*(u_h_gradient[l][0] + u_h_gradient[l][1])*basis_value[j][l]);
				}
			}
		}
		AMGSolver solver(matrix);
		solver.solve(u_h, rhs);
	};
	t += dt;
};



void burgers::getMonitor()
{
	int i, l;
	FEMSpace<double,2>::ElementIterator the_element = fem_space.beginElement();
	FEMSpace<double,2>::ElementIterator end_element = fem_space.endElement();
	for (i = 0;the_element != end_element;++ the_element) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(1);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<AFEPack::Point<2> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
		std::vector<std::vector<double> > u_h_gradient = u_h.gradient(q_point, *the_element);
		float d = 0, area = 0;
		for (l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			area += Jxw;
			d += Jxw*innerProduct(u_h_gradient[l], u_h_gradient[l]);
		}
		monitor(i ++) = d/area;
	}
	std::cout << "max monitor=" << *std::max_element(monitor().begin(), monitor().end())
		<< "\tmin monitor=" << *std::min_element(monitor().begin(), monitor().end())
		<< std::endl;
	smoothMonitor(2);
	for (i = 0;i < n_geometry(2);i ++)
		monitor(i) = 1./sqrt(1. + monitor(i));
};



void burgers::updateSolution()
{
	int i, j, l;
	FEMFunction<double,2> _u_h(u_h);
	const double& msl = moveStepLength();
	MassMatrix<2,double> matrix(fem_space);
	matrix.algebricAccuracy() = 2;
	matrix.build();
	AMGSolver solver(matrix);
	for (i = 1;i > 0;i --) {
		Vector<double> rhs(fem_space.n_dof());
		FEMSpace<double,2>::ElementIterator the_element = fem_space.beginElement();
		FEMSpace<double,2>::ElementIterator end_element = fem_space.endElement();
		for (;the_element != end_element;++ the_element) {
			double volume = the_element->templateElement().volume();
			const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(2);
			std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
			int n_quadrature_point = quad_info.n_quadraturePoint();
			std::vector<AFEPack::Point<2> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
			std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
			std::vector<double> _u_h_value = _u_h.value(q_point, *the_element);
			std::vector<std::vector<double> > u_h_gradient = u_h.gradient(q_point, *the_element);
			std::vector<std::vector<double> > move_vector = moveDirection(q_point, the_element->index());
			int n_element_dof = the_element->n_dof();
			const std::vector<int>& element_dof = the_element->dof();
			for (l = 0;l < n_quadrature_point;l ++) {
				double Jxw = quad_info.weight(l)*jacobian[l]*volume;
				for (j = 0;j < n_element_dof;j ++) {
					rhs(element_dof[j]) += Jxw*basis_value[j][l]*(_u_h_value[l]
							+ (1./i)*msl*innerProduct(move_vector[l], u_h_gradient[l]));
				}
			}
		}
		solver.solve(u_h, rhs);
	};
};



void burgers::outputSolution()
{
	outputPhysicalMesh("E");
	u_h.writeOpenDXData("u_h.dx");
};



//
// end of file
///////////////////////////////////////////////////////////////////////////////

