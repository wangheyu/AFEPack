///////////////////////////////////////////////////////////////////////////////
// burgers.h : by R.Lie
//

#include <AFEPack/MovingMesh2D.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/BilinearOperator.h>


class ArgFunction : public Function<double>
{
private:
	double a;
	double t;
public:
	ArgFunction(const double& _a, const double& _t) :
		a(_a), t(_t)
	{};
	~ArgFunction() {};
public:
	virtual double value(const double *) const;
	virtual std::vector<double> gradient(const double *) const;
};


class burgers : public MovingMesh2D
{
public:
	class Matrix : public StiffMatrix<2,double>
	{
	private:
		double dt, a;
	public:
		Matrix(FEMSpace<double,2>& sp, const double& _dt, const double& _a) :
			dt(_dt), a(_a),
			StiffMatrix<2,double>(sp) {};
		virtual ~Matrix() {};
	public:
		virtual void getElementMatrix(const Element<double,2>& e0,
				const Element<double,2>& e1,
				const ActiveElementPairIterator<2>::State state);
	};
private:
	TemplateGeometry<2> template_geometry;
	CoordTransform<2,2> coord_transform;
	TemplateDOF<2> template_dof;
	BasisFunctionAdmin<double,2,2> basis_function;
	std::vector<TemplateElement<double,2,2> > template_element;
	FEMSpace<double,2> fem_space;

	std::string mesh_file;
	double t;
	double dt;
	double a;
	FEMFunction<double,2> u_h;
public:
	burgers(const std::string& file);
	virtual ~burgers();
public:
	void run();
	void stepForward();
	void initialize();
	void initialValue();
	void boundaryValue();
	void initialMesh();
	virtual void getMonitor();
	virtual void updateSolution();
	virtual void outputSolution();
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////

