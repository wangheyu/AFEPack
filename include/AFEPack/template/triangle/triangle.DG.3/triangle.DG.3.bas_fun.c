/*************************************************************************
*  triangle.DG.3.bas_fun.c
*/

void phi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = lambda[0]*(3*lambda[0] - 1.0)*(3*lambda[0] - 2.0)/2.0;
}

void phi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = lambda[1]*(3*lambda[1] - 1.0)*(3*lambda[1] - 2.0)/2.0;
}

void phi_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = lambda[2]*(3*lambda[2] - 1.0)*(3*lambda[2] - 2.0)/2.0;
}

void phi_12(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 9.0*lambda[0]*lambda[1]*(3*lambda[0] - 1.0)/2.0;
}

void phi_21(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 9.0*lambda[0]*lambda[1]*(3*lambda[1] - 1.0)/2.0;
}

void phi_23(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 9.0*lambda[1]*lambda[2]*(3*lambda[1] - 1.0)/2.0;
}

void phi_32(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 9.0*lambda[1]*lambda[2]*(3*lambda[2] - 1.0)/2.0;
}

void phi_13(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 9.0*lambda[0]*lambda[2]*(3*lambda[0] - 1.0)/2.0;
}

void phi_31(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 9.0*lambda[0]*lambda[2]*(3*lambda[2] - 1.0)/2.0;
}

void phi_123(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[3];
	double area = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[2][1] - p[1])
		- (v[1][1] - p[1])*(v[2][0] - p[0]))/area;
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area;
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area;
	val[0] = 27.0*lambda[0]*lambda[1]*lambda[2];
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
}

void gradient_phi_3(const double * p, const double ** v, void * value)
{
}

void gradient_phi_12(const double * p, const double ** v, void * value)
{
}

void gradient_phi_21(const double * p, const double ** v, void * value)
{
}

void gradient_phi_23(const double * p, const double ** v, void * value)
{
}

void gradient_phi_32(const double * p, const double ** v, void * value)
{
}

void gradient_phi_13(const double * p, const double ** v, void * value)
{
}

void gradient_phi_31(const double * p, const double ** v, void * value)
{
}

void gradient_phi_123(const double * p, const double ** v, void * value)
{
}

/*
*  end of file
**************************************************************************/
