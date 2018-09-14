/*************************************************************************
* 	twin_tirangle.1.D.bas_fun.c : by R.Lie
*/

void lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double area = (v[1][0] - v[0][0])*(v[3][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[3][0] - v[0][0]);
	val[0] = (v[1][0] - p[0])*(v[3][1] - p[1])
		- (v[1][1] - p[1])*(v[3][0] - p[0]);
	val[0] /= area;
};

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double area = (v[1][0] - v[0][0])*(v[3][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[3][0] - v[0][0]);
	val[0] = (v[3][0] - p[0])*(v[0][1] - p[1])
		- (v[3][1] - p[1])*(v[0][0] - p[0]);
	val[0] /= area;
};

void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double area = (v[1][0] - v[0][0])*(v[3][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[3][0] - v[0][0]);
	val[0] = (v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]);
	val[0] /= area;
};

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double area = (v[1][0] - v[0][0])*(v[3][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[3][0] - v[0][0]);
	val[0] = (v[1][1] - v[3][1])/area;
	val[1] = (v[3][0] - v[1][0])/area;
};

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double area = (v[1][0] - v[0][0])*(v[3][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[3][0] - v[0][0]);
	val[0] = (v[3][1] - v[0][1])/area;
	val[1] = (v[0][0] - v[3][0])/area;
};

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double area = (v[1][0] - v[0][0])*(v[3][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[3][0] - v[0][0]);
	val[0] = (v[0][1] - v[1][1])/area;
	val[1] = (v[1][0] - v[0][0])/area;
};

/*
*  end of file
**************************************************************************/
