/*************************************************************************
* 	tirangle.dual.0.bas_fun.c : by R.Lie
*/

#define ZERO (1.0e-8)

#ifdef _cplusplus

extern "C" {
	void chi_0(const double * p, const double ** v, void * value);
	void chi_1(const double * p, const double ** v, void * value);
	void chi_2(const double * p, const double ** v, void * value);
	void gradient_chi_0(const double * p, const double ** v, void * value);
	void gradient_chi_1(const double * p, const double ** v, void * value);
	void gradient_chi_2(const double * p, const double ** v, void * value);
};

#endif

void chi_0(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double l0 = (v[1][0] - p[0])*(v[2][1] - p[1]) - (v[1][1] - p[1])*(v[2][0] - p[0]);
	double l1 = (v[2][0] - p[0])*(v[0][1] - p[1]) - (v[2][1] - p[1])*(v[0][0] - p[0]);
	double l2 = (v[0][0] - p[0])*(v[1][1] - p[1]) - (v[0][1] - p[1])*(v[1][0] - p[0]);
	if (l1 - l0 > ZERO || l2 - l0 > ZERO)
		val[0] = 0.0;
	else if (l1 - l0 >= -ZERO || l2 - l0 >= -ZERO)
		val[0] = 0.5;
	else
		val[0] = 1.0;
};

void chi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double l0 = (v[1][0] - p[0])*(v[2][1] - p[1]) - (v[1][1] - p[1])*(v[2][0] - p[0]);
	double l1 = (v[2][0] - p[0])*(v[0][1] - p[1]) - (v[2][1] - p[1])*(v[0][0] - p[0]);
	double l2 = (v[0][0] - p[0])*(v[1][1] - p[1]) - (v[0][1] - p[1])*(v[1][0] - p[0]);
	if (l0 - l1 > ZERO || l2 - l1 > ZERO)
		val[0] = 0.0;
	else if (l0 - l1 >= -ZERO || l2 - l1 >= -ZERO)
		val[0] = 0.5;
	else
		val[0] = 1.0;
};

void chi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double l0 = (v[1][0] - p[0])*(v[2][1] - p[1]) - (v[1][1] - p[1])*(v[2][0] - p[0]);
	double l1 = (v[2][0] - p[0])*(v[0][1] - p[1]) - (v[2][1] - p[1])*(v[0][0] - p[0]);
	double l2 = (v[0][0] - p[0])*(v[1][1] - p[1]) - (v[0][1] - p[1])*(v[1][0] - p[0]);
	if (l0 - l2 > ZERO || l1 - l2 > ZERO)
		val[0] = 0.0;
	else if (l0 - l2 >= -ZERO || l1 - l2 >= -ZERO)
		val[0] = 0.5;
	else
		val[0] = 1.0;
};

void gradient_chi_0(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 0.0;
};

void gradient_chi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 0.0;
};

void gradient_chi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	val[0] = 0.0;
	val[1] = 0.0;
};

/*
*  end of file
**************************************************************************/
