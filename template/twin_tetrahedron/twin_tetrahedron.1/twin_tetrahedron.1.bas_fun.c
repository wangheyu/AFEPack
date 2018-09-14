/*************************************************************************
* 	twin_tetrahedron.1.bas_fun.c : by R.Lie
*/

double get_volume(const double * v0,
		const double * v1,
		const double * v2,
		const double * v3)
{
	return ((v1[0] - v0[0])*(v2[1] - v0[1])*(v3[2] - v0[2])
		+ (v1[1] - v0[1])*(v2[2] - v0[2])*(v3[0] - v0[0])
		+ (v1[2] - v0[2])*(v2[0] - v0[0])*(v3[1] - v0[1])
		- (v1[0] - v0[0])*(v2[2] - v0[2])*(v3[1] - v0[1])
		- (v1[1] - v0[1])*(v2[0] - v0[0])*(v3[2] - v0[2])
		- (v1[2] - v0[2])*(v2[1] - v0[1])*(v3[0] - v0[0]));
}


#define ZERO 1.0e-12
#define co_det(v, m, n) (\
		((m%2==0)?-1.:1.) * (\
		  (v[(m+2)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
		* (v[(m+3)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
		- (v[(m+2)%4][(n+2)%3] - v[(m+1)%4][(n+2)%3]) \
		* (v[(m+3)%4][(n+1)%3] - v[(m+1)%4][(n+1)%3]) \
	) \
	)

void lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[3], v[4]);
	val[0] = get_volume(p, v[1], v[3], v[4]);
	val[0] /= volume;
}

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[4]);
	val[0] = get_volume(v[0], p, v[2], v[4]);
	val[0] /= volume;
	if (val[0] < 0.) val[0] = 0.;
}

void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double sign = get_volume(v[0], p, v[2], v[4]);
	double volume = 0.5 * get_volume(v[0], v[1], v[3], v[4]);
	if (sign > 0.) {
		val[0] = get_volume(v[0], v[1], p, v[4]);
		val[0] /= volume;
	}
	else {
		val[0] = get_volume(v[0], p, v[3], v[4]);
		val[0] /= volume;
	}
}

void lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[2], v[3], v[4]);
	val[0] = get_volume(v[0], v[2], p, v[4]);
	val[0] /= volume;
	if (val[0] < 0.) val[0] = 0.;
}

void lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[3], v[4]);
	val[0] = get_volume(v[0], v[1], v[3], p);
	val[0] /= volume;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v0[4] = {v[0], v[1], v[3], v[4]};
	double volume = get_volume(v[0], v[1], v[3], v[4]);
	val[0] = co_det(v0, 0, 0)/volume;
	val[1] = co_det(v0, 0, 1)/volume;
	val[2] = co_det(v0, 0, 2)/volume;
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[2], v[4]};
	const double * v2[4] = {v[0], v[2], v[3], v[4]};
	double sign = get_volume(v[0], p, v[2], v[4]);
	double volume = get_volume(v[0], v[1], v[2], v[4]);
	if (sign > ZERO) {
		val[0] = co_det(v1, 1, 0)/volume;
		val[1] = co_det(v1, 1, 1)/volume;
		val[2] = co_det(v1, 1, 2)/volume;
	}
	else if (sign < -ZERO) {
		val[0] = 0.;
		val[1] = 0.;
		val[2] = 0.;
	}
	else {
		val[0] = 0.5*co_det(v1, 1, 0)/volume;
		val[1] = 0.5*co_det(v1, 1, 1)/volume;
		val[2] = 0.5*co_det(v1, 1, 2)/volume;
	}
}

void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[2], v[4]};
	const double * v2[4] = {v[0], v[2], v[3], v[4]};
	double sign = get_volume(v[0], p, v[2], v[4]);
	double volume = 0.5*get_volume(v[0], v[1], v[3], v[4]);
	if (sign > ZERO) {
		val[0] = co_det(v1, 2, 0)/volume;
		val[1] = co_det(v1, 2, 1)/volume;
		val[2] = co_det(v1, 2, 2)/volume;
	}
	else if (sign < -ZERO) {
		val[0] = co_det(v2, 1, 0)/volume;
		val[1] = co_det(v2, 1, 1)/volume;
		val[2] = co_det(v2, 1, 2)/volume;
	}
	else {
		val[0] = 0.5*(co_det(v1, 2, 0) + co_det(v2, 1, 0))/volume;
		val[1] = 0.5*(co_det(v1, 2, 1) + co_det(v2, 1, 1))/volume;
		val[2] = 0.5*(co_det(v1, 2, 2) + co_det(v2, 1, 2))/volume;
	}
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[2], v[4]};
	const double * v2[4] = {v[0], v[2], v[3], v[4]};
	double sign = get_volume(v[0], v[2], p, v[4]);
	double volume = get_volume(v[0], v[1], v[2], v[4]);
	if (sign > ZERO) {
		val[0] = co_det(v2, 2, 0)/volume;
		val[1] = co_det(v2, 2, 1)/volume;
		val[2] = co_det(v2, 2, 2)/volume;
	}
	else if (sign < -ZERO) {
		val[0] = 0.;
		val[1] = 0.;
		val[2] = 0.;
	}
	else {
		val[0] = 0.5*co_det(v2, 2, 0)/volume;
		val[1] = 0.5*co_det(v2, 2, 1)/volume;
		val[2] = 0.5*co_det(v2, 2, 2)/volume;
	}
}



void gradient_lambda_5(const double * p, const double ** v, void * value)
{
        double * val = (double *)value;
        const double * v0[4] = {v[0], v[1], v[3], v[4]};
        double volume = get_volume(v[0], v[1], v[3], v[4]);
        val[0] = co_det(v0, 3, 0)/volume;
        val[1] = co_det(v0, 3, 1)/volume;
        val[2] = co_det(v0, 3, 2)/volume;
}

/*
*  end of file
**************************************************************************/
