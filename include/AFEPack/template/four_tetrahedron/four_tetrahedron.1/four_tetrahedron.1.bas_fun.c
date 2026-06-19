/*************************************************************************
* 	four_tetrahedron.1.bas_fun.c : by R.Lie
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
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = get_volume(p, v[1], v[2], v[3]);
	val[0] /= volume;
}

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[6], v[5]);
	val[0] = get_volume(v[0], p, v[6], v[5]);
	val[0] /= volume;
	if (val[0] < 0.) val[0] = 0.;
}

void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[2], v[4], v[6]);
	val[0] = get_volume(v[0], p, v[4], v[6]);
	val[0] /= volume;
	if (val[0] < 0.) val[0] = 0.;
}

void lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[3], v[5], v[4]);
	val[0] = get_volume(v[0], p, v[5], v[4]);
	val[0] /= volume;
	if (val[0] < 0.) val[0] = 0.;
}

void lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = 0.25 * get_volume(v[0], v[1], v[2], v[3]);
	double d45, d56, d64;
	d56 = get_volume(v[0], p, v[5], v[6]);
	if (d56 <= 0.) {
		val[0] = 0.;
		return;
	}
	d45 = get_volume(v[0], p, v[4], v[5]);
	if (d45 <= 0.) {
		val[0] = get_volume(v[0], p, v[3], v[5]);
		val[0] /= volume;
		return;
	}
	d64 = get_volume(v[0], p, v[6], v[4]);
	if (d64 <= 0.) {
		val[0] = get_volume(v[0], p, v[6], v[2]);
		val[0] /= volume;
		return;
	}
	val[0] = d56/volume;
}

void lambda_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = 0.25 * get_volume(v[0], v[1], v[2], v[3]);
	double d45, d56, d64;
	d64 = get_volume(v[0], p, v[6], v[4]);
	if (d64 <= 0.) {
		val[0] = 0.;
		return;
	}
	d45 = get_volume(v[0], p, v[4], v[5]);
	if (d45 <= 0.) {
		val[0] = get_volume(v[0], p, v[4], v[3]);
		val[0] /= volume;
		return;
	}
	d56 = get_volume(v[0], p, v[5], v[6]);
	if (d56 <= 0.) {
		val[0] = get_volume(v[0], p, v[1], v[6]);
		val[0] /= volume;
		return;
	}
	val[0] = d64/volume;
}

void lambda_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = 0.25 * get_volume(v[0], v[1], v[2], v[3]);
	double d45, d56, d64;
	d45 = get_volume(v[0], p, v[4], v[5]);
	if (d45 <= 0.) {
		val[0] = 0.;
		return;
	}
	d56 = get_volume(v[0], p, v[5], v[6]);
	if (d56 <= 0.) {
		val[0] = get_volume(v[0], p, v[5], v[1]);
		val[0] /= volume;
		return;
	}
	d64 = get_volume(v[0], p, v[6], v[4]);
	if (d64 <= 0.) {
		val[0] = get_volume(v[0], p, v[2], v[4]);
		val[0] /= volume;
		return;
	}
	val[0] = d45/volume;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double volume = get_volume(v[0], v[1], v[2], v[3]);
	val[0] = co_det(v, 0, 0)/volume;
	val[1] = co_det(v, 0, 1)/volume;
	val[2] = co_det(v, 0, 2)/volume;
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[6], v[5]};
	double volume = get_volume(v[0], v[1], v[6], v[5]);
	double sign = get_volume(v[0], p, v[6], v[5])/volume;
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
	const double * v1[4] = {v[0], v[2], v[4], v[6]};
	double volume = get_volume(v[0], v[2], v[4], v[6]);
	double sign = get_volume(v[0], p, v[4], v[6])/volume;
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

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[3], v[5], v[4]};
	double volume = get_volume(v[0], v[3], v[5], v[4]);
	double sign = get_volume(v[0], p, v[5], v[4])/volume;
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

void gradient_lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[6], v[5]};
	const double * v2[4] = {v[0], v[2], v[4], v[6]};
	const double * v3[4] = {v[0], v[3], v[5], v[4]};
	const double * v4[4] = {v[0], v[4], v[5], v[6]};
	double volume = get_volume(v4[0], v4[1], v4[2], v4[3]);
	double d45 = get_volume(v[0], p, v[4], v[5])/volume;
	double d56 = get_volume(v[0], p, v[5], v[6])/volume;
	double d64 = get_volume(v[0], p, v[6], v[4])/volume;
	double count = 0.0;
	val[0] = 0.0; val[1] = 0.0; val[2] = 0.0;
	if (d45 <= ZERO) {
		val[0] += co_det(v3, 3, 0)/volume;
		val[1] += co_det(v3, 3, 1)/volume;
		val[2] += co_det(v3, 3, 2)/volume;
		count += 1.0;
	}
	if (d56 <= ZERO) {
		count += 1.0;
	}
	if (d64 <= ZERO) {
		val[0] += co_det(v2, 2, 0)/volume;
		val[1] += co_det(v2, 2, 1)/volume;
		val[2] += co_det(v2, 2, 2)/volume;
		count += 1.0;
	}
	if (d45 >= -ZERO && d56 >= -ZERO && d64 >= -ZERO) {
		val[0] += co_det(v4, 1, 0)/volume;
		val[1] += co_det(v4, 1, 1)/volume;
		val[2] += co_det(v4, 1, 2)/volume;
		count += 1.0;
	}
	val[0] /= count; val[1] /= count; val[2] /= count;
}



void gradient_lambda_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[6], v[5]};
	const double * v2[4] = {v[0], v[2], v[4], v[6]};
	const double * v3[4] = {v[0], v[3], v[5], v[4]};
	const double * v4[4] = {v[0], v[4], v[5], v[6]};
	double volume = get_volume(v4[0], v4[1], v4[2], v4[3]);
	double d45 = get_volume(v[0], p, v[4], v[5])/volume;
	double d56 = get_volume(v[0], p, v[5], v[6])/volume;
	double d64 = get_volume(v[0], p, v[6], v[4])/volume;
	double count = 0.0;
	val[0] = 0.0; val[1] = 0.0; val[2] = 0.0;
	if (d45 <= ZERO) {
		val[0] += co_det(v3, 2, 0)/volume;
		val[1] += co_det(v3, 2, 1)/volume;
		val[2] += co_det(v3, 2, 2)/volume;
		count += 1.0;
	}
	if (d56 <= ZERO) {
		val[0] += co_det(v1, 3, 0)/volume;
		val[1] += co_det(v1, 3, 1)/volume;
		val[2] += co_det(v1, 3, 2)/volume;
		count += 1.0;
	}
	if (d64 <= ZERO) {
		count += 1.0;
	}
	if (d45 >= -ZERO && d56 >= -ZERO && d64 >= -ZERO) {
		val[0] += co_det(v4, 2, 0)/volume;
		val[1] += co_det(v4, 2, 1)/volume;
		val[2] += co_det(v4, 2, 2)/volume;
		count += 1.0;
	}
	val[0] /= count; val[1] /= count; val[2] /= count;
}



void gradient_lambda_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	const double * v1[4] = {v[0], v[1], v[6], v[5]};
	const double * v2[4] = {v[0], v[2], v[4], v[6]};
	const double * v3[4] = {v[0], v[3], v[5], v[4]};
	const double * v4[4] = {v[0], v[4], v[5], v[6]};
	double volume = get_volume(v4[0], v4[1], v4[2], v4[3]);
	double d45 = get_volume(v[0], p, v[4], v[5])/volume;
	double d56 = get_volume(v[0], p, v[5], v[6])/volume;
	double d64 = get_volume(v[0], p, v[6], v[4])/volume;
	double count = 0.0;
	val[0] = 0.0; val[1] = 0.0; val[2] = 0.0;
	if (d45 <= ZERO) {
		count += 1.0;
	}
	if (d56 <= ZERO) {
		val[0] += co_det(v1, 2, 0)/volume;
		val[1] += co_det(v1, 2, 1)/volume;
		val[2] += co_det(v1, 2, 2)/volume;
		count += 1.0;
	}
	if (d64 <= ZERO) {
		val[0] += co_det(v2, 3, 0)/volume;
		val[1] += co_det(v2, 3, 1)/volume;
		val[2] += co_det(v2, 3, 2)/volume;
		count += 1.0;
	}
	if (d45 >= -ZERO && d56 >= -ZERO && d64 >= -ZERO) {
		val[0] += co_det(v4, 3, 0)/volume;
		val[1] += co_det(v4, 3, 1)/volume;
		val[2] += co_det(v4, 3, 2)/volume;
		count += 1.0;
	}
	val[0] /= count; val[1] /= count; val[2] /= count;
}


/*
*  end of file
**************************************************************************/
