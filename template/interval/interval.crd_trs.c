/**************************************************
* interval.crd_trs.c : by R.Lie
*/

void local_to_global(const double * lp, 
		const double ** lv, 
		const double ** gv, 
		double * gp)
{
	double lambda[2];
	lambda[0] = (lv[1][0] - lp[0])/(lv[1][0] - lv[0][0]);
	lambda[1] = (lp[0] - lv[0][0])/(lv[1][0] - lv[0][0]);
	gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[1][0];
}

void global_to_local(const double * gp, 
		const double ** lv, 
		const double ** gv, 
		double * lp)
{
	double lambda[2];
	lambda[0] = (gv[1][0] - gp[0])/(gv[1][0] - gv[0][0]);
	lambda[1] = (gp[0] - gv[0][0])/(gv[1][0] - gv[0][0]);
	lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[1][0];
}

double local_to_global_jacobian(const double * lp, 
		const double ** lv, 
		const double ** gv)
{
	return (gv[1][0] - gv[0][0])/(lv[1][0] - lv[0][0]);
}

double global_to_local_jacobian(const double * gp, 
		const double ** lv,
		const double ** gv)
{
	return (lv[1][0] - lv[0][0])/(gv[1][0] - gv[0][0]);
}

/*
* end of file
**************************************************/
