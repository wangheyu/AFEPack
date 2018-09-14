///////////////////////////////////////////////////////////////////////////
// four_tetrahedron.crd_trs.cpp : by R.Lie
//

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
		- (v1[2] - v0[2])*(v2[1] - v0[1])*(v3[0] - v0[0]))/6.;
}

void local_to_global(const double * lp, 
		const double ** lv, 
		const double ** gv, 
		double * gp)
{
	double lambda[4];
	double volume = get_volume(lv[0], lv[1], lv[2], lv[3]);
	lambda[0] = get_volume(lp, lv[1], lv[2], lv[3])/volume;
	lambda[1] = get_volume(lv[0], lp, lv[2], lv[3])/volume;
	lambda[2] = get_volume(lv[0], lv[1], lp, lv[3])/volume;
	lambda[3] = get_volume(lv[0], lv[1], lv[2], lp)/volume;
	gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[1][0] + lambda[2]*gv[2][0] + lambda[3]*gv[3][0];
	gp[1] = lambda[0]*gv[0][1] + lambda[1]*gv[1][1] + lambda[2]*gv[2][1] + lambda[3]*gv[3][1];
	gp[2] = lambda[0]*gv[0][2] + lambda[1]*gv[1][2] + lambda[2]*gv[2][2] + lambda[3]*gv[3][2];
};

void global_to_local(const double * gp, 
		const double ** lv, 
		const double ** gv, 
		double * lp)
{
	double lambda[4];
	double volume = get_volume(gv[0], gv[1], gv[2], gv[3]);
	lambda[0] = get_volume(gp, gv[1], gv[2], gv[3])/volume;
	lambda[1] = get_volume(gv[0], gp, gv[2], gv[3])/volume;
	lambda[2] = get_volume(gv[0], gv[1], gp, gv[3])/volume;
	lambda[3] = get_volume(gv[0], gv[1], gv[2], gp)/volume;
	lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[1][0] + lambda[2]*lv[2][0] + lambda[3]*lv[3][0];
	lp[1] = lambda[0]*lv[0][1] + lambda[1]*lv[1][1] + lambda[2]*lv[2][1] + lambda[3]*lv[3][1];
	lp[2] = lambda[0]*lv[0][2] + lambda[1]*lv[1][2] + lambda[2]*lv[2][2] + lambda[3]*lv[3][2];
};

double local_to_global_jacobian(const double * lp, 
		const double ** lv, 
		const double ** gv)
{
	double lvolume = get_volume(lv[0], lv[1], lv[2], lv[3]);
	double gvolume = get_volume(gv[0], gv[1], gv[2], gv[3]);
	return gvolume/lvolume;
};

double global_to_local_jacobian(const double * gp, 
		const double ** lv,
		const double ** gv)
{
	double lvolume = get_volume(lv[0], lv[1], lv[2], lv[3]);
	double gvolume = get_volume(gv[0], gv[1], gv[2], gv[3]);
	return lvolume/gvolume;
};

//
// end of file
///////////////////////////////////////////////////////////////////////////
