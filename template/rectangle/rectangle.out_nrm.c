#include <assert.h>
#include <math.h>

void out_normal(const double * p, const double ** v, int i, double * n)
{
	double l;
	int i1 = (i+1)%4;
	assert(i >= 0 && i <= 3);
	l = sqrt((v[i1][0] - v[i][0])*(v[i1][0] - v[i][0])
		+ (v[i1][1] - v[i][1])*(v[i1][1] - v[i][1]));
	n[0] = (v[i1][1] - v[i][1])/l;
	n[1] = -(v[i1][0] - v[i][0])/l;
}

