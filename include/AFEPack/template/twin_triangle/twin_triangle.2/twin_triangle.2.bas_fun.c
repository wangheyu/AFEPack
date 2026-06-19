/*************************************************************************
*  twin_triangle.2.bas_fun.c
*/

#include "math.h"

#define dlambda0_dx ((v[1][1] - v[2][1])/area)
#define dlambda0_dy ((v[2][0] - v[1][0])/area)
#define dlambda1_dx ((v[2][1] - v[0][1])/area)
#define dlambda1_dy ((v[0][0] - v[2][0])/area)
#define dlambda20_dx ((v[0][1] - v[1][1])/area)
#define dlambda20_dy ((v[1][0] - v[0][0])/area)
#define dlambda21_dx ((v[3][1] - v[0][1])/area)
#define dlambda21_dy ((v[0][0] - v[3][0])/area)
#define dlambda3_dx ((v[0][1] - v[2][1])/area)
#define dlambda3_dy ((v[2][0] - v[0][0])/area)


void get_lambda(const double * p, const double ** v, double * lambda, double * area)
{
	area[0] = (v[1][0] - v[0][0])*(v[2][1] - v[0][1])
		- (v[1][1] - v[0][1])*(v[2][0] - v[0][0]);
	lambda[0] = ((v[1][0] - p[0])*(v[3][1] - p[1])
		- (v[1][1] - p[1])*(v[3][0] - p[0]))/(2.*area[0]);
	lambda[1] = ((v[2][0] - p[0])*(v[0][1] - p[1])
		- (v[2][1] - p[1])*(v[0][0] - p[0]))/area[0];
	lambda[2] = ((v[0][0] - p[0])*(v[1][1] - p[1])
		- (v[0][1] - p[1])*(v[1][0] - p[0]))/area[0];
	lambda[3] = ((v[3][0] - p[0])*(v[0][1] - p[1])
		- (v[3][1] - p[1])*(v[0][0] - p[0]))/area[0];
}

void phi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = lambda[0]*(2*lambda[0] - 1.0);
}

void phi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] < 0) {
		val[0] = 0.;
	}
	else {
		val[0] = lambda[1]*(2*lambda[1] - 1.0);
	}
}

void phi_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] < 0) {
		val[0] = lambda[3]*(2*lambda[3] - 1.0);
	}
	else {
		val[0] = lambda[2]*(2*lambda[2] - 1.0);
	}
}

void phi_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] < 0) {
		lambda[1] = -lambda[1];
		val[0] = lambda[1]*(2*lambda[1] - 1.0);
	}
	else {
		val[0] = 0.;
	}
}

void phi_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] < 0) {
		val[0] = 0.;
	}
	else {
		val[0] = 4.0*lambda[0]*lambda[1];
	}
}

void phi_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] < 0) {
		val[0] = 0.;
	}
	else {
		val[0] = 4.0*lambda[1]*lambda[2];
	}
}

void phi_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] > 0) {
		val[0] = 0.;
	}
	else {
		val[0] = -4.0*lambda[1]*lambda[3];
	}
}

void phi_8(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] > 0) {
		val[0] = 0.;
	}
	else {
		val[0] = -4.0*lambda[0]*lambda[1];
	}
}

void phi_9(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (lambda[1] > 0) {
		val[0] = 4.0*lambda[0]*lambda[2];
	}
	else {
		val[0] = 4.0*lambda[0]*lambda[3];
	}
}

void gradient_phi_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = (4.0*lambda[0] - 1.0)*dlambda0_dx;
	val[1] = (4.0*lambda[0] - 1.0)*dlambda0_dy;
}

void gradient_phi_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = (4.0*lambda[1] - 1.0)*dlambda1_dx;
	val[1] = (4.0*lambda[1] - 1.0)*dlambda1_dy;
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] /= 2.0;
		val[1] /= 2.0;
	}
	else if (lambda[1] < 0) {
		val[0] = 0.0;
		val[1] = 0.0;
	}
}

void gradient_phi_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] = 0.5*((4.0*lambda[2] - 1.0)*dlambda20_dx + (4.0*lambda[3] - 1.0)*dlambda21_dx);
		val[1] = 0.5*((4.0*lambda[2] - 1.0)*dlambda20_dy + (4.0*lambda[3] - 1.0)*dlambda21_dy);
	}
	else if (lambda[1] > 0) {
		val[0] = (4.0*lambda[2] - 1.0)*dlambda20_dx;
		val[1] = (4.0*lambda[2] - 1.0)*dlambda20_dy;
	}
	else {
		val[0] = (4.0*lambda[3] - 1.0)*dlambda21_dx;
		val[1] = (4.0*lambda[3] - 1.0)*dlambda21_dy;
	}
}

void gradient_phi_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = (4.0*lambda[1] + 1.0)*dlambda1_dx;
	val[1] = (4.0*lambda[1] + 1.0)*dlambda1_dy;
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] /= 2.0;
		val[1] /= 2.0;
	}
	else if (lambda[1] > 0) {
		val[0] = 0.0;
		val[1] = 0.0;
	}
}

void gradient_phi_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = 4.0*(lambda[1]*dlambda0_dx + lambda[0]*dlambda1_dx);
	val[1] = 4.0*(lambda[1]*dlambda0_dy + lambda[0]*dlambda1_dy);
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] /= 2.0;
		val[1] /= 2.0;
	}
	else if (lambda[1] < 0) {
		val[0] = 0.0;
		val[1] = 0.0;
	}
}

void gradient_phi_6(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = 4.0*(lambda[2]*dlambda1_dx + lambda[1]*dlambda20_dx);
	val[1] = 4.0*(lambda[2]*dlambda1_dy + lambda[1]*dlambda20_dy);
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] /= 2.0;
		val[1] /= 2.0;
	}
	else if (lambda[1] < 0) {
		val[0] = 0.0;
		val[1] = 0.0;
	}
}

void gradient_phi_7(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = -4.0*(lambda[1]*dlambda21_dx + lambda[3]*dlambda1_dx);
	val[1] = -4.0*(lambda[1]*dlambda21_dy + lambda[3]*dlambda1_dy);
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] /= 2.0;
		val[1] /= 2.0;
	}
	else if (lambda[1] > 0) {
		val[0] = 0.0;
		val[1] = 0.0;
	}
}

void gradient_phi_8(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	val[0] = -4.0*(lambda[0]*dlambda1_dx + lambda[1]*dlambda0_dx);
	val[1] = -4.0*(lambda[0]*dlambda1_dy + lambda[1]*dlambda0_dy);
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] /= 2.0;
		val[1] /= 2.0;
	}
	else if (lambda[1] > 0) {
		val[0] = 0.0;
		val[1] = 0.0;
	}
}

void gradient_phi_9(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
	double lambda[4];
	double area;
	get_lambda(p, v, lambda, &area);
	if (fabs(lambda[1]) < 1.0e-4*area) {
		val[0] = 2.0*(lambda[0]*dlambda20_dx + lambda[2]*dlambda0_dx)
			+ 2.0*(lambda[0]*dlambda21_dx + lambda[3]*dlambda0_dx);
		val[1] = 2.0*(lambda[0]*dlambda20_dy + lambda[2]*dlambda0_dy)
			+ 2.0*(lambda[0]*dlambda21_dy + lambda[3]*dlambda0_dy);
	}
	else if (lambda[1] > 0) {
		val[0] = 4.0*(lambda[0]*dlambda20_dx + lambda[2]*dlambda0_dx);
		val[1] = 4.0*(lambda[0]*dlambda20_dy + lambda[2]*dlambda0_dy);
	}
	else {
		val[0] = 4.0*(lambda[0]*dlambda21_dx + lambda[3]*dlambda0_dx);
		val[1] = 4.0*(lambda[0]*dlambda21_dy + lambda[3]*dlambda0_dy);
	}
}

/*
*  end of file
**************************************************************************/
