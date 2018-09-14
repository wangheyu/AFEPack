/*************************************************************************
* 	interval.4.bas_fun.c : by Huayi Wei
*/

void lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] = (val[0]-1.0)*(2.0*val[0]+1.0)*val[0]*(2.0*val[0]-1.0)/6.0 ;

}

void lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] =  (val[0]+1.0)*(2.0*val[0]+1.0)*val[0]*(2.0*val[0]-1.0)/6.0;
}
void lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] =  -4.0*(val[0]+1.0)*(val[0]-1.0)*val[0]*(2.0*val[0]-1.0)/3.0 ;
}

void lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] = (val[0]+1.0)*(val[0]-1.0)*(2.0*val[0]+1.0)*(2.0*val[0]-1.0) ;
}
void lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] =   -4.0*(val[0]+1.0)*(val[0]-1.0)*(2.0*val[0]+1.0)*val[0]/3.0 ;
}

void gradient_lambda_1(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] = (((16.0*val[0]-12.0)*val[0]-2.0)*val[0]+1.0)/6.0 ;
	val[0] = 2.0*val[0]/(v[1][0]-v[0][0]);
}

void gradient_lambda_2(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] =(((16.0*val[0]+12.0)*val[0]-2.0)*val[0]-1.0)/6.0 ;
	val[0] = 2.0*val[0]/(v[1][0]-v[0][0]);
}
void gradient_lambda_3(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] =(((-32.0*val[0]+12.0)*val[0]+16.0)*val[0]-4.0)/3.0;
	val[0] = 2.0*val[0]/(v[1][0]-v[0][0]);
}

void gradient_lambda_4(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] = (16.0*val[0]*val[0]-10.0)*val[0] ;
	val[0] = 2.0*val[0]/(v[1][0]-v[0][0]);
}

void gradient_lambda_5(const double * p, const double ** v, void * value)
{
	double * val = (double *)value;
        val[0] = (2.0*p[0]-v[0][0]-v[1][0])/(v[1][0]-v[0][0]);
	val[0] = (((-32.0*val[0]-12.0)*val[0]+16.0)*val[0]+4.0)/3.0;
	val[0] = 2.0*val[0]/(v[1][0]-v[0][0]);
}
/*
*  end of file
**************************************************************************/
