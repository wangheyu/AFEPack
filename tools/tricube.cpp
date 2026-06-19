/* tricube.cpp
by Robert Lie, September 28th, 2000

trianglize a cube
	[x0, x1] * [y0, y1] * [z0, z1]
by
	nx * ny * nz
partions.

every cubic cell is divided into 6 simplex.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double x0, x1;
double y0, y1;
double z0, z1;

int nx, ny, nz;

double dx, dy, dz;

int main(int argc, char * argv[])
{
	if (argc != 10) {
		printf("Usage: %s x0 x1 nx y0 y1 ny z0 z1 nz\n",
			argv[0]);
		exit(1);
	}

	x0 = atof(argv[1]);
	x1 = atof(argv[2]);
	nx = atoi(argv[3]);

	y0 = atof(argv[4]);
	y1 = atof(argv[5]);
	ny = atoi(argv[6]);

	z0 = atof(argv[7]);
	z1 = atof(argv[8]);
	nz = atoi(argv[9]);

	dx = (x1 - x0) / nx;
	dy = (y1 - y0) / ny;
	dz = (z1 - z0) / nz;

	FILE * fp0 = fopen("mesh0.dat", "w");
	assert(fp0 != NULL);
	FILE * fp3 = fopen("mesh3.dat", "w");
	assert(fp3 != NULL);

	int i, j, k, n, f;
	fprintf(fp0, "%d\n", (nx+1)*(ny+1)*(nz+1));
	for (i = 0;i <= nx;i ++) {
		for (j = 0;j <= ny;j ++) {
			for (k = 0;k <= nz;k ++) {
				f = i * (i-nx) * j * (j-ny) * k * (k-nz);
				n = (i * (nx+1) + j) * (ny+1) + k + 1;
				fprintf(fp0, "%d\t%f\t%f\t%f\t%d\n", n, x0+i*dx, y0+j*dy, z0+k*dz, !f);
			}
		}
	}
	fclose(fp0);

	int m[2][2][2];
	int i1, j1, k1;
	fprintf(fp3, "%d\n", 6*nx*ny*nz);
	for (i = 0;i < nx;i ++) {
		for(j = 0;j < ny;j ++) {
			for (k = 0;k < nz;k ++) {
				n = (i * nx + j) * ny + k;
				n = 6 * n + 1;
				for (i1 = 0;i1 < 2;i1 ++) {
					for (j1 = 0;j1 < 2;j1 ++) {
						for (k1 = 0;k1 < 2;k1 ++) {
							m[i1][j1][k1] = ((i+i1)*(nx+1)+(j+j1))*(ny+1)+k+k1+1;
						}
					}
				}
				fprintf(fp3,"%d\t%d\t%d\t%d\t%d\n",n++,m[0][0][0],m[0][0][1],m[1][0][0],m[0][1][0]);
				fprintf(fp3,"%d\t%d\t%d\t%d\t%d\n",n++,m[1][0][0],m[0][1][0],m[0][0][1],m[1][0][1]);
				fprintf(fp3,"%d\t%d\t%d\t%d\t%d\n",n++,m[0][1][1],m[0][1][0],m[1][0][1],m[0][0][1]);
				fprintf(fp3,"%d\t%d\t%d\t%d\t%d\n",n++,m[1][0][0],m[1][1][0],m[0][1][0],m[1][0][1]);
				fprintf(fp3,"%d\t%d\t%d\t%d\t%d\n",n++,m[0][1][1],m[1][1][0],m[1][0][1],m[0][1][0]);
				fprintf(fp3,"%d\t%d\t%d\t%d\t%d\n",n++,m[1][1][1],m[1][1][0],m[1][0][1],m[0][1][1]);
			}
		}
	}
	fclose(fp3);

	return 1;
}
