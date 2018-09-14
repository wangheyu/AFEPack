/**
 * 本文件由邓剑提供。
 * 
 */

#include <assert.h>
#include <math.h>

void out_normal(const double * p, const double ** v, int i, double * n)
{
  double area;
  int fidx[4][3] = {
    {1, 3, 4}, {0, 4, 3}, {0, 1, 4}, {0, 3, 1}
  };
  const double * vtx[3] = {v[fidx[i][0]], v[fidx[i][1]], v[fidx[i][2]]};
  double d[2][3] = {
    {vtx[1][0] - vtx[0][0], vtx[1][1] - vtx[0][1], vtx[1][2] - vtx[0][2]},
    {vtx[2][0] - vtx[0][0], vtx[2][1] - vtx[0][1], vtx[2][2] - vtx[0][2]}
  };
  n[0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  n[1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  n[2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
  area = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] /= area; n[1] /= area; n[2] /= area;
}

/**
 * end of file
 * 
 */
