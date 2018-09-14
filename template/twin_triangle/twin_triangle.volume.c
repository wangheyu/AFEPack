#ifdef __cplusplus
extern "C" {
#endif

  double twin_triangle_volume(double ** v)
  {
    return .5*((v[1][0] - v[0][0])*(v[3][1] - v[0][1]) -
               (v[3][0] - v[0][0])*(v[1][1] - v[0][1]));
  }

#ifdef __cplusplus
}
#endif
