/**
 * @file   HGeometry.templates.nd.h
 * @author Robert Lie
 * @date   Sun Apr 29 13:24:24 2007
 * 
 * @brief  
 * 
 * 
 */

template <>
HGeometry<0,DOW>::HGeometry()
  : index(0), bmark(0) {}


template <>
HGeometry<1,DOW>::HGeometry() 
: index(0),
  vertex(n_vertex, NULL),
  parent(NULL),
  child(n_child, NULL),
  bmark(0) {}



template <>
HGeometry<2,DOW>::HGeometry()
: index(0),
  vertex(n_vertex, NULL),
  boundary(n_boundary, NULL),
  parent(NULL),
  child(n_child, NULL),
  bmark(0) {}



template <>
HGeometry<3,DOW>::HGeometry()
: index(0),
  vertex(n_vertex, NULL),
  boundary(n_boundary, NULL),
  parent(NULL),
  child(n_child, NULL),
  bmark(0) {}


template <>
bool HGeometry<1,DOW>::isRefined() const
{
  return (child.size() > 0 && child[0] != NULL);
}



template <>
bool HGeometry<2,DOW>::isRefined() const
{
  return (child.size() > 0 && child[0] != NULL);
}



template <>
bool HGeometry<3,DOW>::isRefined() const
{
  return (child.size() > 0 && child[0] != NULL);
}



template <>
void HGeometry<1,DOW>::refine()
{
  if (isRefined()) return;
	
  // add the midpoint of the edge at first
  HGeometry<0,DOW> * new_point = new HGeometry<0,DOW>();
  Assert(new_point != NULL, ExcOutOfMemory());
  if (mid_point == NULL) {
    *dynamic_cast<Point<DOW> *>(new_point) = midpoint(*vertex[0], *vertex[1]);
  }
  else {
    (*mid_point)(*vertex[0], *vertex[1], bmark, *new_point);
  }
  new_point->bmark = bmark;

  // construct the 0th child for the edge
  child[0] = new HGeometry<1,DOW>();
  Assert(child[0] != NULL, ExcOutOfMemory());
  child[0]->parent = this;
  child[0]->vertex[0] = vertex[0];
  child[0]->vertex[1] = new_point;
  child[0]->bmark = bmark;

  // construct the 1st child for the edge
  child[1] = new HGeometry<1,DOW>();
  Assert(child[1] != NULL, ExcOutOfMemory());
  child[1]->parent = this;
  child[1]->vertex[0] = new_point;
  child[1]->vertex[1] = vertex[1];
  child[1]->bmark = bmark;
}



template <>
void HGeometry<2,DOW>::refine()
{
  int i, ii[]={0, 1, 2, 0, 1, 2};
  if (isRefined()) return;

  // refine its edge at first
  for (i = 0;i < 3;i ++) 
    boundary[i]->refine();
	
  // add the three new edge in the triangle
  HGeometry<0,DOW> * edge_midpoint[3];
  for (i = 0;i < 3;i ++) {
    edge_midpoint[i] = boundary[i]->child[0]->vertex[1];
  }
  HGeometry<1,DOW> * new_edge[3];
  for (i = 0;i < 3;i ++) {
    new_edge[i] = new HGeometry<1,DOW>();
    Assert(new_edge[i] != NULL, ExcOutOfMemory());
    new_edge[i]->vertex[0] = edge_midpoint[ii[i+1]];
    new_edge[i]->vertex[1] = edge_midpoint[ii[i+2]];
    new_edge[i]->bmark = bmark;
  }
	
  // construct the 0th to 2nd child for the triangle
  for (i = 0;i < 3;i ++) {
    if (typeid(*this) == typeid(HGeometry<2, DOW>))
      child[i] = (HGeometry<2,DOW> *) new HGeometry<2, DOW>();
    else
      child[i] = new HGeometry<2,DOW>();
    Assert(child[i] != NULL, ExcOutOfMemory());
    child[i]->parent = this;
    child[i]->vertex[0] = vertex[i];
    child[i]->vertex[1] = edge_midpoint[ii[i+2]];
    child[i]->vertex[2] = edge_midpoint[ii[i+1]];
    child[i]->boundary[0] = new_edge[i];
    if (boundary[ii[i+1]]->vertex[0] == vertex[i]) 
      child[i]->boundary[1] = boundary[ii[i+1]]->child[0];
    else 
      child[i]->boundary[1] = boundary[ii[i+1]]->child[1];
    if (boundary[ii[i+2]]->vertex[0] == vertex[i]) 
      child[i]->boundary[2] = boundary[ii[i+2]]->child[0];
    else 
      child[i]->boundary[2] = boundary[ii[i+2]]->child[1];
    child[i]->bmark = bmark;
  }
	
  // construct the 3rd child, the center child, for the triangle
  child[3] = new HGeometry<2,DOW>();
  Assert(child[3] != NULL, ExcOutOfMemory());
  child[3]->parent = this;
  child[3]->vertex[0] = edge_midpoint[0];
  child[3]->vertex[1] = edge_midpoint[1];
  child[3]->vertex[2] = edge_midpoint[2];
  child[3]->boundary[0] = new_edge[0];
  child[3]->boundary[1] = new_edge[1];
  child[3]->boundary[2] = new_edge[2];
  child[3]->bmark = bmark;
}



template <>
void HGeometry<3,DOW>::refine()
{
  int i, ii[]={0, 1, 2, 0, 1, 2, 0, 1, 2};
  if (isRefined()) return;

  // refine its surfaces at first. After the surfaces of the tetrahedron are refined,
  // all those edges of the tetrahedron are refined, too.
  for (i = 0;i < HGeometry<3,DOW>::n_boundary;i ++) 
    boundary[i]->refine();

  // retrieve the information of those midpoints of the edges of the tetrahedron. 
  // Because the data of the surfaces of the tetrahedron are not unique for the
  // tetrahedron, it's not so simple to get the information of those midpoints.
  // A lot of codes are written to get the relationship detail between the surfaces
  // and the tetrahedron. Be patiently and let's do it step by step,
  HGeometry<0,DOW> * edge_midpoint[6];
  HGeometry<2,DOW> * surface_child[4][4];
	
  {
    // for the first surface of the tetrahedron, which is on the bottom of the tetrahedron,
    // we should judge the pose of the triangle at first.
    // boundary[0]: the first surface
    // boundary[0]->vertex[i]: the i-th vertex of the triangle
    int start, step;
    for (start = 0;start < 3; start ++)
      if (boundary[0]->vertex[start] == vertex[1])
	break;
    Assert (start < 3,ExcInternalError());
    if (boundary[0]->vertex[ii[start+1]] == vertex[2])
      step = 1;
    else {
      Assert(boundary[0]->vertex[ii[start+1]] == vertex[3], ExcInternalError());
      Assert(boundary[0]->vertex[ii[start+2]] == vertex[2], ExcInternalError());
      step = 2;
    }
    edge_midpoint[3] = boundary[0]->boundary[start]->child[0]->vertex[1];
    edge_midpoint[4] = boundary[0]->boundary[ii[start+step]]->child[0]->vertex[1];
    edge_midpoint[5] = boundary[0]->boundary[ii[start+2*step]]->child[0]->vertex[1];
		
    surface_child[0][0] = boundary[0]->child[start];
    surface_child[0][1] = boundary[0]->child[ii[start+step]];
    surface_child[0][2] = boundary[0]->child[ii[start+2*step]];
    surface_child[0][3] = boundary[0]->child[3];

    // for the second surface of the tetrahedron, we should judge its pose at first, too.
    // boundary[1]: the second surface
    // boundary[1]->vertex[i]: the i-th vertex of the triangle
    for (start = 0;start < 3; start ++)
      if (boundary[1]->vertex[start] == vertex[0])
	break;
    Assert (start < 3,ExcInternalError());
    if (boundary[1]->vertex[ii[start+1]] == vertex[2])
      step = 1;
    else {
      Assert(boundary[1]->vertex[ii[start+1]] == vertex[3], ExcInternalError());
      Assert(boundary[1]->vertex[ii[start+2]] == vertex[2], ExcInternalError());
      step = 2;
    }
    Assert(edge_midpoint[3] == boundary[1]->boundary[start]->child[0]->vertex[1], ExcInternalError());
    edge_midpoint[2] = boundary[1]->boundary[ii[start+step]]->child[0]->vertex[1];
    edge_midpoint[1] = boundary[1]->boundary[ii[start+2*step]]->child[0]->vertex[1];

    surface_child[1][0] = boundary[1]->child[start];
    surface_child[1][1] = boundary[1]->child[ii[start+step]];
    surface_child[1][2] = boundary[1]->child[ii[start+2*step]];
    surface_child[1][3] = boundary[1]->child[3];
		
    // for the third surface of the tetrahedron, we should judge its pose at first, too.
    // boundary[2]: the third surface
    // boundary[2]->vertex[i]: the i-th vertex of the triangle
    for (start = 0;start < 3; start ++)
      if (boundary[2]->vertex[start] == vertex[0])
	break;
    Assert (start < 3,ExcInternalError());
    if (boundary[2]->vertex[ii[start+1]] == vertex[1])
      step = 1;
    else {
      Assert(boundary[2]->vertex[ii[start+1]] == vertex[3], ExcInternalError());
      Assert(boundary[2]->vertex[ii[start+2]] == vertex[1], ExcInternalError());
      step = 2;
    }
    Assert(edge_midpoint[4] == boundary[2]->boundary[start]->child[0]->vertex[1], ExcInternalError());
    Assert(edge_midpoint[2] == boundary[2]->boundary[ii[start+step]]->child[0]->vertex[1], ExcInternalError());
    edge_midpoint[0] = boundary[2]->boundary[ii[start+2*step]]->child[0]->vertex[1];

    surface_child[2][0] = boundary[2]->child[start];
    surface_child[2][1] = boundary[2]->child[ii[start+step]];
    surface_child[2][2] = boundary[2]->child[ii[start+2*step]];
    surface_child[2][3] = boundary[2]->child[3];

    // for the fourth surface of the tetrahedron, we should judge its pose at first, too.
    // boundary[3]: the fourth surface
    // boundary[3]->vertex[i]: the i-th vertex of the triangle
    for (start = 0;start < 3; start ++)
      if (boundary[3]->vertex[start] == vertex[0])
	break;
    Assert (start < 3,ExcInternalError());
    if (boundary[3]->vertex[ii[start+1]] == vertex[1])
      step = 1;
    else {
      Assert(boundary[3]->vertex[ii[start+1]] == vertex[2], ExcInternalError());
      Assert(boundary[3]->vertex[ii[start+2]] == vertex[1], ExcInternalError());
      step = 2;
    }
    Assert(edge_midpoint[5] == boundary[3]->boundary[start]->child[0]->vertex[1], ExcInternalError());
    Assert(edge_midpoint[1] == boundary[3]->boundary[ii[start+step]]->child[0]->vertex[1], ExcInternalError());
    Assert(edge_midpoint[0] == boundary[3]->boundary[ii[start+2*step]]->child[0]->vertex[1], ExcInternalError());

    surface_child[3][0] = boundary[3]->child[start];
    surface_child[3][1] = boundary[3]->child[ii[start+step]];
    surface_child[3][2] = boundary[3]->child[ii[start+2*step]];
    surface_child[3][3] = boundary[3]->child[3];
  }

  // we add those geometries should be added for all refine model at first, which
  // include four triangles and four tetrahedrons at the four vertices.
  HGeometry<2,DOW> * triangle[8];
  for (i = 0;i < 8;i ++) {
    triangle[i] = new HGeometry<2,DOW>();
    Assert(triangle[i] != NULL, ExcOutOfMemory());
  }
  triangle[0]->vertex[0] = edge_midpoint[0];
  triangle[0]->vertex[1] = edge_midpoint[1];
  triangle[0]->vertex[2] = edge_midpoint[2];
  triangle[0]->boundary[0] = surface_child[1][0]->boundary[0];
  triangle[0]->boundary[1] = surface_child[2][0]->boundary[0];
  triangle[0]->boundary[2] = surface_child[3][0]->boundary[0];
  triangle[0]->bmark = bmark;
	
  triangle[1]->vertex[0] = edge_midpoint[4];
  triangle[1]->vertex[1] = edge_midpoint[5];
  triangle[1]->vertex[2] = edge_midpoint[0];
  triangle[1]->boundary[0] = surface_child[3][1]->boundary[0];
  triangle[1]->boundary[1] = surface_child[2][1]->boundary[0];
  triangle[1]->boundary[2] = surface_child[0][0]->boundary[0];	
  triangle[1]->bmark = bmark;
	
  triangle[2]->vertex[0] = edge_midpoint[5];
  triangle[2]->vertex[1] = edge_midpoint[3];
  triangle[2]->vertex[2] = edge_midpoint[1];
  triangle[2]->boundary[0] = surface_child[1][1]->boundary[0];
  triangle[2]->boundary[1] = surface_child[3][2]->boundary[0];
  triangle[2]->boundary[2] = surface_child[0][1]->boundary[0];	
  triangle[2]->bmark = bmark;

  triangle[3]->vertex[0] = edge_midpoint[2];
  triangle[3]->vertex[1] = edge_midpoint[3];
  triangle[3]->vertex[2] = edge_midpoint[4];
  triangle[3]->boundary[0] = surface_child[0][2]->boundary[0];
  triangle[3]->boundary[1] = surface_child[2][2]->boundary[0];
  triangle[3]->boundary[2] = surface_child[1][2]->boundary[0];	
  triangle[3]->bmark = bmark;

  // the four child tetrahedrons on the four verteics
  for (i = 0;i < 8;i ++) {
    child[i] = new HGeometry<3,DOW>();
    Assert (child[i] != NULL, ExcOutOfMemory());
  }
	
  child[0]->parent = this;
  child[0]->vertex[0] = vertex[0];
  child[0]->vertex[1] = edge_midpoint[0];
  child[0]->vertex[2] = edge_midpoint[1];
  child[0]->vertex[3] = edge_midpoint[2];
  child[0]->boundary[0] = triangle[0];
  child[0]->boundary[1] = surface_child[1][0];
  child[0]->boundary[2] = surface_child[2][0];
  child[0]->boundary[3] = surface_child[3][0];
  child[0]->bmark = bmark;
	
  child[1]->parent = this;
  child[1]->vertex[0] = vertex[1];
  child[1]->vertex[1] = edge_midpoint[5];
  child[1]->vertex[2] = edge_midpoint[0];
  child[1]->vertex[3] = edge_midpoint[4];
  child[1]->boundary[0] = triangle[1];
  child[1]->boundary[1] = surface_child[2][1];
  child[1]->boundary[2] = surface_child[0][0];
  child[1]->boundary[3] = surface_child[3][1];
  child[1]->bmark = bmark;
	
  child[2]->parent = this;
  child[2]->vertex[0] = vertex[2];
  child[2]->vertex[1] = edge_midpoint[1];
  child[2]->vertex[2] = edge_midpoint[5];
  child[2]->vertex[3] = edge_midpoint[3];
  child[2]->boundary[0] = triangle[2];
  child[2]->boundary[1] = surface_child[0][1];
  child[2]->boundary[2] = surface_child[1][1];
  child[2]->boundary[3] = surface_child[3][2];
  child[2]->bmark = bmark;
	
  child[3]->parent = this;
  child[3]->vertex[0] = vertex[3];
  child[3]->vertex[1] = edge_midpoint[2];
  child[3]->vertex[2] = edge_midpoint[3];
  child[3]->vertex[3] = edge_midpoint[4];
  child[3]->boundary[0] = triangle[3];
  child[3]->boundary[1] = surface_child[0][2];
  child[3]->boundary[2] = surface_child[2][2];
  child[3]->boundary[3] = surface_child[1][2];
  child[3]->bmark = bmark;	

  // choose the best refinement model. There are total 3 models to refine the
  // central part of the tetrahedron, between which the key difference is which
  // line is chosed as the main axis of geometry. To make the obtained tetrahedrons
  // more regular, we will choose the shortest diagonal line as the main axis.
  double l03, l14, l25;
  l03 = (*edge_midpoint[0] - *edge_midpoint[3]).length();
  l14 = (*edge_midpoint[1] - *edge_midpoint[4]).length();
  l25 = (*edge_midpoint[2] - *edge_midpoint[5]).length();
  if (l03 <= l14)
    if (l03 <= l25)
      refine_model = REFINE_MODEL_03;
    else
      refine_model = REFINE_MODEL_25;
  else
    if (l14 <= l25)
      refine_model = REFINE_MODEL_14;
    else
      refine_model = REFINE_MODEL_25;

  HGeometry<1,DOW> * axis;
  switch (refine_model) {
	
  case REFINE_MODEL_03:
    // add the main axis at first
    axis = new HGeometry<1,DOW>();
    Assert (axis != NULL, ExcOutOfMemory());
    axis->vertex[0] = edge_midpoint[0];
    axis->vertex[1] = edge_midpoint[3];
    axis->bmark = bmark;
		
    // add the four triangles
    triangle[4]->vertex[0] = edge_midpoint[4];
    triangle[4]->vertex[1] = edge_midpoint[3];
    triangle[4]->vertex[2] = edge_midpoint[0];
    triangle[4]->boundary[0] = axis;
    triangle[4]->boundary[1] = surface_child[2][1]->boundary[0];
    triangle[4]->boundary[2] = surface_child[0][2]->boundary[0];
    triangle[4]->bmark = bmark;
		
    triangle[5]->vertex[0] = edge_midpoint[2];
    triangle[5]->vertex[1] = edge_midpoint[3];
    triangle[5]->vertex[2] = edge_midpoint[0];
    triangle[5]->boundary[0] = axis;
    triangle[5]->boundary[1] = surface_child[2][0]->boundary[0];
    triangle[5]->boundary[2] = surface_child[1][2]->boundary[0];
    triangle[5]->bmark = bmark;
		
    triangle[6]->vertex[0] = edge_midpoint[1];
    triangle[6]->vertex[1] = edge_midpoint[3];
    triangle[6]->vertex[2] = edge_midpoint[0];
    triangle[6]->boundary[0] = axis;
    triangle[6]->boundary[1] = surface_child[3][0]->boundary[0];
    triangle[6]->boundary[2] = surface_child[1][1]->boundary[0];
    triangle[6]->bmark = bmark;

    triangle[7]->vertex[0] = edge_midpoint[5];
    triangle[7]->vertex[1] = edge_midpoint[3];
    triangle[7]->vertex[2] = edge_midpoint[0];
    triangle[7]->boundary[0] = axis;
    triangle[7]->boundary[1] = surface_child[3][1]->boundary[0];
    triangle[7]->boundary[2] = surface_child[0][1]->boundary[0];
    triangle[7]->bmark = bmark;
		
    // add the four tetrahedrons
    child[4]->parent = this;
    child[4]->vertex[0] = edge_midpoint[0];
    child[4]->vertex[1] = edge_midpoint[3];
    child[4]->vertex[2] = edge_midpoint[4];
    child[4]->vertex[3] = edge_midpoint[5];
    child[4]->boundary[0] = surface_child[0][3];
    child[4]->boundary[1] = triangle[1];
    child[4]->boundary[2] = triangle[7];
    child[4]->boundary[3] = triangle[4];
    child[4]->bmark = bmark;

    child[5]->parent = this;
    child[5]->vertex[0] = edge_midpoint[0];
    child[5]->vertex[1] = edge_midpoint[3];
    child[5]->vertex[2] = edge_midpoint[1];
    child[5]->vertex[3] = edge_midpoint[2];
    child[5]->boundary[0] = surface_child[1][3];
    child[5]->boundary[1] = triangle[0];
    child[5]->boundary[2] = triangle[5];
    child[5]->boundary[3] = triangle[6];
    child[5]->bmark = bmark;

    child[6]->parent = this;
    child[6]->vertex[0] = edge_midpoint[3];
    child[6]->vertex[1] = edge_midpoint[0];
    child[6]->vertex[2] = edge_midpoint[4];
    child[6]->vertex[3] = edge_midpoint[2];
    child[6]->boundary[0] = surface_child[2][3];
    child[6]->boundary[1] = triangle[3];
    child[6]->boundary[2] = triangle[5];
    child[6]->boundary[3] = triangle[4];
    child[6]->bmark = bmark;

    child[7]->parent = this;
    child[7]->vertex[0] = edge_midpoint[3];
    child[7]->vertex[1] = edge_midpoint[0];
    child[7]->vertex[2] = edge_midpoint[1];
    child[7]->vertex[3] = edge_midpoint[5];
    child[7]->boundary[0] = surface_child[3][3];
    child[7]->boundary[1] = triangle[2];
    child[7]->boundary[2] = triangle[7];
    child[7]->boundary[3] = triangle[6];
    child[7]->bmark = bmark;

    break;

  case REFINE_MODEL_14:
    // add the main axis at first
    axis = new HGeometry<1,DOW>();
    Assert (axis != NULL, ExcOutOfMemory());
    axis->vertex[0] = edge_midpoint[1];
    axis->vertex[1] = edge_midpoint[4];
    axis->bmark = bmark;
		
    // add the four triangles
    triangle[4]->vertex[0] = edge_midpoint[5];
    triangle[4]->vertex[1] = edge_midpoint[4];
    triangle[4]->vertex[2] = edge_midpoint[1];
    triangle[4]->boundary[0] = axis;
    triangle[4]->boundary[1] = surface_child[3][2]->boundary[0];
    triangle[4]->boundary[2] = surface_child[0][0]->boundary[0];
    triangle[4]->bmark = bmark;
		
    triangle[5]->vertex[0] = edge_midpoint[0];
    triangle[5]->vertex[1] = edge_midpoint[4];
    triangle[5]->vertex[2] = edge_midpoint[1];
    triangle[5]->boundary[0] = axis;
    triangle[5]->boundary[1] = surface_child[3][0]->boundary[0];
    triangle[5]->boundary[2] = surface_child[2][1]->boundary[0];
    triangle[5]->bmark = bmark;
		
    triangle[6]->vertex[0] = edge_midpoint[2];
    triangle[6]->vertex[1] = edge_midpoint[4];
    triangle[6]->vertex[2] = edge_midpoint[1];
    triangle[6]->boundary[0] = axis;
    triangle[6]->boundary[1] = surface_child[1][0]->boundary[0];
    triangle[6]->boundary[2] = surface_child[2][2]->boundary[0];
    triangle[6]->bmark = bmark;

    triangle[7]->vertex[0] = edge_midpoint[3];
    triangle[7]->vertex[1] = edge_midpoint[4];
    triangle[7]->vertex[2] = edge_midpoint[1];
    triangle[7]->boundary[0] = axis;
    triangle[7]->boundary[1] = surface_child[1][1]->boundary[0];
    triangle[7]->boundary[2] = surface_child[0][2]->boundary[0];
    triangle[7]->bmark = bmark;
		
    // add the four tetrahedrons
    child[4]->parent = this;
    child[4]->vertex[0] = edge_midpoint[1];
    child[4]->vertex[1] = edge_midpoint[4];
    child[4]->vertex[2] = edge_midpoint[5];
    child[4]->vertex[3] = edge_midpoint[3];
    child[4]->boundary[0] = surface_child[0][3];
    child[4]->boundary[1] = triangle[2];
    child[4]->boundary[2] = triangle[7];
    child[4]->boundary[3] = triangle[4];
    child[4]->bmark = bmark;

    child[5]->parent = this;
    child[5]->vertex[0] = edge_midpoint[4];
    child[5]->vertex[1] = edge_midpoint[1];
    child[5]->vertex[2] = edge_midpoint[2];
    child[5]->vertex[3] = edge_midpoint[3];
    child[5]->boundary[0] = surface_child[1][3];
    child[5]->boundary[1] = triangle[3];
    child[5]->boundary[2] = triangle[7];
    child[5]->boundary[3] = triangle[6];
    child[5]->bmark = bmark;

    child[6]->parent = this;
    child[6]->vertex[0] = edge_midpoint[1];
    child[6]->vertex[1] = edge_midpoint[4];
    child[6]->vertex[2] = edge_midpoint[2];
    child[6]->vertex[3] = edge_midpoint[0];
    child[6]->boundary[0] = surface_child[2][3];
    child[6]->boundary[1] = triangle[0];
    child[6]->boundary[2] = triangle[5];
    child[6]->boundary[3] = triangle[6];
    child[6]->bmark = bmark;

    child[7]->parent = this;
    child[7]->vertex[0] = edge_midpoint[4];
    child[7]->vertex[1] = edge_midpoint[1];
    child[7]->vertex[2] = edge_midpoint[5];
    child[7]->vertex[3] = edge_midpoint[0];
    child[7]->boundary[0] = surface_child[3][3];
    child[7]->boundary[1] = triangle[1];
    child[7]->boundary[2] = triangle[5];
    child[7]->boundary[3] = triangle[4];
    child[7]->bmark = bmark;

    break;

  case REFINE_MODEL_25:
    // add the main axis at first
    axis = new HGeometry<1,DOW>();
    Assert (axis != NULL, ExcOutOfMemory());
    axis->vertex[0] = edge_midpoint[2];
    axis->vertex[1] = edge_midpoint[5];
    axis->bmark = bmark;
		
    // add the four triangles
    triangle[4]->vertex[0] = edge_midpoint[4];
    triangle[4]->vertex[1] = edge_midpoint[5];
    triangle[4]->vertex[2] = edge_midpoint[2];
    triangle[4]->boundary[0] = axis;
    triangle[4]->boundary[1] = surface_child[2][2]->boundary[0];
    triangle[4]->boundary[2] = surface_child[0][0]->boundary[0];
    triangle[4]->bmark = bmark;
		
    triangle[5]->vertex[0] = edge_midpoint[0];
    triangle[5]->vertex[1] = edge_midpoint[5];
    triangle[5]->vertex[2] = edge_midpoint[2];
    triangle[5]->boundary[0] = axis;
    triangle[5]->boundary[1] = surface_child[2][0]->boundary[0];
    triangle[5]->boundary[2] = surface_child[3][1]->boundary[0];
    triangle[5]->bmark = bmark;
		
    triangle[6]->vertex[0] = edge_midpoint[1];
    triangle[6]->vertex[1] = edge_midpoint[5];
    triangle[6]->vertex[2] = edge_midpoint[2];
    triangle[6]->boundary[0] = axis;
    triangle[6]->boundary[1] = surface_child[1][0]->boundary[0];
    triangle[6]->boundary[2] = surface_child[3][2]->boundary[0];
    triangle[6]->bmark = bmark;

    triangle[7]->vertex[0] = edge_midpoint[3];
    triangle[7]->vertex[1] = edge_midpoint[5];
    triangle[7]->vertex[2] = edge_midpoint[2];
    triangle[7]->boundary[0] = axis;
    triangle[7]->boundary[1] = surface_child[1][2]->boundary[0];
    triangle[7]->boundary[2] = surface_child[0][1]->boundary[0];
    triangle[7]->bmark = bmark;
		
    // add the four tetrahedrons
    child[4]->parent = this;
    child[4]->vertex[0] = edge_midpoint[2];
    child[4]->vertex[1] = edge_midpoint[5];
    child[4]->vertex[2] = edge_midpoint[3];
    child[4]->vertex[3] = edge_midpoint[4];
    child[4]->boundary[0] = surface_child[0][3];
    child[4]->boundary[1] = triangle[3];
    child[4]->boundary[2] = triangle[4];
    child[4]->boundary[3] = triangle[7];
    child[4]->bmark = bmark;

    child[5]->parent = this;
    child[5]->vertex[0] = edge_midpoint[5];
    child[5]->vertex[1] = edge_midpoint[2];
    child[5]->vertex[2] = edge_midpoint[3];
    child[5]->vertex[3] = edge_midpoint[1];
    child[5]->boundary[0] = surface_child[1][3];
    child[5]->boundary[1] = triangle[2];
    child[5]->boundary[2] = triangle[6];
    child[5]->boundary[3] = triangle[7];
    child[5]->bmark = bmark;

    child[6]->parent = this;
    child[6]->vertex[0] = edge_midpoint[5];
    child[6]->vertex[1] = edge_midpoint[2];
    child[6]->vertex[2] = edge_midpoint[0];
    child[6]->vertex[3] = edge_midpoint[4];
    child[6]->boundary[0] = surface_child[2][3];
    child[6]->boundary[1] = triangle[1];
    child[6]->boundary[2] = triangle[4];
    child[6]->boundary[3] = triangle[5];
    child[6]->bmark = bmark;

    child[7]->parent = this;
    child[7]->vertex[0] = edge_midpoint[2];
    child[7]->vertex[1] = edge_midpoint[5];
    child[7]->vertex[2] = edge_midpoint[0];
    child[7]->vertex[3] = edge_midpoint[1];
    child[7]->boundary[0] = surface_child[3][3];
    child[7]->boundary[1] = triangle[0];
    child[7]->boundary[2] = triangle[6];
    child[7]->boundary[3] = triangle[5];
    child[7]->bmark = bmark;

    break;

  default:
    Assert(false, ExcInternalError());
  }	
}


template <>
void HGeometry<1,DOW>::checkIntegrity() const
{
  if (!isRefined()) return;
  for (int i = 0;i < n_child;i ++) {
    child[i]->checkIntegrity();
  }
}



template <>
void HGeometry<2,DOW>::checkIntegrity() const
{
  int i, j, k;
  for (i = 0;i < n_boundary;i ++) {
    HGeometry<1,DOW> * b = boundary[i];
    for (j = 0;j < b->n_vertex;j ++) {
      for (k = 0;k < n_vertex;k ++) {
	if (k == i) continue;
	if (b->vertex[j] == vertex[k])
	  break;
      }
      Assert(k < n_vertex, ExcInternalError());
    }
  }
  if (!isRefined()) return;
  for (int i = 0;i < n_child;i ++) {
    child[i]->checkIntegrity();	
  }
}



template <>
void HGeometry<3,DOW>::checkIntegrity() const
{
  int i, j, k;
  for (i = 0;i < n_boundary;i ++) {
    HGeometry<2,DOW> * b = boundary[i];
    for (j = 0;j < b->n_vertex;j ++) {
      for (k = 0;k < n_vertex;k ++) {
	if (k == i) continue;
	if (b->vertex[j] == vertex[k])
	  break;
      }
      Assert(k < n_vertex, ExcInternalError());
    }
  }
  if (!isRefined()) return;
  for (int i = 0;i < n_child;i ++) {
    child[i]->checkIntegrity();	
  }
}



template <>
bool HGeometry<1,DOW>::isIncludePoint(const Point<DOW>& p) const
{
  std::cerr << "warning: HGeometry<1,DOW>::isIncludePoint not implemented, return true" << std::endl;
  return true;
}


template <>
bool HGeometry<2,DOW>::isIncludePoint(const Point<DOW>& p) const
{
  double lambda[3];
  const std::vector<vertex_t *>& v = vertex;
  double volume = (((*v[1])[0] - (*v[0])[0])*((*v[2])[1] - (*v[0])[1]) - 
                   ((*v[1])[1] - (*v[0])[1])*((*v[2])[0] - (*v[0])[0]));
  double mzero = -1.0e-08;
  lambda[0] = (((*v[1])[0] - p[0])*((*v[2])[1] - p[1]) - 
               ((*v[1])[1] - p[1])*((*v[2])[0] - p[0]))/volume;
  lambda[1] = (((*v[2])[0] - p[0])*((*v[0])[1] - p[1]) - 
               ((*v[2])[1] - p[1])*((*v[0])[0] - p[0]))/volume;
  lambda[2] = (((*v[0])[0] - p[0])*((*v[1])[1] - p[1]) -
               ((*v[0])[1] - p[1])*((*v[1])[0] - p[0]))/volume;
  return (lambda[0] >= mzero && lambda[1] >= mzero && lambda[2] >= mzero);
}



template <>
bool HGeometry<3,DOW>::isIncludePoint(const Point<DOW>& p) const
{
  double ZERO = 1.0e-06; /// 似乎不能设置得太接近机器精度
  const std::vector<vertex_t *>& v = vertex;
  double A[4] = {
    tetrahedron_volume(    p, *v[1], *v[2], *v[3]),
    tetrahedron_volume(*v[0],     p, *v[2], *v[3]),
    tetrahedron_volume(*v[0], *v[1],     p, *v[3]),
    tetrahedron_volume(*v[0], *v[1], *v[2],     p)
  };
  double AZERO = ZERO*(A[0] + A[1] + A[2] + A[3]);
  return ((A[0] >= -AZERO) && (A[1] >= -AZERO) && 
          (A[2] >= -AZERO) && (A[3] >= -AZERO));
}



template <>
void HElement<1, DOW>::refine()
{
  if (isRefined()) return;

  h_element->refine();

  // construct its own child
  for (int i = 0;i < n_child;i ++) {
    child[i] = new HElement<1, DOW>();
    Assert (child[i] != NULL, ExcOutOfMemory());
    child[i]->h_element = h_element->child[i];
    child[i]->parent = this;
  }
}



template <>
void HElement<2, DOW>::refine()
{
  if (isRefined()) return;

  h_element->refine();

  // construct its own child
  for (int i = 0;i < n_child;i ++) {
    child[i] = new HElement<2, DOW>();
    Assert (child[i] != NULL, ExcOutOfMemory());
    child[i]->h_element = h_element->child[i];
    child[i]->parent = this;
  }
}



template <>
void HElement<3, DOW>::refine()
{
  if (isRefined()) return;
	
  // refine the h geometry tree at first
  h_element->refine();

  // construct its own child
  for (int i = 0;i < n_child;i ++) {
    child[i] = new HElement<3, DOW>();
    Assert (child[i] != NULL, ExcOutOfMemory());
    child[i]->h_element = h_element->child[i];
    child[i]->parent = this;
  }
}



template <>
std::ostream& operator<<(std::ostream& os, const HGeometry<1,DOW>& geometry)
{
  os << (*geometry.vertex[0])[0] << "\t"
     << (*geometry.vertex[1])[0] << "\t"
     << (*geometry.vertex[0])[1] << "\t"
     << (*geometry.vertex[1])[1] << "\n";
  return os;
}



template <>
std::ostream& operator<<(std::ostream& os, const HGeometry<0,DOW>& geometry)
{
  //	for (int i = 0;i < DOW;i ++)
  //		os << geometry[i] << "\t";
  return os;
}



template <>
void IrregularMesh<1, DOW>::regularize(bool renumerate)
{
  if (regular_mesh != NULL)
    delete regular_mesh;
  std::cerr << "Generating regular mesh from the semiregular mesh ..." << std::flush; 
  int i, j, k, n_node, n_element;

  std::cerr << "\n\tpreparing data ..." << std::flush;
  Tools tools;
  ActiveElementIterator<1, DOW> 
    the_ele = beginActiveElement(),
    end_ele = endActiveElement();
  // 为对几何体进行排序进行准备
  n_element = 0;
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<1, DOW> *& h_element = the_ele->h_element;
    /**
     * 半正则化的时候对点的index未做设置，而每个active的单元顶点必将在
     * 最后的网格中出现，所以在此将其设置为active状态。
     */
    for (i = 0;i < h_element->n_vertex;i ++) {
      tools.setGeometryActive(*(h_element->vertex[i]));
    }
    assert (tools.isGeometryUsed(*h_element));
    tools.setParentInactive(*h_element);
    n_element += 1;
  }

  // 对所有几何体进行排序
  n_node = 0;
  n_element = 0;
  the_ele = beginActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<1, DOW> *& h_element = the_ele->h_element;
    /**
     * 对顶点进行编号，已经编号的顶点将具有非负的index，从而不会满足
     * active判据，这样使得每个顶点只是被编号一次。
     */
    for (i = 0;i < the_ele->n_vertex;i ++) { /// 排顶点
      HGeometry<0, DOW>& node = *(h_element->vertex[i]);
      if (tools.isGeometryActive(node)) {
	node.index = n_node ++;
      }
    }
    assert (tools.isGeometryActive(*h_element));
    h_element->index = n_element ++; /// 排单元
    the_ele->index = h_element->index;
  }
	
  std::cerr << "\n\tbuilding the regular mesh ..." << std::flush;
  // 建立正则化的网格
  regular_mesh = new RegularMesh<1, DOW>(this);
  Assert (regular_mesh != NULL, ExcOutOfMemory());
  regular_mesh->point().resize(n_node);
  regular_mesh->geometry(0).resize(n_node);
  regular_mesh->geometry(1).resize(n_element);

#ifdef __SERIALIZATION__
  std::vector<std::vector<HGeometryBase*> >& h_geometry = regular_mesh->h_geometry();
  h_geometry.clear();
  h_geometry.resize(2);
  h_geometry[0].resize(n_node);
  h_geometry[1].resize(n_element);
#endif

  the_ele = beginActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<1, DOW> * h_element = the_ele->h_element;
    // 先加入线段的顶点
    for (i = 0;i < the_ele->n_vertex;i ++) {
      HGeometry<0, DOW>& vtx = *(h_element->vertex[i]);
      j = vtx.index;
      GeometryBM& pnt = regular_mesh->geometry(0,j);
      tools.regularize_add_node(vtx, 
                                pnt,
                                regular_mesh->point(j));
#ifdef __SERIALIZATION__
      h_geometry[0][j] = &vtx;
#endif
    }

    // 最后是线段自身
    j = h_element->index;
    assert (j == the_ele->index);
    GeometryBM& element = regular_mesh->geometry(1,j);
    Assert(element.index() == -1, ExcInternalError());
#ifdef __SERIALIZATION__
    h_geometry[1][j] = h_element;
#endif
    tools.regularize_add_side(*h_element, element);
  }
	
  geometry_tree->unlock(); /// 对几何遗传树解锁
  std::cerr << "\n\tnodes: " << n_node 
	    << "; elements: " << n_element << std::endl;
	
  if (renumerate) renumerateElement();
}



template <>
void IrregularMesh<2, DOW>::regularize(bool renumerate)
{
  if (regular_mesh != NULL)
    delete regular_mesh;
  std::cerr << "Generating regular mesh from the semiregular mesh ..." << std::flush; 
  int i, j, k, n_node, n_side, n_element;
  int ii[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};

  std::cerr << "\n\tpreparing data ..." << std::flush;
  Tools tools;
  ActiveElementIterator<2, DOW> 
    the_ele = beginActiveElement(),
    end_ele = endActiveElement();
  // 为对几何体进行排序进行准备
  n_element = 0;
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<2, DOW> *& h_element = the_ele->h_element;
    /**
     * 半正则化的时候对点的index未做设置，而每个active的单元顶点必将在
     * 最后的网格中出现，所以在此将其设置为active状态。
     */
    for (i = 0;i < h_element->n_vertex;i ++) {
      tools.setGeometryActive(*(h_element->vertex[i]));
    }
    /**
     * 对于active的单元的每条边，如果该边被细分，则将其设置为inactive
     * 状态，过去此处通过检测邻居是否被细分来判定，考虑到并行的时候邻
     * 居可能在另外的分区，从而修改为直接判断边是否已经被细分。
     */
    int n_refined_edge = 0;
    for (i = 0;i < the_ele->n_boundary;i ++) { /// 边
      HGeometry<1, DOW>& side = *(h_element->boundary[i]);
      /**
       * 处理双生三角形的被剖分的边。在并行的情况，我们需要在这里处理这
       * 个情况，以应对邻居在另外的分区的情况。
       */      
      if (tools.isRefined(side)) { 
        tools.setGeometryActive(*side.child[0]->vertex[1]); /// 边中点
        tools.setGeometryActive(*side.child[0]);
        tools.setGeometryActive(*side.child[1]);
        tools.setGeometryInactive(side);
        n_refined_edge += 1;
      } else { /// 如果是正常的边
        tools.setGeometryActive(side);
      }
    }
    assert (n_refined_edge <= 1);
    assert (tools.isGeometryUsed(*h_element));
    tools.setParentInactive(*h_element);
    n_element += 1;
  }

  // 对所有几何体进行排序
  n_node = 0;
  n_side = 0;
  n_element = 0;
  the_ele = beginActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<2, DOW> *& h_element = the_ele->h_element;
    /**
     * 对顶点进行编号，已经编号的顶点将具有非负的index，从而不会满足
     * active判据，这样使得每个顶点只是被编号一次。
     */
    for (i = 0;i < the_ele->n_vertex;i ++) { /// 排顶点
      HGeometry<0, DOW>& node = *(h_element->vertex[i]);
      if (tools.isGeometryActive(node)) {
	node.index = n_node ++;
      }
    }
    /**
     * 每条边的编号可能已经被处理成为正的编号，也可能该边处于inactive
     * 状态，所以我们只是对满足active判据的边进行编号。
     */
    int n_refined_edge = 0;
    for (i = 0;i < the_ele->n_boundary;i ++) { /// 排边
      HGeometry<1, DOW>& side = *(h_element->boundary[i]);
      if (tools.isGeometryActive(side)) { /// 正常边
        side.index = n_side ++;
      } else if (tools.isGeometryInactive(side)) { /// 细分的边
        HGeometry<0, DOW>& mp = *(side.child[0]->vertex[1]);
        if (tools.isGeometryActive(mp)) { /// 先排中点
          mp.index = n_node ++;
        }
        for (j = 0;j < side.n_child;++ j) { /// 然后排两个孩子
          HGeometry<1, DOW> *& chd = side.child[j];
          if (tools.isGeometryActive(*chd)) {
            chd->index = n_side ++;
          }
        }
        n_refined_edge += 1;
      }
    }
    assert (n_refined_edge <= 1);
    assert (tools.isGeometryActive(*h_element));
    h_element->index = n_element ++; /// 排单元
    the_ele->index = h_element->index;
  }
	
  std::cerr << "\n\tbuilding the regular mesh ..." << std::flush;
  // 建立正则化的网格
  regular_mesh = new RegularMesh<2, DOW>(this);
  Assert (regular_mesh != NULL, ExcOutOfMemory());
  regular_mesh->point().resize(n_node);
  regular_mesh->geometry(0).resize(n_node);
  regular_mesh->geometry(1).resize(n_side);
  regular_mesh->geometry(2).resize(n_element);

#ifdef __SERIALIZATION__
  std::vector<std::vector<HGeometryBase*> >& h_geometry = regular_mesh->h_geometry();
  h_geometry.clear();
  h_geometry.resize(3);
  h_geometry[0].resize(n_node);
  h_geometry[1].resize(n_side);
  h_geometry[2].resize(n_element);
#endif

  int n_triangle = 0, n_twin_triangle = 0;
  the_ele = beginActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<2, DOW> * h_element = the_ele->h_element;
    // 先加入三角形的顶点
    for (i = 0;i < the_ele->n_vertex;i ++) {
      HGeometry<0, DOW>& vtx = *(h_element->vertex[i]);
      j = vtx.index;
      GeometryBM& pnt = regular_mesh->geometry(0,j);
      tools.regularize_add_node(vtx, 
                                pnt,
                                regular_mesh->point(j));
#ifdef __SERIALIZATION__
      h_geometry[0][j] = &vtx;
#endif
    }

    // 然后是三角形的三条边
    int m_refined_edge = -1;
    for (i = 0;i < the_ele->n_boundary;i ++) {
      HGeometry<1, DOW>& bnd = *(h_element->boundary[i]);
      if (tools.isGeometryInactive(bnd)) { // 如果是细分的边
        assert (m_refined_edge == -1);
        m_refined_edge = i;

        /// 先加入边的中点
        HGeometry<0, DOW>& mp = *(bnd.child[0]->vertex[1]);
        j = mp.index;
        tools.regularize_add_node(mp, 
                                  regular_mesh->geometry(0,j),
                                  regular_mesh->point(j));
#ifdef __SERIALIZATION__
        h_geometry[0][j] = &mp;
#endif

        HGeometry<1, DOW>& bnd0 = *(bnd.child[0]);
        j = bnd0.index;
        tools.regularize_add_side(bnd0,
                                  regular_mesh->geometry(1,j));
#ifdef __SERIALIZATION__
        h_geometry[1][j] = &bnd0;
#endif

        HGeometry<1, DOW>& bnd1 = *(bnd.child[1]);
        j = bnd1.index;
        tools.regularize_add_side(bnd1,
                                  regular_mesh->geometry(1,j));
#ifdef __SERIALIZATION__
        h_geometry[1][j] = &bnd1;
#endif
      } else { /// 正常的边的情况

        j = bnd.index;
        tools.regularize_add_side(bnd,
                                  regular_mesh->geometry(1,j));
#ifdef __SERIALIZATION__
        h_geometry[1][j] = &bnd;
#endif
      }
    }

    // 最后是三角形自身，分为三角形和双生三角形两种情况
    j = h_element->index;
    assert (j == the_ele->index);
    GeometryBM& element = regular_mesh->geometry(2,j);
    Assert(element.index() == -1, ExcInternalError());
#ifdef __SERIALIZATION__
    h_geometry[2][j] = h_element;
#endif
    element.boundaryMark() = h_element->bmark;
    if (m_refined_edge == -1) { /// 正常三角形
      n_triangle ++;
      tools.regularize_add_triangle(*h_element, element);
    } else { // 双生三角形的情形
      n_twin_triangle ++;
      tools.regularize_add_twin_triangle(*h_element, element, m_refined_edge);
    }
  }
	
  geometry_tree->unlock(); /// 对几何遗传树解锁
  std::cerr << "\n\tnodes: " << n_node 
	    << "; sides: " << n_side
	    << "; elements: " << n_element
	    << " (" << n_triangle
	    << ", " << n_twin_triangle
	    << ")" << std::endl;
	
  if (renumerate) renumerateElement();
}



template <>
void IrregularMesh<3, DOW>::regularize(bool renumerate)
{
  if (regular_mesh != NULL)
    delete regular_mesh;
  std::cerr << "Generating regular mesh from the semiregular mesh ..." << std::flush; 
  int ii[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  Tools tools;

  std::cerr << "\n\tpreparing data ..." << std::flush;
  ActiveElementIterator<3, DOW> 
    the_ele = beginActiveElement(),
    end_ele = endActiveElement();
  // 为所有几何体的排序进行准备
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<3, DOW> * h_element = the_ele->h_element;
    for (int i = 0;i < the_ele->n_vertex;++ i) {
      tools.setGeometryActive(*(h_element->vertex[i]));
    }
    for (int i = 0;i < h_element->n_boundary;++ i) { /// 面
      HGeometry<2, DOW> * surface = h_element->boundary[i];
      /**
       * 通过直接检查面是不是被细分来确定该面是否出现在最后的网格中，
       * 这和原先版本的检查邻居是不同的，这样可以使得程序在分区并行状
       * 态下有正确的行为。否则，如果邻居在其它分区上，此处将有错误的
       * 判断。下面对边的判断也是同样的方式。
       */
      if (tools.isRefined(*surface)) { /// 细分过的面
        /// 将其第三个孩子的三条边设置为active
        for (int j = 0;j < surface->child[3]->n_boundary;++ j) {
          assert ((! tools.isRefined(*surface->child[3]->boundary[j])));
          tools.setGeometryActive(*surface->child[3]->boundary[j]);
        }

        /// 将其所有孩子设置为active
        for (int j = 0;j < surface->n_child;++ j) {
          assert ((! tools.isRefined(*surface->child[j])));
          tools.setGeometryActive(*surface->child[j]);
        }

        /// 将自身设置为不是active
        tools.setGeometryInactive(*surface);
      } else { /// 正常面
        tools.setGeometryActive(*surface);
      }
      int n_refined_edge = 0;
      for (int j = 0;j < surface->n_boundary;++ j) { /// 边
        HGeometry<1, DOW> * edge = surface->boundary[j];
        if (tools.isRefined(*edge)) {
          tools.setGeometryActive(*edge->child[0]->vertex[1]);

          for (int k = 0;k < edge->n_child;++ k) {
            assert ((! tools.isRefined(*edge->child[k])));
            tools.setGeometryActive(*edge->child[k]);
          }

          tools.setGeometryInactive(*edge);
          n_refined_edge += 1;
        } else {
          tools.setGeometryActive(*edge);
        }
      }
      if (!( (n_refined_edge == 0) ||
             (n_refined_edge == 1) ||
             ((n_refined_edge == 3) &&
              tools.isRefined(*surface)))) {
        abort();
      }
    }
    Assert(tools.isGeometryActive(*h_element), ExcInternalError());
    if (! tools.isGeometryActive(*h_element)) abort();
    tools.setParentInactive(*h_element);
  }

  // 对所有几何体进行排序
  int n_node = 0;
  int n_edge = 0;
  int n_surface = 0;
  int n_element = 0;
  the_ele = beginActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    HGeometry<3, DOW> * h_element = the_ele->h_element;
    // 先排顶点
    for (int i = 0;i < the_ele->n_vertex;++ i) {
      HGeometry<0, DOW> * node = h_element->vertex[i];
      if (tools.isGeometryActive(*node)) {
	node->index = n_node ++;
      }
    }

    // 再排面和棱
    for (int i = 0;i < the_ele->n_boundary;++ i) {
      HGeometry<2, DOW> * surface = h_element->boundary[i];
      if (tools.isGeometryActive(*surface)) {
        surface->index = n_surface ++;

        for (int j = 0;j < surface->n_boundary;++ j) {
          HGeometry<1, DOW> * edge = surface->boundary[j];
          if (tools.isGeometryActive(*edge)) {
            edge->index = n_edge ++;
          } else if (tools.isGeometryInactive(*edge)) {
            /**
             * 这一段代码是为了并行的情况加上的。在串行的情形下，肯定会
             * 有一个普通三角形会将 edge 的孩子作为一条边，但是并行的时
             * 候，这个普通三角形可能位于其它分区上了，可能并不一定能够
             * 有一个普通三角形会将 edge 的孩子作为一条边。
             *
             * 在这样的情况下，我们加入这条边的两个孩子和中点。
             */
            for (int k = 0;k < edge->n_child;++ k) {
              HGeometry<1, DOW> * chd = edge->child[k];
              if (tools.isGeometryActive(*chd)) {
                chd->index = n_edge ++; /// 孩子
              }
            }
            HGeometry<0,DOW> * vtx = edge->child[0]->vertex[1];
            if (tools.isGeometryActive(*vtx)) {
              vtx->index = n_node ++; /// 中点
            }
          }
        }
      } else if (tools.isGeometryInactive(*surface)) {
        for (int j = 0;j < surface->n_child;++ j) {
          HGeometry<2,DOW> * chd = surface->child[j];
          if (tools.isGeometryActive(*chd)) {
            chd->index = n_surface ++;
          }
          for (int k = 0;k < chd->n_vertex;++ k) {
            HGeometry<0,DOW> * vtx = chd->vertex[k];
            if (tools.isGeometryActive(*vtx)) {
              vtx->index = n_node ++;
            }
          }
          for (int k = 0;k < chd->n_boundary;++ k) {
            HGeometry<1,DOW> * bnd = chd->boundary[k];
            if (tools.isGeometryActive(*bnd)) {
              bnd->index = n_edge ++;
            }
          }
        }
      }
    }
    Assert(tools.isGeometryActive(*h_element),  ExcInternalError());
    if (! tools.isGeometryActive(*h_element)) abort();
    h_element->index = n_element ++;
    the_ele->index = h_element->index;
  }
	
  std::cerr << "\n\tbuilding the regular mesh ..." << std::flush;
  // 建立正则化网格
  regular_mesh = new RegularMesh<3, DOW>(this);
  regular_mesh->point().resize(n_node);
  regular_mesh->geometry(0).resize(n_node);
  regular_mesh->geometry(1).resize(n_edge);
  regular_mesh->geometry(2).resize(n_surface);
  regular_mesh->geometry(3).resize(n_element);
	
#ifdef __SERIALIZATION__
  std::vector<std::vector<HGeometryBase*> >& h_geometry = regular_mesh->h_geometry();
  h_geometry.clear();
  h_geometry.resize(4);
  h_geometry[0].resize(n_node);
  h_geometry[1].resize(n_edge);
  h_geometry[2].resize(n_surface);
  h_geometry[3].resize(n_element);
#endif

  int n_tetrahedron = 0;
  int n_twin_tetrahedron = 0;
  int n_four_tetrahedron = 0;
  the_ele = beginActiveElement();
  for (;the_ele != end_ele;++ the_ele) {

    HGeometry<3, DOW> * h_element = the_ele->h_element;
    // 先加入顶点
    for (int i = 0;i < h_element->n_vertex;++ i) {
      HGeometry<0,DOW> * vtx = h_element->vertex[i];
      int j = vtx->index;
      tools.regularize_add_node(*vtx, 
                                regular_mesh->geometry(0, j),
                                regular_mesh->point(j));
#ifdef __SERIALIZATION__
      h_geometry[0][j] = vtx;
#endif
    }

    // 然后处理单元的表面和棱
    int n_twin_triangle_surface = 0, twin_triangle_surface[3];
    int n_nonactive_neighbour = 0, nonactive_neighbour;
    HGeometry<1, DOW> * refined_edge = NULL;
    for (int i = 0;i < the_ele->n_boundary;++ i) {
      HGeometry<2, DOW> * bnd = h_element->boundary[i];
      if (tools.isGeometryIndexed(*bnd)) { // 是一个active面
	/**
         * 检测此三角形是一个正常三角形还是双生三角形，如果是双生三角
         * 形，就将细分的那个边记在 refined_edge 中。并且此检测循环跳
         * 出的时候，k的值比 HGeometry<2, DOW>::n_boundary 要小。如果是
         * 正常三角形，则 k==HGeometry<2, DOW>::n_boundary。
         */
        int k = 0;
	for (;k < bnd->n_boundary;k ++) {
          HGeometry<1,DOW> * edge = bnd->boundary[k];
	  if (tools.isGeometryInactive(*edge)) {
            refined_edge = edge; break;
	  }
	}
        if (k < bnd->n_boundary) {
          twin_triangle_surface[n_twin_triangle_surface ++] = i;
        }

        int j = bnd->index;
        GeometryBM& surface = regular_mesh->geometry(2,j);
	if (surface.index() != -1) continue; // 如果这个面处理过了，就记一下数

	if (k == bnd->n_boundary) { // 如果是三角形
	  for (int l = 0;l < bnd->n_boundary;l ++) { // 先加三条边
	    HGeometry<1,DOW> * edge = bnd->boundary[l];
	    int m = edge->index;
            tools.regularize_add_side(*edge, regular_mesh->geometry(1,m));
#ifdef __SERIALIZATION__
            h_geometry[1][m] = edge;
#endif
	  }
	  // 再加这个三角形自身
          tools.regularize_add_triangle(*bnd, surface);
#ifdef __SERIALIZATION__
          h_geometry[2][j] = bnd;
#endif
	} else { // 是双生三角形，且是第 k 条边被细分了
          /**
           * 加入这个双生三角形两个没有细分的边和细分的边的两个孩子。
           */
          for (int m = 1;m < 3;++ m) { /// 加入两个没有细分的边
            HGeometry<1,DOW> * edge = bnd->boundary[ii[k+m]];
            int n = edge->index;
            tools.regularize_add_side(*edge, regular_mesh->geometry(1,n));
#ifdef __SERIALIZATION__
            h_geometry[1][n] = edge;
#endif
          }

          { /// 加入细分的边的两个孩子
            for (int m = 0;m < refined_edge->n_child;++ m) {
              HGeometry<1,DOW> * edge = refined_edge->child[m];
              int n = edge->index;
              tools.regularize_add_side(*edge, regular_mesh->geometry(1,n));
#ifdef __SERIALIZATION__
              h_geometry[1][n] = edge;
#endif
            }
          }

          { /// 加入细分的边的中点
            HGeometry<0,DOW> * vtx = refined_edge->child[0]->vertex[1];
            int n = vtx->index;
            tools.regularize_add_node(*vtx,
                                      regular_mesh->geometry(0, n),
                                      regular_mesh->point(n));
#ifdef __SERIALIZATION__
            h_geometry[0][n] = vtx;
#endif
          }

	  // 然后加入这个双生三角形自身
          tools.regularize_add_twin_triangle(*bnd, surface, k);
#ifdef __SERIALIZATION__
          h_geometry[2][j] = bnd;
#endif
	}
      } else { /**
                * 否则就是一个细分过的三角形，我们将细分邻居计数加一，
                * 尽管不一定是这样(细分的邻居可能在其它分区上)。
                */
	n_nonactive_neighbour ++;
	Assert(n_nonactive_neighbour == 1, ExcInternalError());
	nonactive_neighbour = i; // 记下数就可以了

        for (int j = 0;j < bnd->n_child;++ j) {
          HGeometry<2,DOW> * chd = bnd->child[j];
          {
            int n = chd->index;
            tools.regularize_add_triangle(*chd, regular_mesh->geometry(2,n));
#ifdef __SERIALIZATION__
            h_geometry[2][n] = chd;
#endif
          }

          for (int k = 0;k < chd->n_vertex;++ k) {
            HGeometry<0,DOW> * vtx = chd->vertex[k];
            int n = vtx->index;
            tools.regularize_add_node(*vtx,
                                      regular_mesh->geometry(0, n),
                                      regular_mesh->point(n));
#ifdef __SERIALIZATION__
            h_geometry[0][n] = vtx;
#endif
          }

          for (int k = 0;k < chd->n_boundary;++ k) {
            HGeometry<1,DOW> * edge = chd->boundary[k];
            int n = edge->index;
            tools.regularize_add_side(*edge,
                                      regular_mesh->geometry(1,n));
#ifdef __SERIALIZATION__
            h_geometry[1][n] = edge;
#endif
          }
        }
        
      }
    }

    // 最后加入四面体自身
    int j = h_element->index;
    GeometryBM& tetra = regular_mesh->geometry(3,j);
    tetra.index() = j;
#ifdef __SERIALIZATION__
    h_geometry[3][j] = h_element;
#endif
    if (n_nonactive_neighbour == 1) { // 四胞胎的情形
      Assert((n_twin_triangle_surface == 3), ExcInternalError());

      int tt[4][4] = {{0, 1, 2, 3}, {1, 2, 0, 3}, {2, 3, 0, 1}, {3, 0, 2, 1}};
      int * t = &tt[nonactive_neighbour][0]; /// 四面体的姿态
      HGeometry<2,DOW> *& refined_face = h_element->boundary[t[0]];
      int q[3]; /// 细分的面的姿态
      for (int k = 0;k < 3;++ k) {
        for (int l = 0;l < 3;++ l) {
          if (h_element->vertex[t[k+1]] == refined_face->vertex[l]) {
            q[k] = l; break;
          }
        }
      }

      tetra.vertex().resize(7);
      for (int k = 0;k < 4;++ k) {
	tetra.vertex(k) = h_element->vertex[t[k]]->index;
      }
      for (int k = 0;k < 3;++ k) {
        tetra.vertex(k+4) = refined_face->boundary[q[k]]->child[0]->vertex[1]->index;
      }

      tetra.boundary().resize(7);
      tetra.boundary(0) = refined_face->child[3]->index;
      for (int k = 0;k < 3;++ k) {
	tetra.boundary(k+1) = refined_face->child[q[k]]->index;
      }
      for (int k = 0;k < 3;++ k) {
	tetra.boundary(k+4) = h_element->boundary[t[k+1]]->index;
      }

      tetra.boundaryMark() = h_element->bmark;
      n_four_tetrahedron ++;
    } else if (n_twin_triangle_surface > 0) { // 双胞胎的情形
      Assert (n_twin_triangle_surface == 2, ExcInternalError());

      // 先判断是哪条边细分了
      int rei = (1<<twin_triangle_surface[0]) + (1<<twin_triangle_surface[1]);
      int tt[16] = {
        -1, -1, -1,  3, -1,  4,  2, -1,
        -1,  5,  1, -1,  0, -1, -1, -1,
      };
      int refined_edge_idx = tt[rei];
      Assert ((refined_edge_idx >= 0), ExcInternalError());

      tetra.vertex().resize(5);
      tetra.boundary().resize(4);
      int kk[6][4] = {
        {3, 1, 0, 2}, {1, 2, 0, 3}, {1, 0, 3, 2},
        {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 1, 2, 3}
      };
      for (int k = 0;k < 4;++ k) {
        int k0 = kk[refined_edge_idx][k];
        int k1 = (k > 1)?(k + 1):k;
        tetra.vertex(k1) = h_element->vertex[k0]->index;
        tetra.boundary(k) = h_element->boundary[k0]->index;
      }
      tetra.vertex(2) = refined_edge->child[0]->vertex[1]->index;

      tetra.boundaryMark() = h_element->bmark;
      n_twin_tetrahedron ++;	
    } else { // 最后是正常的四面体情形
      Assert((n_nonactive_neighbour == 0 && n_twin_triangle_surface == 0), ExcInternalError());
		
      tetra.vertex().resize(4);
      tetra.boundary().resize(4);
      for (int k = 0; k < 4;++ k) {
        tetra.vertex(k) = h_element->vertex[k]->index;
        tetra.boundary(k) = h_element->boundary[k]->index;
      }

      tetra.boundaryMark() = h_element->bmark;
      n_tetrahedron ++;
    }
  }
	
  geometry_tree->unlock(); /// 对几何遗传树解锁
  std::cerr << "\n\tnodes: " << n_node 
	    << "; edges: " << n_edge
	    << "; surface: " << n_surface
	    << "; elements: " << n_element
	    << " (" << n_tetrahedron
	    << ", " << n_twin_tetrahedron
	    << ", " << n_four_tetrahedron
	    << ")" << std::endl;
	
  if (renumerate) renumerateElement();
}


template <>
void RegularMesh<2, DOW>::writeEasyMesh(const std::string& filename) const
{
  unsigned int i, j, k, ii[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  unsigned int n_node = n_geometry(0);
  unsigned int n_side = n_geometry(1);
  unsigned int n_element = n_geometry(2);
  unsigned int n_triangle = 0;
  unsigned int n_twin_triangle = 0;
  std::cerr << "Writing easymesh data file ..." << std::endl;	
  std::cerr << "  preparing data ..." << std::flush;
  std::vector<int> index(n_element, -1);
  for (i = 0;i < n_element;i ++) {
    j = geometry(2, i).n_vertex();
    if (j == 3)
      n_triangle ++;
    else if (j == 4)
      index[i] = n_twin_triangle ++;
    else
      Assert(false, ExcInternalError());
  }
	
  /**< 准备边的邻居单元的数据 */
  std::vector<std::vector<int> > side_neighbour
    (2, std::vector<int>(n_side, -1));
  for (i = 0;i < n_element;i ++) {
    const GeometryBM& the_ele = geometry(2, i);
    if (index[i] == -1) { // this is a triangle
      for (j = 0;j < 3;j ++) {
	k = the_ele.boundary(j);
	const GeometryBM& the_side = geometry(1, k);
	if (the_side.vertex(0) == the_ele.vertex(ii[j+1])) {
	  Assert(the_side.vertex(1) == the_ele.vertex(ii[j+2]), ExcInternalError());
	  side_neighbour[1][k] = i;
	}
	else if (the_side.vertex(0) == the_ele.vertex(ii[j+2])) {
	  Assert(the_side.vertex(1) == the_ele.vertex(ii[j+1]), ExcInternalError());
	  side_neighbour[0][k] = i;
	}
	else {
	  Assert(false, ExcInternalError()); // something must be wrong!
	}
      }
    }
    else { // this is a twin triangle
      // the 0-th side
      k = the_ele.boundary(0);
      const GeometryBM& the_side_0 = geometry(1, k);
      if (the_side_0.vertex(0) == the_ele.vertex(0)) {
	Assert(the_side_0.vertex(1) == the_ele.vertex(1), ExcInternalError());
	side_neighbour[1][k] = i;
      }
      else if (the_side_0.vertex(0) == the_ele.vertex(1)) {
	Assert(the_side_0.vertex(1) == the_ele.vertex(0), ExcInternalError());
	side_neighbour[0][k] = i;
      }
      else {
	Assert(false, ExcInternalError());
      }
      // the 1-st side
      k = the_ele.boundary(1);
      const GeometryBM& the_side_1 = geometry(1, k);
      if (the_side_1.vertex(0) == the_ele.vertex(1)) {
	Assert(the_side_1.vertex(1) == the_ele.vertex(2), ExcInternalError());
	side_neighbour[1][k] = i;
      }
      else if (the_side_1.vertex(0) == the_ele.vertex(2)) {
	Assert(the_side_1.vertex(1) == the_ele.vertex(1), ExcInternalError());
	side_neighbour[0][k] = i;
      }
      else {
	Assert(false, ExcInternalError());
      }
      // the 2-nd side
      k = the_ele.boundary(2);
      const GeometryBM& the_side_2 = geometry(1, k);
      if (the_side_2.vertex(0) == the_ele.vertex(2)) {
	Assert(the_side_2.vertex(1) == the_ele.vertex(3), ExcInternalError());
	side_neighbour[1][k] = n_element + index[i];
      }
      else if (the_side_2.vertex(0) == the_ele.vertex(3)) {
	Assert(the_side_2.vertex(1) == the_ele.vertex(2), ExcInternalError());
	side_neighbour[0][k] = n_element + index[i];
      }
      else {
	Assert(false, ExcInternalError());
      }
      // the 3th side
      k = the_ele.boundary(3);
      const GeometryBM& the_side_3 = geometry(1, k);
      if (the_side_3.vertex(0) == the_ele.vertex(3)) {
	Assert(the_side_3.vertex(1) == the_ele.vertex(0), ExcInternalError());
	side_neighbour[1][k] = n_element + index[i];
      }
      else if (the_side_3.vertex(0) == the_ele.vertex(0)) {
	Assert(the_side_3.vertex(1) == the_ele.vertex(3), ExcInternalError());
	side_neighbour[0][k] = n_element + index[i];
      }
      else {
	Assert(false, ExcInternalError());
      }
    }
  }

  /**< 准备单元的邻居的数据 */
  std::vector<std::vector<int> > element_neighbour
    (3, std::vector<int>(n_element + n_twin_triangle, -1));
  for (i = 0;i < n_element;i ++) {
    const GeometryBM& the_ele = geometry(2, i);
    if (index[i] == -1) {
      for (j = 0;j < 3;j ++) {
	k = the_ele.boundary(j);
	const GeometryBM& the_side = geometry(1, k);
	if (the_side.vertex(0) == the_ele.vertex(ii[j+1])) {
	  element_neighbour[j][i] = side_neighbour[0][k];
	}
	else { // if (the_side->vertex[0] == the_ele->vertex[ii[j+2]])
	  element_neighbour[j][i] = side_neighbour[1][k];
	}
      }
    }
    else {
      element_neighbour[1][i] = n_element + index[i];
      element_neighbour[2][n_element + index[i]] = i;
      k = the_ele.boundary(0);
      const GeometryBM& the_side_0 = geometry(1, k);
      if (the_side_0.vertex(0) == the_ele.vertex(0))
	element_neighbour[2][i] = side_neighbour[0][k];
      else
	element_neighbour[2][i] = side_neighbour[1][k];
      k = the_ele.boundary(1);
      const GeometryBM& the_side_1 = geometry(1, k);
      if (the_side_1.vertex(0) == the_ele.vertex(1))
	element_neighbour[0][i] = side_neighbour[0][k];
      else
	element_neighbour[0][i] = side_neighbour[1][k];
      k = the_ele.boundary(2);
      const GeometryBM& the_side_2 = geometry(1, k);
      if (the_side_2.vertex(0) == the_ele.vertex(2))
	element_neighbour[0][n_element + index[i]] = side_neighbour[0][k];
      else
	element_neighbour[0][n_element + index[i]] = side_neighbour[1][k];
      k = the_ele.boundary(3);
      const GeometryBM& the_side_3 = geometry(1, k);
      if (the_side_3.vertex(0) == the_ele.vertex(3))
	element_neighbour[1][n_element + index[i]] = side_neighbour[0][k];
      else
	element_neighbour[1][n_element + index[i]] = side_neighbour[1][k];
    }
  }
  std::cerr << " OK!" << std::endl;

  std::cerr << "  writing node data ..." << std::flush;
  std::ofstream os((filename + ".n").c_str());
  os << n_node << "\t"
     << n_element + n_twin_triangle << "\t"
     << n_side + n_twin_triangle << " **(Nnd, Nee, Nsd)**\n";
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  for (i = 0;i < n_node;i ++) {
    j = geometry(0, i).vertex(0);
    os << i << "\t"
       << point(j)[0] << "\t"
       << point(j)[1] << "\t"
       << boundaryMark(0, i) << "\n";
  }
  os << "---------------------------------------------------------------\n"
     << "   n:  x                          y                            \n";
  os.close();
  std::cerr << " OK!" << std::endl;

  std::cerr << "  writing side data ..." << std::flush;
  os.open((filename + ".s").c_str());
  os << n_side + n_twin_triangle << "\n";
  for (i = 0;i < n_side;i ++) {
    const GeometryBM& the_side = geometry(1, i);
    os << i << "\t"
       << the_side.vertex(0) << "\t"
       << the_side.vertex(1) << "\t"
       << side_neighbour[0][i] << "\t"
       << side_neighbour[1][i] << "\t"
       << boundaryMark(1, i) << "\n";
  }
  for (i = 0;i < n_element;i ++) {
    if (index[i] == -1) continue;
    os << n_side + index[i] << "\t"
       << geometry(2,i).vertex(0) << "\t"
       << geometry(2,i).vertex(2) << "\t"
       << i << "\t"
       << n_element + index[i] << "\t"
       << boundaryMark(2, i) << "\n";
  }
  os << "---------------------------------------------------------------\n"
     << "   s:    c      d       ea       eb                            \n";
  os.close();
  std::cerr << " OK!" << std::endl;

  std::cerr << "  writing element data ..." << std::flush;
  os.open((filename+".e").c_str());
  os << n_element + n_twin_triangle << "\t"
     << n_node << "\t"
     << n_side + n_twin_triangle << " **(Nee, Nnd, Nsd)**\n";
  for (i = 0;i < n_element;i ++) {
    const GeometryBM& the_ele = geometry(2, i);
    if (index[i] == -1) {
      os << i << "\t"
	 << the_ele.vertex(0) << "\t"
	 << the_ele.vertex(1) << "\t"
	 << the_ele.vertex(2) << "\t"
	 << element_neighbour[0][i] << "\t"
	 << element_neighbour[1][i] << "\t"
	 << element_neighbour[2][i] << "\t"
	 << the_ele.boundary(0) << "\t"
	 << the_ele.boundary(1) << "\t"
	 << the_ele.boundary(2) << "\n";
    }
    else {
      os << i << "\t"
	 << the_ele.vertex(0) << "\t"
	 << the_ele.vertex(1) << "\t"
	 << the_ele.vertex(2) << "\t"
	 << element_neighbour[0][i] << "\t"
	 << element_neighbour[1][i] << "\t"
	 << element_neighbour[2][i] << "\t"
	 << the_ele.boundary(1) << "\t"
	 << n_side + index[i] << "\t"
	 << the_ele.boundary(0) << "\n";
    }
  }
  for (i = 0;i < n_element;i ++) {
    if (index[i] == -1) continue;
    const GeometryBM& the_ele = geometry(2, i);
    j = n_element + index[i];
    os << j << "\t"
       << geometry(2,i).vertex(0) << "\t"
       << geometry(2,i).vertex(2) << "\t"
       << geometry(2,i).vertex(3) << "\t"
       << element_neighbour[0][j] << "\t"
       << element_neighbour[1][j] << "\t"
       << element_neighbour[2][j] << "\t"
       << the_ele.boundary(2) << "\t"
       << the_ele.boundary(3) << "\t"
       << n_side + index[i] << "\n";
  }
  os << "---------------------------------------------------------------\n"
     << "   e:  i,   j,   k,   ei,   ej,   ek,   si,   sj,   sk         \n";
  os.close();
  std::cerr << " OK!" << std::endl;
}



template <>
void RegularMesh<2, DOW>::writeTecplotData(const std::string& filename) const
{
  int i;
  std::cerr << "Write mesh data into Tecplot data file " 
	    << filename << " ... " << std::flush;
  std::ofstream os(filename.c_str());
  os << "TITLE = \x22" << "2D mesh data generated by AFEPack" << "\x22\n"
     << "VARIABLES = \x22" << "X\x22, \x22" << "Y\x22\n";
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);

  int n_node = n_point();
  int n_element = n_geometry(2);
  os << "ZONE N=" << n_node
     << ",E=" << n_element
     << ",F=FEPOINT ET=TRIANGLE\n";
  for (i = 0;i < n_node;i ++)
    os << point(i) << "\n";
  os << "\n";
  for (i = 0;i < n_element;i ++) {
    int n_vertex = geometry(2,i).n_vertex();
    switch (n_vertex) {
    case 3:
      os << geometry(0,geometry(2,i).vertex(0)).vertex(0)+1 << "\t"
	 << geometry(0,geometry(2,i).vertex(1)).vertex(0)+1 << "\t"
	 << geometry(0,geometry(2,i).vertex(2)).vertex(0)+1 << "\n";
      break;
    case 4:
      os << geometry(0,geometry(2,i).vertex(0)).vertex(0)+1 << "\t"
	 << geometry(0,geometry(2,i).vertex(1)).vertex(0)+1 << "\t"
	 << geometry(0,geometry(2,i).vertex(2)).vertex(0)+1 << "\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os.close();
  std::cerr << "OK!" << std::endl;
}



template <>
void RegularMesh<3, DOW>::writeTecplotData(const std::string& filename) const
{
  int i;
  std::cerr << "Write mesh data into Tecplot data file " 
	    << filename << " ... " << std::flush;
  std::ofstream os(filename.c_str());
  os << "TITLE = \x22" << "3D mesh data generated by AFEPack" << "\x22\n"
     << "VARIABLES = \x22" << "X\x22, \x22" << "Y\x22, \x22" << "Z\x22\n";
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);

  int n_node = n_point();
  int n_element = n_geometry(3);
  os << "ZONE N=" << n_node
     << ",E=" << n_element
     << ",F=FEPOINT ET=TETRAHEDRON\n";
  for (i = 0;i < n_node;i ++)
    os << point(i) << "\n";
  os << "\n";
  for (i = 0;i < n_element;i ++) {
    int n_vertex = geometry(3,i).n_vertex();
    switch (n_vertex) {
    case 4:
      os << geometry(0,geometry(3,i).vertex(0)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(1)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(2)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(3)).vertex(0) + 1 << "\n";
      break;
    case 5:
      os << geometry(0,geometry(3,i).vertex(0)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(1)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(3)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(4)).vertex(0) + 1 << "\n";
      break;
    case 7:
      os << geometry(0,geometry(3,i).vertex(0)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(1)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(2)).vertex(0) + 1 << "\t"
	 << geometry(0,geometry(3,i).vertex(3)).vertex(0) + 1 << "\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os.close();
  std::cerr << "OK!" << std::endl;
}

template <>
void RegularMesh<1, DOW>::writeOpenDXData(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << "object 1 class array type float rank 1 shape " << DOW + 1 << " item " 
     << 2*n_node << " data follows\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\t0.0\n"
       << point(geometry(0,i).vertex(0)) << "\t1.0\n";
  }
	
  int n_element = n_geometry(1);
  os << "\nobject 2 class array type int rank 1 shape 4 item "
     << n_element << " data follows\n";
  for (i = 0;i < n_element;i ++) {
    os << 2*geometry(1,i).vertex(0) << "\t"
       << 2*geometry(1,i).vertex(0) + 1 << "\t"
       << 2*geometry(1,i).vertex(1) << "\t"
       << 2*geometry(1,i).vertex(1) + 1 << "\n";
  }
  os << "attribute \"element type\" string \"quads\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction-1d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "end\n";
  os.close();
}

template <>
void RegularMesh<2, DOW>::writeOpenDXData(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << "object 1 class array type float rank 1 shape " << DOW << " item " 
     << n_node << " data follows\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }
	
  int n_element = n_geometry(2);
  for (i = 0, j = 0;i < n_element;i ++) {
    switch (geometry(2,i).n_vertex()) {
    case 3:
      j += 1;
      break;
    case 4:
      j += 2;
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << "\nobject 2 class array type int rank 1 shape 3 item "
     << j << " data follows\n";
  for (i = 0;i < n_element;i ++) {
    switch (geometry(2,i).n_vertex()) {
    case 3:
      os << geometry(2,i).vertex(0) << "\t"
	 << geometry(2,i).vertex(1) << "\t"
	 << geometry(2,i).vertex(2) << "\t\n";
      break;
    case 4:
      os << geometry(2,i).vertex(0) << "\t"
	 << geometry(2,i).vertex(1) << "\t"
	 << geometry(2,i).vertex(2) << "\t\n";
      os << geometry(2,i).vertex(0) << "\t"
	 << geometry(2,i).vertex(2) << "\t"
	 << geometry(2,i).vertex(3) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << "attribute \"element type\" string \"triangles\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction-2d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "end\n";
  os.close();
}



template <>
void RegularMesh<3, DOW>::writeOpenDXData(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << "object 1 class array type float rank 1 shape " << DOW << " item " 
     << n_node << " data follows\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }
	
  int n_element = n_geometry(3);
  for (i = 0, j = 0;i < n_element;i ++) {
    switch (geometry(3,i).n_vertex()) {
    case 4:
      j += 1;
      break;
    case 5:
      j += 2;
      break;
    case 7:
      j += 4;
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << "\nobject 2 class array type int rank 1 shape 4 item "
     << j << " data follows\n";
  for (i = 0;i < n_element;i ++) {
    switch (geometry(3,i).n_vertex()) {
    case 4:
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(1) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(3) << "\t\n";
      break;
    case 5:
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(1) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(4) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(3) << "\t"
	 << geometry(3,i).vertex(4) << "\t\n";
      break;
    case 7:
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(1) << "\t"
	 << geometry(3,i).vertex(6) << "\t"
	 << geometry(3,i).vertex(5) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(4) << "\t"
	 << geometry(3,i).vertex(6) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(3) << "\t"
	 << geometry(3,i).vertex(5) << "\t"
	 << geometry(3,i).vertex(4) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(4) << "\t"
	 << geometry(3,i).vertex(5) << "\t"
	 << geometry(3,i).vertex(6) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << "attribute \"element type\" string \"tetrahedra\"\n"
     << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction-3d\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "end\n";
  os.close();
}



template <>
void HGeometryTree<2, DOW>::readEasyMesh(const std::string& filename)
{
  std::cerr << "Reading easymesh data file ..." << std::endl;
	
  std::ifstream is((filename + ".n").c_str());
  int i, j, k, l;
  int n_node, n_side, n_element;
  char dummy[64];
  is >> n_node >> n_element >> n_side;
  is.getline(dummy, 64);
  std::vector<HGeometry<0, DOW> *> node(n_node);
  std::vector<HGeometry<1, DOW> *> side(n_side);
  std::vector<HGeometry<2, DOW> *> element(n_element);
  std::cerr << "\treading the nodes data ..." << std::flush;
  for (i = 0;i < n_node;i ++) {
    node[i] = new HGeometry<0, DOW>();
    Assert(node[i] != NULL, ExcOutOfMemory());
    is >> j
       >> *(dynamic_cast<Point<DOW> *>(node[i]))
       >> node[i]->bmark;
  }
  is.close();
  std::cerr << " OK!" << std::endl;
	
  is.open((filename + ".s").c_str());
  is >> i;
  Assert(i == n_side, ExcInternalError());
  std::cerr << "\treading the sides data ..." << std::flush;	
  for (i = 0;i < n_side;i ++) {
    side[i] = new HGeometry<1, DOW>();
    Assert(side[i] != NULL, ExcOutOfMemory());
    is >> l >> j >> k;
    side[i]->vertex[0] = node[j];
    side[i]->vertex[1] = node[k];
    is >> j >> k >> side[i]->bmark;
  }
  is.close();
  std::cerr << " OK!" << std::endl;

  is.open((filename + ".e").c_str());
  is >> i >> j >> k;
  Assert(i == n_element, ExcInternalError());
  Assert(j == n_node, ExcInternalError());
  Assert(k == n_side, ExcInternalError());
  is.getline(dummy, 64);
  std::cerr << "\treading the elements data ..." << std::flush;	
  for (i = 0;i < n_element;i ++) {
    element[i] = new HGeometry<2, DOW>();
    Assert(element[i] != NULL, ExcOutOfMemory());
    is >> l >> j >> k >> l;
    element[i]->vertex[0] = node[j];
    element[i]->vertex[1] = node[k];
    element[i]->vertex[2] = node[l];
    element[i]->bmark = 0;
    is >> j >> k >> l;
    is >> j >> k >> l;
    element[i]->boundary[0] = side[j];
    element[i]->boundary[1] = side[k];
    element[i]->boundary[2] = side[l];
    root_element.push_back(element[i]);
  }
  is.close();
  std::cerr << " OK!" << std::endl;
}


template <>
void RegularMesh<1, DOW>::writeSimplestSimplexMesh(const std::string& filename) const
{
  std::cout << "Not implemented." << std::endl;
}

template <>
void RegularMesh<2, DOW>::writeSimplestSimplexMesh(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << n_node << "\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }
	
  int n_element = n_geometry(2);
  for (i = 0, j = 0;i < n_element;i ++) {
    switch (geometry(2,i).n_vertex()) {
    case 3:
      j += 1;
      break;
    case 4:
      j += 2;
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << j << "\n";
  for (i = 0;i < n_element;i ++) {
    switch (geometry(2,i).n_vertex()) {
    case 3:
      os << geometry(2,i).vertex(0) << "\t"
	 << geometry(2,i).vertex(1) << "\t"
	 << geometry(2,i).vertex(2) << "\t\n";
      break;
    case 4:
      os << geometry(2,i).vertex(0) << "\t"
	 << geometry(2,i).vertex(1) << "\t"
	 << geometry(2,i).vertex(2) << "\t\n";
      os << geometry(2,i).vertex(0) << "\t"
	 << geometry(2,i).vertex(2) << "\t"
	 << geometry(2,i).vertex(3) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os.close();
}



template <>
void RegularMesh<3, DOW>::writeSimplestSimplexMesh(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << n_node << "\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }
	
  int n_element = n_geometry(3);
  for (i = 0, j = 0;i < n_element;i ++) {
    switch (geometry(3,i).n_vertex()) {
    case 4:
      j += 1;
      break;
    case 5:
      j += 2;
      break;
    case 7:
      j += 4;
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os << j << "\n";
  for (i = 0;i < n_element;i ++) {
    switch (geometry(3,i).n_vertex()) {
    case 4:
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(1) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(3) << "\t\n";
      break;
    case 5:
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(1) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(4) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(3) << "\t"
	 << geometry(3,i).vertex(4) << "\t\n";
      break;
    case 7:
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(1) << "\t"
	 << geometry(3,i).vertex(6) << "\t"
	 << geometry(3,i).vertex(5) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(2) << "\t"
	 << geometry(3,i).vertex(4) << "\t"
	 << geometry(3,i).vertex(6) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(3) << "\t"
	 << geometry(3,i).vertex(5) << "\t"
	 << geometry(3,i).vertex(4) << "\t\n";
      os << geometry(3,i).vertex(0) << "\t"
	 << geometry(3,i).vertex(4) << "\t"
	 << geometry(3,i).vertex(5) << "\t"
	 << geometry(3,i).vertex(6) << "\t\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os.close();
}

template <>
void RegularMesh<1, DOW>::writeSimplexMesh(const std::string& filename) const
{
  std::cout << "Not implemented." << std::endl;
}

template <>
void RegularMesh<2, DOW>::writeSimplexMesh(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << n_node << "\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }

  os << "\n" << n_node << "\n";
  for (i = 0;i < n_node;i ++) {
    os << i << "\n"
       << "1\t" << i << "\n"
       << "1\t" << i << "\n"
       << geometry(0,i).boundaryMark() << "\n";
  }

  int n_twin_triangle = 0;
  int n_element = n_geometry(2);
  for (i = 0;i < n_element;i ++) {
    if (geometry(2,i).n_vertex() == 4) {
      n_twin_triangle += 1;
    }
  }

  int n_side = n_geometry(1);
  os << "\n" << n_side + n_twin_triangle << "\n";
  for (i = 0;i < n_side;i ++) {
    os << i << "\n"
       << "2\t" << geometry(1,i).vertex(0) << " " << geometry(1,i).vertex(1) << "\n"
       << "2\t" << geometry(1,i).vertex(0) << " " << geometry(1,i).vertex(1) << "\n"
       << geometry(1,i).boundaryMark() << "\n";  
  }
  for (i = 0, j = 0;i < n_element;i ++) {
    if (geometry(2,i).n_vertex() == 4) {
      os << n_side + j << "\n"
         << "2\t" << geometry(2,i).vertex(0) << " " << geometry(2,i).vertex(2) << "\n"
         << "2\t" << geometry(2,i).vertex(0) << " " << geometry(2,i).vertex(2) << "\n"
         << geometry(2,i).boundaryMark() << "\n";
      j += 1;
    }
  }

  os << "n" << n_element + n_twin_triangle << "\n";
  for (i = 0, j = 0;i < n_element;i ++) {
    const GeometryBM& geo = geometry(2,i);
    switch (geo.n_vertex()) {
    case 3:
      os << i + j << "\n"
         << "3\t" << geo.vertex(0) << " " << geo.vertex(1) << " " << geo.vertex(2) << "\n"
         << "3\t" << geo.boundary(0) << " " << geo.boundary(1) << " " << geo.boundary(2) << "\n"
         << geo.boundaryMark() << "\n";
      break;
    case 4:
      os << i + j << "\n"
         << "3\t" << geo.vertex(0) << " " << geo.vertex(1) << " " << geo.vertex(2) << "\n"
         << "3\t" << geo.boundary(1) << " " << n_side + j << " " << geo.boundary(0) << "\n"
         << geo.boundaryMark() << "\n";
      j += 1;
      os << i + j << "\n"
         << "3\t" << geo.vertex(0) << " " << geo.vertex(2) << " " << geo.vertex(3) << "\n"
         << "3\t" << geo.boundary(2) << " " << geo.boundary(3) << " " << n_side + j << "\n"
         << geo.boundaryMark() << "\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  os.close();
}



template <>
void RegularMesh<3, DOW>::writeSimplexMesh(const std::string& filename) const
{
  std::ofstream os(filename.c_str());
	
  os.precision(8);
  os.setf(std::ios::scientific, std::ios::floatfield);
  int n_node = n_point();
  int i, j;
	
  os << n_node << "\n";
  for (i = 0;i < n_node;i ++) {
    os << point(geometry(0,i).vertex(0)) << "\n";
  }

  os << "\n" << n_node << "\n";
  for (i = 0;i < n_node;i ++) {
    os << i << "\n"
       << "1\t" << i << "\n"
       << "1\t" << i << "\n"
       << geometry(0,i).boundaryMark() << "\n";
  }

  int n_twin_triangle = 0;
  int n_side = n_geometry(1);
  int n_face = n_geometry(2);
  for (i = 0;i < n_face;++ i) {
    if (geometry(2,i).n_vertex() == 4) {
      n_twin_triangle += 1;
    }
  }

  int n_twin_tetrahedron = 0;
  int n_four_tetrahedron = 0;
  int n_element = n_geometry(3);
  for (i = 0, j = 0;i < n_element;i ++) {
    switch (geometry(3,i).n_vertex()) {
    case 5:
      n_twin_tetrahedron += 1;
      break;
    case 7:
      n_four_tetrahedron += 1;
      break;
    }
  }

  os << "\n" << n_side + n_twin_triangle << "\n";
  for (i = 0;i < n_side;i ++) {
    os << i << "\n"
       << "2\t" << geometry(1,i).vertex(0) << " " << geometry(1,i).vertex(1) << "\n"
       << "2\t" << geometry(1,i).vertex(0) << " " << geometry(1,i).vertex(1) << "\n"
       << geometry(1,i).boundaryMark() << "\n";  
  }
  for (i = 0, j = 0;i < n_face;i ++) {
    if (geometry(2,i).n_vertex() == 4) {
      os << n_side + j << "\n"
         << "2\t" << geometry(2,i).vertex(0) << " " << geometry(2,i).vertex(2) << "\n"
         << "2\t" << geometry(2,i).vertex(0) << " " << geometry(2,i).vertex(2) << "\n"
         << geometry(2,i).boundaryMark() << "\n";
      j += 1;
    }
  }
    
  os << "n" << (n_face + 
                n_twin_triangle + 
                n_twin_tetrahedron +
                3*n_four_tetrahedron) << "\n";
  for (i = 0;i < n_face;i ++) {
    const GeometryBM& geo = geometry(2,i);
    switch (geo.n_vertex()) {
    case 3:
      os << i + j << "\n"
         << "3\t" << geo.vertex(0) << " " << geo.vertex(1) << " " << geo.vertex(2) << "\n"
         << "3\t" << geo.boundary(0) << " " << geo.boundary(1) << " " << geo.boundary(2) << "\n"
         << geo.boundaryMark() << "\n";
      break;
    case 4:
      os << i + j << "\n"
         << "3\t" << geo.vertex(0) << " " << geo.vertex(1) << " " << geo.vertex(2) << "\n"
         << "3\t" << geo.boundary(1) << " " << n_side + j << " " << geo.boundary(0) << "\n"
         << geo.boundaryMark() << "\n";
      break;
    default:
      Assert(false, ExcInternalError());
    }
  }
  for (i = 0, j = 0;i < n_face;i ++) {
    const GeometryBM& geo = geometry(2,i);
    if (geo.n_vertex() == 3) continue;
    os << n_face + j << "\n"
       << "3\t" << geo.vertex(0) << " " << geo.vertex(2) << " " << geo.vertex(3) << "\n"
       << "3\t" << geo.boundary(2) << " " << geo.boundary(3) << " " << n_side + j << "\n"
       << geo.boundaryMark() << "\n";
    j += 1;
  }
  for (i = 0, j = 0;i < n_element;i ++) {
    const GeometryBM& geo = geometry(3,i);
    if (geo.n_vertex() != 5) continue;
    os << n_face + n_twin_triangle + j << "\n"
       << "3\t" << geo.vertex(0) << " " << geo.vertex(2) << " " << geo.vertex(4) << "\n"
       << "3\t" << n_side + geo.boundary(0) << " "
              //<< geo.???
                << n_side + geo.boundary(3) << "\n"
       << geo.boundaryMark() << "\n";
    j += 1;
  }

  /// Todo: ...

  os.close();
}

/**
 * end of file
 * 
 */


