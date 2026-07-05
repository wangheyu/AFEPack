/**
 * @file   GmshMesh.cpp
 * @author Robert Lie and Guanghui Hu
 * @date   Sat Oct  6 11:48:47 2007
 * 
 * @brief  本程序作用: 读取 Gmsh 的 网格文件 "*.msh", 并将其转化成 AFEPack 适用的网格文件 "*.mesh".
 *         注意: Gmsh 有两种文件格式, Version 1.0 和 Version 2.0. 本程序适用于 Version 2.0.
 *              两种格式的区别参见 Gmsh 的 Reference manual. 该文档可以从 http://geuz.org/gmsh 下载.
 * 
 * 
 */

#include <AFEPack/GmshMesh.h>

AFEPACK_OPEN_NAMESPACE

#define N_POINT_NODE		1
#define N_LINE_NODE		2
#define N_TRIANGLE_NODE		3
#define N_QUADRANGLE_NODE	4
#define N_TETRAHEDRON_NODE	4
#define N_HEXAHEDRON_NODE	8
#define N_PRISM_NODE		6
#define N_PYRAMID_NODE		5


GmshMesh::GmshMesh()
{
  setTemplateGeometry(*(new std::vector<TemplateGeometry<3> >(1)));
  templateGeometry(0).readData("tetrahedron.tmp_geo");
  //   templateGeometry(1).readData("hexahedron.tmp_geo");
  //   templateGeometry(2).readData("prism.tmp_geo");
  //   templateGeometry(3).readData("pyramid.tmp_geo");
};

GmshMesh::~GmshMesh()
{
  delete &(templateGeometry());
};



#if 1 /// for Gmsh data version 4.1

/*
  1. In .msh file, the counting for everything is from 1. In .mesh, it is from 0.
  2. 

*/

void GmshMesh::readData(const std::string& filename)
{
  std::cerr << "Reading in gmsh data ..." << std::endl;

  std::string dummy; /// for all useless information  
  double gmsh_version = 0;
  int gmsh_data_format = 0;
  int gmsh_data_precision = 0;

  std::ifstream is(filename.c_str());
  is >> dummy;/// $MeshFormat
  is >> gmsh_version;/// integer : for example, 2 means the version of Gmsh is 2.0
  is >> gmsh_data_format;/// integer : for example, 0 means the type of data file is ASCII
  is >> gmsh_data_precision;/// integer : the size of the floating point numbers used in the file.
  is >> dummy;/// $EndMeshFormat$
  
  GmshMesh::GeometryType gtype;


  if(int(gmsh_version) == 4){/// if this is 4.x data format
  
    int n_0D_entity, n_1D_entity, n_2D_entity, n_3D_entity;
    int n_bounding_points = 0, n_tags=0;
    
    /// for reading entities
    std::cerr << "\tReading entity data ..." << std::endl;  
    is >> dummy;/// $Entities
    is >> n_0D_entity >> n_1D_entity >> n_2D_entity >> n_3D_entity;

    std::vector<std::vector<int> > tags_0D_entity(n_0D_entity, std::vector<int>(2, 0));
    for(int i = 0;i < n_0D_entity;++ i){
      is >> tags_0D_entity[i][0] >> dummy >> dummy >> dummy >> n_tags;
      if(n_tags != 0){
	is >> tags_0D_entity[i][1];/// the first tag is boundary mark
	for(int j = 1;j < n_tags;++ j){
	  is >> dummy;
	}
      }
    }

    std::vector<std::vector<int> > tags_1D_entity(n_1D_entity, std::vector<int>(2, 0));
    for(int i = 0;i < n_1D_entity;++ i){
      is >> tags_1D_entity[i][0] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> n_tags;
      if(n_tags != 0){
	is >> tags_1D_entity[i][1];/// the first tag is boundary mark
	for(int j = 1;j < n_tags;++ j){
	  is >> dummy;
	}
      }

      is >> n_bounding_points; /// no. of bounding points
      for(int j = 0;j < n_bounding_points;++ j){
	is >> dummy;/// tags of those bounding points. Currenlty we do not need it.
      }
    }

    std::vector<std::vector<int> > tags_2D_entity(n_2D_entity, std::vector<int>(2, 0));
    for(int i = 0;i < n_2D_entity;++ i){
      is >> tags_2D_entity[i][0] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> n_tags;
      if(n_tags != 0){
	is >> tags_2D_entity[i][1];/// the first tag is boundary mark
	for(int j = 1;j < n_tags;++ j){
	  is >> dummy;
	}
      }

      is >> n_bounding_points; /// no. of bounding points
      for(int j = 0;j < n_bounding_points;++ j){
	is >> dummy;
      }
    }

    std::vector<std::vector<int> > tags_3D_entity(n_3D_entity, std::vector<int>(2, 0));
    for(int i = 0;i < n_3D_entity;++ i){
      is >> tags_3D_entity[i][0] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> n_tags;
      if(n_tags != 0){
	is >> tags_3D_entity[i][1];/// the first tag is boundary mark
	for(int j = 1;j < n_tags;++ j){
	  is >> dummy;
	}
      }

      is >> n_bounding_points; /// no. of bounding points
      for(int j = 0;j < n_bounding_points;++ j){
	is >> dummy;
      }
    }
  
    is >> dummy;/// $EndEntities
    is >> dummy;/// $Nodes
    /// reading all nodes information from msh file

    int n_entity_blocks = 0, n_nodes = 0, min_node_tag = 0, max_node_tag = 0;
    is >> n_entity_blocks >> n_nodes >> min_node_tag >> max_node_tag;

    assert(n_nodes == max_node_tag);///

    pnt.resize(n_nodes);
    std::vector<int> ind(n_nodes + 1), tmp_ind(n_nodes);

    int idx_node = 0;
    int idx_coord = 0;
    int n_nodes_inBlock = 0;  
    for(int i = 0;i < n_entity_blocks;++ i){
      is >> dummy >> dummy >> dummy >> n_nodes_inBlock;
      for(int j = 0;j < n_nodes_inBlock;++ j){
	is >> tmp_ind[idx_node];
	idx_node ++;
      }
      for(int j = 0;j < n_nodes_inBlock;++ j){
	is >> pnt[idx_coord][0] >> pnt[idx_coord][1] >> pnt[idx_coord][2];
	idx_coord ++;
      }
    }
    assert(idx_node == n_nodes);
  
    /// to build an index record in vector ind.
    for(int i = 0;i < n_nodes;++ i){
      ind[tmp_ind[i]] = i;
    }
  
    is >> dummy;/// $EndNodes
    is >> dummy;/// $Elements

    int n_elements = 0, min_element_tag = 0, max_element_tag = 0;

    is >> n_entity_blocks >> n_elements >> min_element_tag >> max_element_tag;
    assert(n_elements == max_element_tag);

    boundary_ele.resize(n_elements);/// to save all boundary elements
				    /// with boundary mark
				    /// (SimplestGeometryBM)

    ele.resize(n_elements);

    int idx_tet = 0, idx_bnd_ele = 0;
    int dim = 0, tag_entity = 0, type_element = 0, n_records = 0;
    int n_vtx = 0, idx_vtx = 0;
    /// this (dim + tag_entity) would correspond to the boundary mark given in tags_[0123]D_entity
  
    for(int i = 0;i < n_entity_blocks;++ i){

      is >> dim >> tag_entity >> type_element >> n_records;

      if (dim == 0){
	for(int j = 0;j < n_records;++ j){
	  is >> dummy;/// element_tag
	  if(type_element == POINT){
	    n_vtx = N_POINT_NODE;
	    boundary_ele[idx_bnd_ele].vertex.resize(n_vtx);
	    for(int l = 0;l < n_vtx;++ l){
	      is >> idx_vtx;
	      boundary_ele[idx_bnd_ele].vertex[l] = ind[idx_vtx];
	    }
	    boundary_ele[idx_bnd_ele].dim = dim;

	    /// for boundaryMark,  std::vector<std::vector<double> > tags_0D_entity;
	    //int n_0D_entity = tags_0D_entity.size();
	    for(int k = 0;k < n_0D_entity;++ k){
	      if(tags_0D_entity[k][0] == tag_entity){
		boundary_ele[idx_bnd_ele].boundaryMark = tags_0D_entity[k][1];
		break;/// jump out this for loop
	      }
	    }
	  }
	  idx_bnd_ele ++;
	}
      }
      else if(dim == 1){
	for(int j = 0;j < n_records;++ j){
	  is >> dummy;
	  if(type_element == LINE){
	    n_vtx = N_LINE_NODE;
	    boundary_ele[idx_bnd_ele].vertex.resize(n_vtx);
	    for(int l = 0;l < n_vtx;++ l){
	      is >> idx_vtx;
	      boundary_ele[idx_bnd_ele].vertex[l] = ind[idx_vtx];
	    }
	    boundary_ele[idx_bnd_ele].dim = dim;

	    /// for boundaryMark,  std::vector<std::vector<double> > tags_0D_entity;
	    //int n_1D_entity = tags_1D_entity.size();
	    for(int k = 0;k < n_1D_entity;++ k){
	      if(tags_1D_entity[k][0] == tag_entity){
		boundary_ele[idx_bnd_ele].boundaryMark = tags_1D_entity[k][1];
		break;/// jump out this for loop
	      }
	    }
	  }
	  idx_bnd_ele ++;	  
	}
      }
      else if(dim == 2){
	for(int j = 0;j < n_records;++ j){
	  is >> dummy;
	  if(type_element == TRIANGLE){
	    n_vtx = N_TRIANGLE_NODE;
	    boundary_ele[idx_bnd_ele].vertex.resize(n_vtx);
	    for(int l = 0;l < n_vtx;++ l){
	      is >> idx_vtx;
	      boundary_ele[idx_bnd_ele].vertex[l] = ind[idx_vtx];
	    }
	    boundary_ele[idx_bnd_ele].dim = dim;

	    /// for boundaryMark,  std::vector<std::vector<double> > tags_0D_entity;
	    //int n_2D_entity = tags_2D_entity.size();
	    for(int k = 0;k < n_2D_entity;++ k){
	      if(tags_2D_entity[k][0] == tag_entity){
		boundary_ele[idx_bnd_ele].boundaryMark = tags_2D_entity[k][1];
		break;/// jump out this for loop
	      }
	    }
	  }
	  idx_bnd_ele ++;	  	  
	}
      }
      else if(dim == 3){ /// we need this information, for tetrahedron
	for(int j = 0;j < n_records;++ j){
	  is >> dummy;/// tag of element
	  if(type_element == TETRAHEDRON){
	    n_vtx = N_TETRAHEDRON_NODE;
	    ele[idx_tet].vertex.resize(n_vtx);
	    for(int l = 0;l < n_vtx;++ l){
	      //int idx_vtx = 0;
	      is >> idx_vtx;
	      ele[idx_tet].vertex[l] = ind[idx_vtx];	      
	    }
	    /**
	     * Maybe we need to exchange the order of any two pnts, to
	     * make the volume be Positive.
	     * 
	     */
	    // int l = ele[idx_tet].vertex[1]; 
	    // ele[idx_tet].vertex[1] = ele[idx_tet].vertex[0];
	    // ele[idx_tet].vertex[0] = l;
	    // ele[idx_tet].template_element = 0;
	    
	  }
	  idx_tet ++;
	}
      }  
    }
    ele.resize(idx_tet);
    boundary_ele.resize(idx_bnd_ele);
    is >> dummy;/// $EndElements

  }/// if this is 4.1 version of Gmsh data
  else if(int(gmsh_version == 2)){
    std::cerr << "Reading in gmsh data ..." << std::endl;
    int i, j, k, l, m, n, n_point, n_geometry, n_tags;
    GmshMesh::GeometryType gtype;
    std::string dummy;
    std::ifstream is(filename.c_str());
    is >> dummy;/// $MeshFormat
    is >> dummy;/// integer : for example, 2 means the version of Gmsh is 2.0
    is >> dummy;/// integer : for example, 0 means the type of data file is ASCII
    is >> dummy;/// integer : the size of the floating point numbers used in the file.
    is >> dummy;/// $EndMeshFormat$
    is >> dummy;/// $Nodes
    std::cerr << "\tReading nodes data ..." << std::endl;
    is >> n_point;
    pnt.resize(n_point);
    std::vector<int> ind0(n_point);
    for (i = 0, k = 0;i < n_point;i ++) {
      is >> j;
      if (j > k) k = j;
      ind0[i] = j;
      is >> pnt[i];
    }
    std::vector<int> ind(k+1, -1);
    for (i = 0;i < n_point;i ++) ind[ind0[i]] = i;
    is >> dummy;/// $EndNodes
    std::cerr << "\tReading geometry data ..." << std::endl;
    is >> dummy;/// $Elements
    is >> n_geometry;
    ele.resize(n_geometry);
    for (i = 0, j = 0;i < n_geometry;i ++) {
      is >> k; // geometry number
      is >> gtype; // geometry type
      is >> n_tags; // number of tags;
      is >> k; // reading the tags, the first tag is the boundary mark
      for (l = 1;l < n_tags;++ l) is >> m;
      GeometryBM geo;
      switch (gtype) {
      case POINT:
	geo.boundaryMark() = k; // geometry region(material marker)
	k = N_POINT_NODE;
	geo.vertex().resize(k);
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  geo.vertex(l) = ind[m];
	}
	node.push_back(geo);
	break;
      case LINE:
	geo.boundaryMark() = k; // geometry region(material marker)
	k = N_LINE_NODE;
	geo.vertex().resize(k);
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  geo.vertex(l) = ind[m];
	}
	line.push_back(geo);
	break;
      case TRIANGLE:
	geo.boundaryMark() = k; // geometry region(material marker)
	k = N_TRIANGLE_NODE;
	geo.vertex().resize(k);
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  geo.vertex(l) = ind[m];
	}
	surface.push_back(geo);
	break;
      case QUADRANGLE:
	geo.boundaryMark() = k; // geometry region(material marker)
	k = N_QUADRANGLE_NODE;
	geo.vertex().resize(k);
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  geo.vertex(l) = ind[m];
	}
	surface.push_back(geo);
	break;
      case TETRAHEDRON:
	k = N_TETRAHEDRON_NODE;
	ele[j].vertex.resize(k);			
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  ele[j].vertex[l] = ind[m];
	}
	/**
	 * 重要：新版的 Gmsh 的四面体的定向和 AFEPack 的四面体的定向是相
	 * 反的，我们通过将 0, 1 两个点的顺序掉换来改变四面体的定向！
	 */
	l = ele[j].vertex[1]; 
	ele[j].vertex[1] = ele[j].vertex[0];
	ele[j].vertex[0] = l;
	ele[j ++].template_element = 0;
	break;
      case HEXAHEDRON:
	k = N_HEXAHEDRON_NODE;
	ele[j].vertex.resize(k);			
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  ele[j].vertex[l] = ind[m];
	}
	ele[j ++].template_element = 1;
	break;
      case PRISM:
	k = N_PRISM_NODE;
	ele[j].vertex.resize(k);			
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  ele[j].vertex[l] = ind[m];
	}
	ele[j ++].template_element = 2;
	break;
      case PYRAMID:
	k = N_PYRAMID_NODE;
	ele[j].vertex.resize(k);			
	for (l = 0;l < k;l ++) { // index of node
	  is >> m;
	  ele[j].vertex[l] = ind[m];
	}
	ele[j ++].template_element = 3;
	break;
      default:
	assert(false);
      }
    }
    ele.resize(j);
    is >> dummy;/// $EndElements

  }
};

#else
void GmshMesh::readData(const std::string& filename)
{
  std::cerr << "Reading in gmsh data ..." << std::endl;
  int i, j, k, l, m, n, n_point, n_geometry, n_tags;
  GmshMesh::GeometryType gtype;
  std::string dummy;
  std::ifstream is(filename.c_str());
  is >> dummy;/// $MeshFormat
  is >> dummy;/// integer : for example, 2 means the version of Gmsh is 2.0
  is >> dummy;/// integer : for example, 0 means the type of data file is ASCII
  is >> dummy;/// integer : the size of the floating point numbers used in the file.
  is >> dummy;/// $EndMeshFormat$
  is >> dummy;/// $Nodes
  std::cerr << "\tReading nodes data ..." << std::endl;
  is >> n_point;
  pnt.resize(n_point);
  std::vector<int> ind0(n_point);
  for (i = 0, k = 0;i < n_point;i ++) {
    is >> j;
    if (j > k) k = j;
    ind0[i] = j;
    is >> pnt[i];
  }
  std::vector<int> ind(k+1, -1);
  for (i = 0;i < n_point;i ++) ind[ind0[i]] = i;
  is >> dummy;/// $EndNodes
  std::cerr << "\tReading geometry data ..." << std::endl;
  is >> dummy;/// $Elements
  is >> n_geometry;
  ele.resize(n_geometry);
  for (i = 0, j = 0;i < n_geometry;i ++) {
    is >> k; // geometry number
    is >> gtype; // geometry type
    is >> n_tags; // number of tags;
    is >> k; // reading the tags, the first tag is the boundary mark
    for (l = 1;l < n_tags;++ l) is >> m;
    GeometryBM geo;
    switch (gtype) {
    case POINT:
      geo.boundaryMark() = k; // geometry region(material marker)
      k = N_POINT_NODE;
      geo.vertex().resize(k);
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	geo.vertex(l) = ind[m];
      }
      node.push_back(geo);
      break;
    case LINE:
      geo.boundaryMark() = k; // geometry region(material marker)
      k = N_LINE_NODE;
      geo.vertex().resize(k);
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	geo.vertex(l) = ind[m];
      }
      line.push_back(geo);
      break;
    case TRIANGLE:
      geo.boundaryMark() = k; // geometry region(material marker)
      k = N_TRIANGLE_NODE;
      geo.vertex().resize(k);
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	geo.vertex(l) = ind[m];
      }
      surface.push_back(geo);
      break;
    case QUADRANGLE:
      geo.boundaryMark() = k; // geometry region(material marker)
      k = N_QUADRANGLE_NODE;
      geo.vertex().resize(k);
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	geo.vertex(l) = ind[m];
      }
      surface.push_back(geo);
      break;
    case TETRAHEDRON:
      k = N_TETRAHEDRON_NODE;
      ele[j].vertex.resize(k);			
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	ele[j].vertex[l] = ind[m];
      }
      /**
       * 重要：新版的 Gmsh 的四面体的定向和 AFEPack 的四面体的定向是相
       * 反的，我们通过将 0, 1 两个点的顺序掉换来改变四面体的定向！
       */
      l = ele[j].vertex[1]; 
      ele[j].vertex[1] = ele[j].vertex[0];
      ele[j].vertex[0] = l;
      ele[j ++].template_element = 0;
      break;
    case HEXAHEDRON:
      k = N_HEXAHEDRON_NODE;
      ele[j].vertex.resize(k);			
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	ele[j].vertex[l] = ind[m];
      }
      ele[j ++].template_element = 1;
      break;
    case PRISM:
      k = N_PRISM_NODE;
      ele[j].vertex.resize(k);			
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	ele[j].vertex[l] = ind[m];
      }
      ele[j ++].template_element = 2;
      break;
    case PYRAMID:
      k = N_PYRAMID_NODE;
      ele[j].vertex.resize(k);			
      for (l = 0;l < k;l ++) { // index of node
	is >> m;
	ele[j].vertex[l] = ind[m];
      }
      ele[j ++].template_element = 3;
      break;
    default:
      assert(false);
    }
  }
  ele.resize(j);
  is >> dummy;/// $EndElements
};
#endif


bool are_two_pntVecs_equal(const std::vector<AFEPack::Point<3> >& pv1,
			   const std::vector<AFEPack::Point<3> >& pv2)
{
  const int& n_pv1 = pv1.size();
  const int& n_pv2 = pv2.size();

  //bool is_equal = false;

  assert(n_pv1 == n_pv2);

  for(int i = 0;i < n_pv1;++ i){
    const AFEPack::Point<3>& p1 = pv1[i];
    bool is_found = false;
    for(int j = 0;j < n_pv1;++ j){
      const AFEPack::Point<3>& p2 = pv2[j];

      if(distance(p1, p2) < 1.0e-06){/// p1 and p2 are two same points
	is_found = true;
	break;
      }
    }

    if(!is_found){/// if p1 is not found from pv2
      return false;
    }
  }

  return true;
}


void GmshMesh::addBoundaryMark(Mesh<3,3>& m)
{
  /// now we have already got the mesh m, we need to add the boundary
  /// mark according to the information in boundary_ele

  /// The logic is we loop all element in boundary_ele, and find its
  /// corresponding element in m, and use boundary_ele.boundaryMark to
  /// assign the boundary mark from m.

  /// Plz be noted that, when using generateMesh(), the indices of
  /// those nodes have been revised. So we dont have the
  /// correspondence anymore. To check if two geometries are equal to
  /// each other, we need to verify that the COORDINATES of those
  /// nodes are equals.

  std::cerr << "Assigning boundary mark in mesh data ..." << std::endl;
  
  const int & n_bnd_ele = boundary_ele.size();

  for(int i = 0;i < n_bnd_ele;++ i){
    const int & dim = boundary_ele[i].dim;
    const int & bm = boundary_ele[i].boundaryMark;
    const std::vector<int>& bnd_ele_vtx = boundary_ele[i].vertex;
    const int & n_bnd_ele_vtx = boundary_ele[i].vertex.size();
    std::vector<AFEPack::Point<3> > boundary_ele_pnts(n_bnd_ele_vtx);
    for(int j = 0;j < n_bnd_ele_vtx;++ j){
      boundary_ele_pnts[j] = pnt[bnd_ele_vtx[j]];
    }

    const int & n_dim_element = m.n_geometry(dim);

    for(int j = 0;j < n_dim_element;++ j){
      GeometryBM& the_geo = m.geometry(dim, j);/// or j, dim ?
      const std::vector<int>& geo_vtx = the_geo.vertex();
      const int& n_geo_vtx = geo_vtx.size();

      assert(n_geo_vtx == n_bnd_ele_vtx);/// otherwise, there is
					 /// something wrong.
      
      std::vector<AFEPack::Point<3> > the_geo_pnts(n_geo_vtx);
      for(int l = 0;l < n_bnd_ele_vtx;++ l){
	the_geo_pnts[l] = m.point(geo_vtx[l]);
      }
      
      if(are_two_pntVecs_equal(the_geo_pnts, boundary_ele_pnts)){
	the_geo.boundaryMark() = bm;
	break;	
      }

      if(j == n_dim_element -1){
	std::cout << "There is something wrong on finding matching geo with "
		  << dim << "-dimensional Physical Entity tag :"
		  << bm << std::endl;
      }
    }
  }

  /// at this stage, all entities' tags listed in x.msh file have been
  /// added on x.mesh. However, there is chance that some low
  /// dimension geometries will be missed (those geometries generated
  /// by meshing). So in the following, we loop all (DIM - 1)
  /// geometries, and use the tag from those geometries to assign all
  /// its lower geometries.

  const int& n_surface = m.n_geometry(3 - 1);
  for(int i = 0;i < n_surface;++ i){
    GeometryBM& the_geo = m.geometry(3 - 1, i);
    const int& bm = the_geo.boundaryMark();

    if(bm == 0) continue;

    std::vector<int>& the_bnd = the_geo.boundary();
    const int& n_bnd = the_geo.n_boundary();
    for(int j = 0;j < n_bnd;++ j){
      GeometryBM& the_bnd_geo = m.geometry(3 - 2, the_bnd[j]); /// for curve's boundary mark
      the_bnd_geo.boundaryMark() = bm;

      std::vector<int>& the_bnd_bnd = the_bnd_geo.boundary();
      const int& n_bnd_bnd = the_bnd_geo.n_boundary();
      for(int k = 0;k < n_bnd_bnd;++ k){
	GeometryBM& the_bnd_bnd_geo = m.geometry(3 - 3, the_bnd_bnd[k]);/// for point's boundary mark
	the_bnd_bnd_geo.boundaryMark() = bm;
      }
    }
  }
}


void GmshMesh::generateMesh(Mesh<3,3>& m)
{
  SimplestMesh<3,3>::generateMesh(m);

  int i, j;

  return;
  /**
   * 我们在generateMesh中对顶点的指标进行了整理，下面的代码现在暂时没
   * 有法子用了。
   */

  /// 重新整理顶点坐标
  int n_vtx = m.n_geometry(0);
  std::vector<int> index(m.n_geometry(0));
  for (i = 0;i < n_vtx;i ++) {
    index[m.geometry(0,i).vertex(0)] = i;
  }

  /// 整理面和线段的顶点指标	
  std::list<GeometryBM>::iterator the_geo, end_geo;
  the_geo = surface.begin();
  end_geo = surface.end();
  for (;the_geo != end_geo;the_geo ++) {
    for (i = 0;i < the_geo->n_vertex();i ++) {
      j = the_geo->vertex(i);
      the_geo->vertex(i) = index[j];
    }
  }
  the_geo = line.begin();
  end_geo = line.end();
  for (;the_geo != end_geo;the_geo ++) {
    for (i = 0;i < the_geo->n_vertex();i ++) {
      j = the_geo->vertex(i);
      the_geo->vertex(i) = index[j];
    }
  }
  the_geo = node.begin();
  end_geo = node.end();
  for (;the_geo != end_geo;the_geo ++) {
    for (i = 0;i < the_geo->n_vertex();i ++) {
      j = the_geo->vertex(i);
      the_geo->vertex(i) = index[j];
    }
  }	

  the_geo = surface.begin();
  end_geo = surface.end();
  for (;the_geo != end_geo;the_geo ++) {
    for (i = 0;i < m.n_geometry(2);i ++) {
      if (isSame(m.geometry(2,i), *the_geo)) {
	m.boundaryMark(2,i) = the_geo->boundaryMark();
	break;
      }
    }
    assert(i < m.n_geometry(2));
  }
  for (i = 0;i < m.n_geometry(2);i ++) {
    if (m.boundaryMark(2,i) == 0) continue;
    for (j = 0;j < m.geometry(2,i).n_boundary();j ++) {
      m.boundaryMark(1,m.geometry(2,i).boundary(j)) = m.boundaryMark(2,i);
    }
  }
  the_geo = line.begin();
  end_geo = line.end();
  for (;the_geo != end_geo;the_geo ++) {
    for (i = 0;i < m.n_geometry(1);i ++) {
      if (isSame(m.geometry(1,i), *the_geo)) {
	m.boundaryMark(1,i) = the_geo->boundaryMark();
	break;
      }
    }
    assert(i < m.n_geometry(1));
  }
  for (i = 0;i < m.n_geometry(1);i ++) {
    if (m.boundaryMark(1,i) == 0) continue;
    for (j = 0;j < m.geometry(1,i).n_boundary();j ++) {
      m.boundaryMark(0,m.geometry(1,i).boundary(j)) = m.boundaryMark(1,i);
    }
  }
  the_geo = node.begin();
  end_geo = node.end();
  for (;the_geo != end_geo;the_geo ++) {
    for (i = 0;i < m.n_geometry(0);i ++) {
      if (isSame(m.geometry(0,i), *the_geo)) {
	m.boundaryMark(0,i) = the_geo->boundaryMark();
	break;
      }
    }
    assert(i < m.n_geometry(0));
  }
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

