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

#include "GmshMesh.h"

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
      Assert(false, ExcMeshData("Unknown geometry type."));
    }
  }
  ele.resize(j);
  is >> dummy;/// $EndElements
};

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
    Assert(i < m.n_geometry(2), ExcInternalError());
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
    Assert(i < m.n_geometry(1), ExcInternalError());
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
    Assert(i < m.n_geometry(0), ExcInternalError());
  }
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */

