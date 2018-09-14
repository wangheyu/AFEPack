///////////////////////////////////////////////////////////////////////////
// DBMesh.cpp : by R.Lie
//

#include "DBMesh.h"

AFEPACK_OPEN_NAMESPACE

DBMesh::DBMesh()
{
  setTemplateGeometry(*(new std::vector<TemplateGeometry<2> >(2)));
  templateGeometry(0).readData("triangle.tmp_geo");
  templateGeometry(1).readData("rectangle.tmp_geo");
};

DBMesh::~DBMesh()
{
  delete &(templateGeometry());
};

void DBMesh::readData(const std::string& filename)
{
  int i;
  std::string dummy;
	
  std::ifstream is(filename.c_str());
  element().clear();
  do {
    getline(is, dummy);
    if (dummy == "Vertices")
      readNode(is);
    else if (dummy == "Triangles")
      readTriangleElement(is);
    else if (dummy == "Quadrilaterals")
      readQuadrilateralElement(is);
    else if (dummy == "Dimension") {
      is >> i;
      Assert(i == 2, ExcMeshData("the dimension of the mesh should be 2"));
    }
    else if (dummy == "End")
      break;
  } while(!is.eof());
  is.close();
};

void DBMesh::generateMesh(Mesh<2,2>& m)
{
  SimplestMesh<2,2>::generateMesh(m);
	
  int i, j, k;
  for (i = 0;i < m.n_geometry(0);i ++) {
    m.boundaryMark(0, i) = bm[m.geometry(0, i).vertex(0)];
  }
  for (i = 0;i < m.n_geometry(1);i ++) {
    j = m.boundaryMark(0, m.geometry(1,i).vertex(0));
    k = m.boundaryMark(0, m.geometry(1,i).vertex(1));
    if (j == 0 && k == 0)
      m.boundaryMark(1, i) = 0;
    else
      m.boundaryMark(1, i) = std::min(j,k);
  }
  for (i = 0;i < m.n_geometry(2);i ++)
    m.boundaryMark(2, i) = 0;
};

void DBMesh::readNode(std::istream& is)
{
  int i, n_node;

  std::cout << "Reading node data ..." << std::endl;
  is >> n_node;
  point().resize(n_node);
  bm.resize(n_node);
  for (i = 0;i < n_node;i ++) {
    is >> point(i)[0]
       >> point(i)[1];
    is >> bm[i];
  }
};

void DBMesh::readTriangleElement(std::istream& is)
{
  int i, j, l, v[3], n_element;
	
  std::cout << "Reading triangle element data ..." << std::endl;
  is >> n_element;
  l = element().size();
  element().resize(l + n_element);
  for (i = 0;i < n_element;i ++) {
    element(l+i).template_element = 0;
    elementVertex(l+i).resize(3);
    for (j = 0;j < 3;j ++) {
      is >> v[j];
      v[j] --;
      elementVertex(l+i,j) = v[j];
    }
    is >> j;
  }
};

void DBMesh::readQuadrilateralElement(std::istream& is)
{
  int i, j, l, v[4], n_element;
	
  std::cout << "Reading quadrilateral element data ..." << std::endl;
  is >> n_element;
  l = element().size();
  element().resize(l + n_element);
  for (i = 0;i < n_element;i ++) {
    element(l+i).template_element = 1;
    elementVertex(l+i).resize(4);
    for (j = 0;j < 4;j ++) {
      is >> v[j];
      v[j] --;
      elementVertex(l+i,j) = v[j];
    }
    is >> j;
  }
};

AFEPACK_CLOSE_NAMESPACE

//
// end of file
///////////////////////////////////////////////////////////////////////////
