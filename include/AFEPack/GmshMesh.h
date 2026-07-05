///////////////////////////////////////////////////////////////////////////
// GmshMesh.h : by R.Lie
//

#ifndef _GmshMesh_h_
#define _GmshMesh_h_

#include <iostream>
#include <string>
#include <fstream>

#include <AFEPack/Geometry.h>

AFEPACK_OPEN_NAMESPACE

/**
 * This class provides facilities to asscess the mesh data file generated
 * by the mesh generator \p{gmsh}. For 3-dimensional only. Though we can
 * read in very flexible data format, we currently only use it to read in
 * pure tetrahedron mesh.
 */
#if 1/// add a structure SimplestGeometryBM, i.e., add boundary mark info.
class GmshMesh : public SimplestMesh<3,3>
{
public:
  typedef int GeometryType;
  static const GeometryType POINT		= 15;
  static const GeometryType LINE		= 1;
  static const GeometryType TRIANGLE		= 2;
  static const GeometryType QUADRANGLE		= 3;
  static const GeometryType TETRAHEDRON		= 4;
  static const GeometryType HEXAHEDRON		= 5;
  static const GeometryType PRISM		= 6;
  static const GeometryType PYRAMID		= 7;
protected:
  struct SimplestGeometryBM : SimplestGeometry {
  public:
    //int template_element; /**< Index of template element for the simplest geometry. */
    //std::vector<int> vertex; /**< The vertex array of this element. */
    int dim;
    int boundaryMark;
  };
  std::vector<SimplestGeometryBM> boundary_ele;
private:
  std::list<GeometryBM> node;
  std::list<GeometryBM> line;
  std::list<GeometryBM> surface;
public:
  GmshMesh();
  virtual ~GmshMesh();
public:
  void readData(const std::string&);
  virtual void generateMesh(Mesh<3,3>& m);
  virtual void addBoundaryMark(Mesh<3, 3>& m);
  // DeclException1(ExcMeshData, 
  // 		 char *, 
  // 		 << "Mesh data error: " << arg1);
};
#else
class GmshMesh : public SimplestMesh<3,3>
{
public:
  typedef int GeometryType;
  static const GeometryType POINT		= 15;
  static const GeometryType LINE		= 1;
  static const GeometryType TRIANGLE		= 2;
  static const GeometryType QUADRANGLE		= 3;
  static const GeometryType TETRAHEDRON		= 4;
  static const GeometryType HEXAHEDRON		= 5;
  static const GeometryType PRISM		= 6;
  static const GeometryType PYRAMID		= 7;
private:
  std::list<GeometryBM> node;
  std::list<GeometryBM> line;
  std::list<GeometryBM> surface;
public:
  GmshMesh();
  virtual ~GmshMesh();
public:
  void readData(const std::string&);
  virtual void generateMesh(Mesh<3,3>& m);
  DeclException1(ExcMeshData, 
		 char *, 
		 << "Mesh data error: " << arg1);
};
#endif
AFEPACK_CLOSE_NAMESPACE

#endif // _GmshMesh_h_

//
// end of file
///////////////////////////////////////////////////////////////////////////
