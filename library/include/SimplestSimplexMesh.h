/**
 * @file   SimplestSimplexMesh.h
 * @author Ruo Li <rli@aztec>
 * @date   Sat Mar 22 11:24:27 2008
 * 
 * @brief  
 * 
 * 
 */

#ifndef __SimplestSimplexMesh_h__
#define __SimplestSimplexMesh_h__

#include "Geometry.h"

template <int DIM, int DOW=DIM>
  class SimplestSimplexMesh : public SimplestMesh<DIM,DOW>
  {
  public:
  SimplestSimplexMesh() {}
  virtual ~SimplestSimplexMesh() {}

  public:
  void readData(const std::string& file) {}
  };

template <int DOW>
class SimplestSimplexMesh<2,DOW> : public SimplestMesh<2,DOW>
{
 public:
  SimplestSimplexMesh<2,DOW>() 
    {
      this->setTemplateGeometry(*(new std::vector<TemplateGeometry<2> >(1)));
      this->templateGeometry(0).readData("triangle.tmp_geo");
    }
  virtual ~SimplestSimplexMesh()
  {
    delete &(this->templateGeometry());
  }

 public:
  void readData(const std::string& file)
  {
    std::cerr << "Reading in simplest triangle mesh data ..." << std::endl;

    u_int n_point, n_element;
    std::ifstream is(file.c_str());
    is >> n_point;
    this->point().resize(n_point);
    for (u_int i = 0;i < n_point;++ i) {
      is >> this->point(i);
    }

    is >> n_element;
    this->element().resize(n_element);
    for (u_int i = 0;i < n_element;++ i) {
      this->element(i).vertex.resize(4);
      for (u_int k = 0;k < 3;++ k) {
        is >> this->element(i).vertex[k];
      }
      this->element(i).template_element = 0;
    }

    is.close();
  }
};

template <int DOW>
class SimplestSimplexMesh<3,DOW> : public SimplestMesh<3,DOW>
{
 public:
  SimplestSimplexMesh<3,DOW>() 
    {
      this->setTemplateGeometry(*(new std::vector<TemplateGeometry<3> >(1)));
      this->templateGeometry(0).readData("tetrahedron.tmp_geo");
    }
  virtual ~SimplestSimplexMesh()
  {
    delete &(this->templateGeometry());
  }

 public:
  void readData(const std::string& file)
  {
    std::cerr << "Reading in simplest tetrahedron mesh data ..." << std::endl;

    u_int n_point, n_element;
    std::ifstream is(file.c_str());
    is >> n_point;
    this->point().resize(n_point);
    for (u_int i = 0;i < n_point;++ i) {
      is >> this->point(i);
    }

    is >> n_element;
    this->element().resize(n_element);
    for (u_int i = 0;i < n_element;++ i) {
      this->element(i).vertex.resize(4);
      for (u_int k = 0;k < 4;++ k) {
        is >> this->element(i).vertex[k];
      }
      this->element(i).template_element = 0;
    }

    is.close();
  }
};

#endif // __SimplestSimplexMesh_h__

/**
 * end of file
 * 
 */
