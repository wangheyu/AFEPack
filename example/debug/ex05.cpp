/**
 * @file   ex05.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct  6 13:16:13 2009
 * 
 * @brief 测试有限元函数序列化函数 export_fe_func 和 import_fe_func。
 * 
 * 
 */

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <AFEPack/Serialization.h>
#include <AFEPack/Operator.h>
#include <AFEPack/mpi/MPI_Migration.h>

#define DOW 2
#define DIM DOW

double _u_(const double * p)
{
  return sin(12*p[0])*sin(12*p[1]);
}

int main()
{
  HGeometryTree<DIM,DOW> h_tree;
  h_tree.readEasyMesh("D");
  IrregularMesh<DIM,DOW> ir_mesh(h_tree);
  ir_mesh.globalRefine(2);
  ir_mesh.semiregularize();
  ir_mesh.regularize();

  TemplateGeometry<DIM> template_geometry;
  template_geometry.readData("triangle.tmp_geo");
  CoordTransform<DIM,DIM> coord_transform;
  coord_transform.readData("triangle.crd_trs");
  TemplateDOF<DIM> template_dof(template_geometry);
  template_dof.readData("triangle.2.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> basis_function(template_dof);
  basis_function.readData("triangle.2.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > template_element(1);
  template_element[0].reinit(template_geometry,
			     template_dof,
			     coord_transform,
			     basis_function);

  RegularMesh<DIM,DOW>& mesh = ir_mesh.regularMesh();
  FEMSpace<double,DIM,DOW> fem_space(mesh, template_element);
	
  int n_element = mesh.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space, i, 0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  FEMFunction<double,DIM,DOW> u_h(fem_space);
  Operator::L2Interpolate(_u_, u_h);

  Migration::initialize();
  Migration::clear_data_buffer(h_tree);
  const Migration::data_id_t& id = Migration::register_data_name("u_h", true);
  Migration::export_fe_func(u_h, id);

  u_int n_root_ele = h_tree.n_rootElement();
  u_int n_output_ele = 0.3*n_root_ele;
  std::ofstream os("htree1.dat");
  boost::archive::text_oarchive oa1(os);
  IrregularMesh<DIM,DOW>::RootIterator
    the_root_ele = ir_mesh.beginRootElement(),
    end_root_ele = ir_mesh.endRootElement();
  for (u_int i = 0;i < n_output_ele;++ i) {
    oa1 & *the_root_ele;
    ++ the_root_ele;
  }
  os.close();

  ir_mesh.clear();
  h_tree.clear();
  ir_mesh.reinit(h_tree, true);

  std::ifstream is("htree1.dat");
  boost::archive::text_iarchive ia1(is);
  for (u_int i = 0;i < n_output_ele;++ i) {
    HElement<DIM,DOW> * ele = new HElement<DIM,DOW>();
    ia1 & *ele;
    ir_mesh.rootElement().push_back(ele);
    h_tree.rootElement().push_back(ele->h_element);
  }
  is.close();

  ir_mesh.semiregularize();
  ir_mesh.regularize();

  RegularMesh<DIM,DOW>& mesh1 = ir_mesh.regularMesh();
  fem_space.reinit(mesh1, template_element);
	
  n_element = mesh1.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space, i, 0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  FEMFunction<double,DIM,DOW> u1_h(fem_space);
  Migration::import_fe_func(u1_h, id);  

  u1_h.writeOpenDXData("u1_h.dx");

  return 0;
}

/**
 * end of file
 * 
 */
