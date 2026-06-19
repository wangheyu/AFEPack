/**
 * @file   ex04.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct  6 13:16:33 2009
 * 
 * @brief  测试 IrregularMesh 的序列化函数。
 * 
 * 
 */

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <AFEPack/Serialization.h>

#define DOW 2
#define DIM DOW

int main()
{
  HGeometryTree<DIM,DOW> h_tree;
  h_tree.readEasyMesh("D");
  IrregularMesh<DIM,DOW> ir_mesh(h_tree);
  ir_mesh.globalRefine(1);
  ir_mesh.semiregularize();
  ir_mesh.regularize();

  u_int n_root_ele = h_tree.n_rootElement();
  std::ofstream os("htree1.dat");
  boost::archive::text_oarchive oa1(os);
  boost::serialization::save(oa1, h_tree, 0, n_root_ele/2);
  boost::serialization::save(oa1, ir_mesh, 0, n_root_ele/2);
  os.close();

  os.open("htree2.dat");
  boost::archive::text_oarchive oa2(os);
  boost::serialization::save(oa2, h_tree, n_root_ele/2, n_root_ele);
  boost::serialization::save(oa2, ir_mesh, n_root_ele/2, n_root_ele);
  os.close();

  ir_mesh.clear();
  h_tree.clear();

  HGeometryTree<DIM,DOW> sub_tree_1, sub_tree_2;
  IrregularMesh<DIM,DOW> sub_mesh_1, sub_mesh_2;
  std::ifstream is("htree1.dat");
  boost::archive::text_iarchive ia1(is);
  boost::serialization::load(ia1, sub_tree_1);
  boost::serialization::load(ia1, sub_mesh_1);
  is.close();

  is.open("htree2.dat");
  boost::archive::text_iarchive ia2(is);
  boost::serialization::load(ia2, sub_tree_2);
  boost::serialization::load(ia2, sub_mesh_2);
  is.close();

  sub_mesh_1.reinit(sub_tree_1, true);
  sub_mesh_1.globalRefine(1);
  sub_mesh_1.semiregularize();
  sub_mesh_1.regularize();
  sub_mesh_1.regularMesh().writeOpenDXData("mesh1.dx");

  sub_mesh_2.reinit(sub_tree_2, true);
  sub_mesh_2.globalRefine(2);
  sub_mesh_2.semiregularize();
  sub_mesh_2.regularize();
  sub_mesh_2.regularMesh().writeOpenDXData("mesh2.dx");

  return 0;
}

/**
 * end of file
 * 
 */
