/**
 * @file   ex03.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Tue Oct  6 16:16:47 2009
 * 
 * @brief  测试用 HGeometryTree 的序列化函数，来完成一棵树的拆分。
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

  int n_root_ele = h_tree.n_rootElement();
  std::ofstream os("htree1.dat");
  boost::archive::text_oarchive oa1(os);
  boost::serialization::save(oa1, h_tree, 0, n_root_ele/2);
  os.close();

  os.open("htree2.dat");
  boost::archive::text_oarchive oa2(os);
  boost::serialization::save(oa2, h_tree, n_root_ele/2, n_root_ele);
  os.close();

  HGeometryTree<DIM,DOW> sub_tree_1, sub_tree_2;
  std::ifstream is("htree1.dat");
  boost::archive::text_iarchive ia1(is);
  boost::serialization::load(ia1, sub_tree_1);
  is.close();

  is.open("htree2.dat");
  boost::archive::text_iarchive ia2(is);
  boost::serialization::load(ia2, sub_tree_2);
  is.close();

  return 0;
}

/**
 * end of file
 * 
 */
