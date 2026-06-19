/**
 * @file   main.cpp
 * @author Robert Lie
 * @date   Wed Feb 21 16:44:21 2007
 *
 * @brief  主文件，我们在这里开始运行，其中的这几行程序的意义是显而易
 *         见的，就不多解释了。
 *
 */

#include "EllipticEquation.h"

int main(int argc, char * argv[])
{
  try {
    EllipticEquation the_app;
    std::string filename = argv[1];
    the_app.initialize(filename);
    the_app.buildSpace();
    the_app.solve();
  }
  catch(std::exception& e) {
    std::cerr << "Exception caughted:" << std::endl
              << e.what ()
              << std::endl;
  }
  catch(...) {
    std::cerr << "Exception caughted:" << std::endl
              << "unknown exception caughted."
              << std::endl;
  }

  return 0;
}

/**
 * end of file
 *
 */
