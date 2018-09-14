/**
 * @file   main.cpp
 * @author Ruo Li <rli@aztec>
 * @date   Fri Oct 23 23:22:55 2009
 * 
 * @brief  
 * 
 * 
 */

#include "SCL.h"

extern void init_data(SCL& scl, int argc, char * argv[]);
SCL app;

void save_data(const char * dummy) {
  app.save_data();
}

int main(int argc, char * argv[])
{
  try {
    MPI_Init(&argc, &argv);

    AFEPack::MPI::registerController("save", &save_data, "save solution.");
    AFEPack::MPI::controllerScript("run.sh", MPI_COMM_WORLD);

    init_data(app, argc, argv);
    app.run();

    MPI_Finalize();
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
};

/**
 * end of file
 * 
 */

