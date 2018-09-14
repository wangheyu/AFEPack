/**
 * 将 Tecplot 格式的数据文件转换为 OpenDX 格式的数据文件
 * 
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char * argv[])
{
  if (argc != 5) {
    std::cout << "Usage: " << argv[0]
	      << " input_file output_file n_coord_dim n_var_dim\n"
	      << "\tinput_file: Input file with Open DX data file format\n"
	      << "\toutput_file: Output file with Tecplot data file format\n"
	      << "\tn_coord_dim: Dimension of coordinates\n"
	      << "\tn_var_dim: Dimension of variables\n"
	      << std::endl;
    return 1;
  }
  std::string infile(argv[1]);
  std::string outfile(argv[2]);
  int n_coord_dim = atoi(argv[3]);
  int n_var_dim = atoi(argv[4]);

  int n_point, n_element;
  char buffer[256];
  std::ifstream is(infile.c_str());
  is.getline(buffer, 256, '\n');
  std::cout << buffer << std::endl;
  is.getline(buffer, 256, '=');
  std::cout << buffer << "=" << std::flush;
  is >> n_point; // read the first line of the information
  is.getline(buffer, 256, '=');
  std::cout << n_point << buffer << "=" << std::flush;
  is >> n_element;
  is.getline(buffer, 256, '\n');
  std::cout << n_element << buffer << std::endl;
  std::vector<std::vector<double> > point(n_point,
					  std::vector<double>(n_coord_dim));
  std::vector<std::vector<double> > var(n_point,
					std::vector<double>(n_var_dim));
  std::cout << "Reading in points and variables ..." << std::flush;
  for (int i = 0;i < n_point;i ++) {
    for (int j = 0;j < n_coord_dim;j ++)
      is >> point[i][j];
    for (int j = 0;j < n_var_dim;j ++)
      is >> var[i][j];
  }

  std::cout << " OK!\nReading in elements ... " << std::flush;
  std::vector<std::vector<int> > element(n_element,
					 std::vector<int>(n_coord_dim + 1));
  for (int i = 0;i < n_element;i ++)
    for (int j = 0;j <= n_coord_dim;j ++)
      is >> element[i][j];
  
  std::cout << " OK!\nWriting points ... " << std::flush;
  std::ofstream os(outfile.c_str());
  os << "object 1 class array type float rank 1 shape " 
     << n_coord_dim << " item " 
     << n_point << " data follows\n";
  for (int i = 0;i < n_point;i ++) {
    for (int j = 0;j < n_coord_dim;j ++)
      os << point[i][j] << "\t";
    os << "\n";
  }

  std::cout << " OK!\nWriting connections ... " << std::flush;
  os << "\nobject 2 class array type int rank 1 shape "
     << n_coord_dim + 1 << " item "   
     << n_element << " data follows\n";
  for (int i = 0;i < n_element;i ++) {
    for (int j = 0;j <= n_coord_dim;j ++)
      os << element[i][j]-1 << "\t";
    os << "\n";
  }

  std::cout << " OK!\nWriting variables ... " << std::flush;
  os << "attribute \"element type\" string \"";
  switch(n_coord_dim) {
  case 2:
    os << "triangles\"\n";
    break;
  case 3:
    os << "tetrahedra\"\n";
    break;
  default:
    std::cout << "Can't handle such case.\n";
    exit(1);
  }
  os << "attribute \"ref\" string \"positions\"\n\n";
	
  os << "object 3 class array type float rank 1 shape "
     << n_var_dim << " item "
     << n_point << " data follows\n";
  for (int i = 0;i < n_point;i ++) {
    for (int j = 0;j < n_var_dim;j ++)
      os << var[i][j] << "\t";
    os << "\n";
  }
  os << "attribute \"dep\" string \"positions\"\n\n";
	
  os << "object \"FEMFunction\" class field\n"
     << "component \"positions\" value 1\n"
     << "component \"connections\" value 2\n"
     << "component \"data\" value 3\n"
     << "end\n";
  os.close();
  std::cout << " OK!" << std::endl;

  return 0;
};

/**
 * end of file
 * 
 */
