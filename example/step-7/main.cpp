

// main.cpp

#include "EllipticSystem.h"

int main(int argc, char * argv[])
{
	if (argc != 2) {
		std::cout << "Usage: " 
			<< argv[0]
			<< " mesh_file_name"
			<< std::endl;
		return 1;
	}
	try {
		EllipticSystem elliptic_system(argv[1]);
		elliptic_system.run();
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
};

// end of file

