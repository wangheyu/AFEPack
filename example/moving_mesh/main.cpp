////////////////////////////////////////////////////////////////////////////////////////////
// main.cpp :
//

#include "burgers.h"

int main(int argc, char * argv[])
{
	try {
		if (argc != 2) {
			std::cout << "Usage: " << argv[0] 
				<< " meshfile" << std::endl;
			return 1;
		}
		burgers the_app(argv[1]);
		the_app.run();
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

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
