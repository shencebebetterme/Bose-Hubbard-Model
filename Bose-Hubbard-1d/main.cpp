/*
in order to use armadillo hdf5 save/load functionality
hdf5.h must be available ==> additional  include dir hdf5/include/
link with the hdf5 library:
		static linking ==> libszip.lib;libzlib.lib;libhdf5.lib;libhdf5_cpp.lib;
		dynamic linking ==> szip.lib;zlib.lib;hdf5.lib;hdf5_cpp.lib; 
*/
#include "pch.h"
#include "BHModel.h"

int numSites = 1;
int numParticles = 1;
double intStrength; //the interaction strength U (with J set to 1).{

void argParser(int argc, char* argv[]);

int main(int argc, char* argv[]) {
	
	argParser(argc, argv);
	std::cout << "number of sites = " << numSites << "\n";
	std::cout << "number of Particles = " << numParticles << "\n";
	/*
	std::cout << "1d Bose-Hubbard model.\n\nPlease input the number of sites M:" << std::endl;
	std::cin >> numSites;
	std::cout << "\nPlease enter the number of total particles N:" << std::endl;
	std::cin >> numParticles;
	*/


	BHModel bh(numSites, numParticles);
	bh.generateHDF5DataFile();
	bh.readh5dims();

	int nBasis = bh.nBasis;
	int basisLen = bh.basisLen;
	//arma::Mat<int> data = bh.readData();
	//data.print("data:");
	arma::Mat<int> data;
	data.load(arma::hdf5_name("data.h5", "Dataset1"));
	//because h5 is column based, we need to take the transpose
	data = data.t();//now data is a nBasis*basisLen matrix
	data.print("data matrix is:");

}


void argParser(int argc, char* argv[]) {
	//no arguments, use the default value
	if (argc == 1) {
		return;
	}

	//help
	if (std::string(argv[1]) == "-h") {
		std::cout << "pass command line arguments\n"
			<< "-nS\tnumber of Sites, default value = 1\n"
			<< "-nP\tnumber of Particles, default value = 1"
			<< std::endl;
		std::exit(EXIT_FAILURE);
	}

	enum argStr
	{
		num_Sites,
		num_Particles,
		inter_Strength
	};

	std::map<std::string, argStr> argMap;
	//initialize the map
	argMap["-nS"] = num_Sites;
	argMap["-nP"] = num_Particles;

	// the first arg is the funcion name
	// start from the second arg
	for (int i = 1; i < argc; i++) {
		std::string arg_i(argv[i]);
		// other command arguments have structure arg=val
		int equalSign = arg_i.find("=");
		// obtain the arg-val pair
		std::string arg_ic = arg_i.substr(0, equalSign);//arg_i content
		std::string arg_iv = arg_i.substr(equalSign + 1, -1);//arg_i value
		switch (argMap[arg_ic])
		{
		case num_Sites: numSites = std::stoi(arg_iv); break;
		case num_Particles: numParticles = std::stoi(arg_iv); break;
		}
	}
}
