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
double intStrength; //the interaction strength U (with J set to 1)


void argParser(int argc, char* argv[]);
void getModel(int ns, int np);

int main(int argc, char* argv[]) {
	argParser(argc, argv);
	std::cout << "number of Sites = " << numSites << "\n";
	std::cout << "number of Particles = " << numParticles << "\n";

	getModel(numSites, numParticles);
	
	arma::Mat<int> data;
	data.load(arma::hdf5_name(h5_file_name, "Dataset1"));
	//because h5 is column based, we need to take the transpose
	data = data.t();//now data is a nBasis*basisLen matrix
	data.print("data matrix is:");
}

void inline getModel(int ns, int np) {
	BHModel bh(ns, np);
	bh.generateHDF5DataFile();
	bh.readh5dims();

	int nBasis = bh.nBasis;
	int basisLen = bh.basisLen;
	std::cout << "\nnumber of basis vectors = " << nBasis << "\n";
	std::cout << "each vector has dim = " << basisLen << "\n\n";
}