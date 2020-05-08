/*
in order to use armadillo hdf5 save/load functionality
hdf5.h must be available ==> additional  include dir hdf5/include/
link with the hdf5 library:
		static linking ==> libszip.lib;libzlib.lib;libhdf5.lib;libhdf5_cpp.lib;
		dynamic linking ==> szip.lib;zlib.lib;hdf5.lib;hdf5_cpp.lib;
*/
#include "pch.h"
#include "BHModel.h"

using namespace std::chrono;
namespace fs = std::filesystem;

int numSites = 1;
int numParticles = 1;
double intStrength; //the interaction strength U (with J set to 1)

void argParser(int argc, char* argv[]);
void getModel(int ns, int np);

int main(int argc, char* argv[]) {
	argParser(argc, argv);
	std::cout << "\nnumber of Sites = " << numSites << "\n";
	std::cout << "number of Particles = " << numParticles << "\n\n";

#if 1
	getModel(numSites, numParticles);

	//TODO: write the Hamiltonian matrix

	//arma::Mat<int> data;
	//data.load(arma::hdf5_name(h5_file_name, "Dataset1"));
	////because h5 is column based, we need to take the transpose
	//data = data.t();//now data is a nBasis*basisLen matrix
	//std::cout << "\nnumber of basis vectors = " << data.n_rows << "\n";
	//std::cout << "each vector has dim = " << data.n_cols << "\n\n";
	//data.print("basis matrix is:");
#else
	arma::mat A = arma::randu(5, 5);
	arma::mat B = arma::randu(5, 5);
	arma::umat C = (A >= B);
	C.print();
#endif
}

void inline getModel(int ns, int np) {
	BHModel bh(ns, np);
	auto start = high_resolution_clock::now();
	bh.mkBasisMatrix();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() / 1000.0 << " seconds has elapsed" << std::endl;
	//bh.generateH5();
	//bh.readh5dims();
}