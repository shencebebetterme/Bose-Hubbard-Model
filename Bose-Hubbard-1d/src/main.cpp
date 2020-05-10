/*
in order to use armadillo hdf5 save/load functionality
hdf5.h must be available ==> additional  include dir hdf5/include/
link with the hdf5 library:
		static linking ==> libszip.lib;libzlib.lib;libhdf5.lib;libhdf5_cpp.lib;
		dynamic linking ==> szip.lib;zlib.lib;hdf5.lib;hdf5_cpp.lib;
*/
#include "pch.h"
#include "BHModel.h"
#include "Hamiltonian.h"

using namespace std::chrono;
namespace fs = std::filesystem;

//communicate with arg parser
int numSites = 1;
int numParticles = 1;
double intStrength = 1.0; //the interaction strength U/2 (with J set to 1)
int numEig = 1;
double tol = 0.0001;


/////////////////////////////////////////////////////////////////////////////////////////



void argParser(int argc, char* argv[]);
void getModel(int ns, int np);


int main(int argc, char* argv[]) {
	argParser(argc, argv);
	std::cout << "\nnumber of Sites = " << numSites << "\n";
	std::cout << "number of Particles = " << numParticles << "\n";
	std::cout << "interaction strength U/2J = " << intStrength << "\n\n";


	//getModel(numSites, numParticles);

	//getHamiltonian(numSites, numParticles, intStrength);

	Hamiltonian ham(numSites, numParticles, 1.0);
	ham.getHamiltonianMatrix();
	std::cout << "\ncalculating eigenvalues...\n";
	arma::vec eigval;
	arma::mat eigvec;
	arma::eigs_sym(eigval, eigvec, ham.H, numEig, "sa", tol);
	eigval.print("\nthe smallest eigenvalues are");
	//eigvec.print("\nthe corresponding eigenvectors are");
}


//log the created h5 filename, update numBasis
void getModel(int ns, int np) {
	BHModel bh(ns, np);

	auto start = high_resolution_clock::now();
	bh.mkBasisMatrix();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() / 1000.0 << " seconds has elapsed in generating basis matrix\n\n" << std::endl;
}
