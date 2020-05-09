#include "pch.h"
#include "Hamiltonian.h"
#include "BHModel.h"


extern int numSites;
extern int numParticles;
extern int numBasis;
extern double intStrength; //the interaction strength U/2 (with J set to 1)
extern std::string basis_h5_file_name;
extern std::string hamil_bin_file_name;
extern const char* dataset_name;


Hamiltonian::Hamiltonian() {
	this->loadBasisMat();
	dim = basisMat.n_rows;
	numBasis = dim;
	H = arma::sp_mat(dim, dim);
	//H1(d, d);
}

void Hamiltonian::loadBasisMat() {
	if (fs::exists(basis_h5_file_name)) {
		basisMat.load(arma::hdf5_name(basis_h5_file_name, dataset_name));
		std::cout << basis_h5_file_name << " ==> " << dataset_name << "\thas been loaded!\n\n";
	}
	else {
		std::cout << basis_h5_file_name << " doesn't exist!\n Exiting.\n\n";
		std::exit(EXIT_FAILURE);
	}
}

void Hamiltonian::getH0() {
	for (int i = 0; i < dim; i++) {
		//take the i-th row of basisMat
		int iM = 0;
		for (int j = 0; j < numSites; j++) {
			iM += basisMat(i, j)*(basisMat(i, j) - 1);
		}
		H(i, i) = iM * intStrength;
	}
}

void Hamiltonian::getH1() {

}

//set hamiltonian matrix name
//sparse matrix can't be saved to h5 format
std::string Hamiltonian::bin_name() {
	std::string str = "nS=" + std::to_string(numSites) + "_nP=" + std::to_string(numParticles) + "_hamiltonian.bin";
	return str;
}

void Hamiltonian::createHamiltonianMatrix() {
	H.save(bin_name(), arma::arma_binary);//use the default dataset name = "dataset"
	//H.save(arma::hdf5_name(hamil_h5_file_name, dataset_name));
	std::cout << "\nHamiltonian matrix has dimension " << dim << " * " << dim << "\n";
	std::cout << "Hamiltonian matrix is successfully created and saved to " << bin_name() << "\n\n";
}