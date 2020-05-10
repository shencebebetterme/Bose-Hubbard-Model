#pragma once
#include "pch.h"
#include "BHModel.h"



class Hamiltonian {
public:
	int nSites;
	int nParticles;
	double intStr;
	int nBasis = 1;
	unsigned int dim;
	arma::sp_mat H;
	//arma::sp_mat H1;//the off diagonal part
	basisMatType basisMat;
	std::vector<float> T;
	std::vector<int> ind;

	//constructor
	Hamiltonian(int, int, double);
	void loadBasisMat();
	
	void getH0();
	void getH1();

	void calculateH() {
		this->getH0();
		this->getH1();
		float density = 100.0 * H.n_nonzero / H.n_elem;
		std::cout << "\nHamiltonian matrix calculated!\n";
		std::cout << "Hamiltonian matrix has dimension " << dim << " * " << dim << "\n"
			<< "Hamiltonian matrix has density " << density << "%\n\n";
	}


	basisVecType extractRow(const basisMatType& bm, const int i); 
	float calculateHash(const basisVecType& bv);
	float calculateHash(const basisMatType& bm, int i); 
	float hashAfterHop(const basisVecType& bv, int i, int j);
	inline bool existHopping(const basisVecType& bv, int i, int j);

	//set hamiltonian matrix name
	static std::string bin_name(int, int, double);

	//if .bin file exists, load it from disk
	//otherwise create the Hamiltonian matrix and save it to disk
	void getHamiltonianMatrix() {
		std::string ham_name = bin_name(nSites, nParticles, intStr);
		if (fs::exists(ham_name)) {
			//std::cout << ham_name << " file already exists!\n";
			H.load(ham_name, arma::arma_binary);
			dim = H.n_rows;
			float density = 100.0 * H.n_nonzero / H.n_elem;
			std::cout << ham_name << "\thas been loaded from disk\n"
				<< "Hamiltonian matrix has dimension " << dim << " * " << dim << "\n"
				<< "Hamiltonian matrix has density " << density << "%\n\n";
			return;
		}
		else {
			std::cout << "\ncalculating Hamiltonian matrix...\n";
			this->calculateH();
			saveHamiltonianMatrix();
		}
	}

	void saveHamiltonianMatrix();
};