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

	void getH() {
		this->getH0();
		this->getH1();
		std::cout << "\nHamiltonian matrix obtained!\n\n";
	}


	basisVecType extractRow(const basisMatType& bm, const int i); 
	float calculateHash(const basisVecType& bv);
	float calculateHash(const basisMatType& bm, int i); 
	float hashAfterHop(const basisVecType& bv, int i, int j);
	inline bool existHopping(const basisVecType& bv, int i, int j);
	//set hamiltonian matrix name
	static std::string bin_name(int, int, double);

	void mkHamiltonianMatrix() {
		if (fs::exists(bin_name(nSites, nParticles, intStr))) {
			std::cout << bin_name(nSites, nParticles, intStr) << " file already exists!\n";
			return;
		}
		else {
			std::cout << "\ncreating Hamiltonian matrix\n";
			this->getH();
			saveHamiltonianMatrix();
		}
	}

	void saveHamiltonianMatrix();
};