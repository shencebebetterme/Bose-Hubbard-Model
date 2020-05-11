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
	void getHamiltonianMatrix(bool save = false);

	void saveHamiltonianMatrix();
};