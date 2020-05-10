#pragma once
#include "pch.h"
#include "BHModel.h"



class Hamiltonian {
public:
	unsigned int dim;
	arma::sp_mat H;
	//arma::sp_mat H1;//the off diagonal part
	basisMatType basisMat;
	std::vector<float> T;
	std::vector<int> ind;

	//constructor
	Hamiltonian();
	void loadBasisMat();
	
	void getH0();
	void getH1();

	void getH() {
		this->getH0();
		this->getH1();
	}

	//set hamiltonian matrix name
	std::string bin_name();

	void mkHamiltonianMatrix() {
		if (fs::exists(bin_name())) {
			std::cout << "bin file already exists! Exiting.\n";
			return;
		}
		else {
			this->getH();
			saveHamiltonianMatrix();
		}
	}

	void saveHamiltonianMatrix();
};