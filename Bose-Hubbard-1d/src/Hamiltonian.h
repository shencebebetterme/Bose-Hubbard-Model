#pragma once
#include "pch.h"
#include "BHModel.h"



class Hamiltonian {
public:
	unsigned int dim;
	arma::sp_mat H;
	//arma::sp_mat H1;//the off diagonal part
	basisMatType basisMat;

	//constructor
	Hamiltonian();
	void loadBasisMat();
	
	void getH0();
	void getH1();

	//set hamiltonian matrix name
	std::string bin_name();


	void createHamiltonianMatrix();
};