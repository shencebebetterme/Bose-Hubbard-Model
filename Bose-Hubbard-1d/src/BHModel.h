#pragma once
#include "pch.h"


using basisVecType = arma::Row<uint8_t>;
using basisMatType = arma::Mat<uint8_t>;
namespace fs = std::filesystem;

class BHModel {
public:
	int nSites;
	int nParticles;
	int nBasis = 1; // dim0 in .h5
	//int basisLen = 1; // equan to nSites
	basisVecType firstVec;
	basisVecType lastVec;
	basisVecType currentVec;
	std::vector<basisVecType> matVec = {};

	//constructor of BHModel class
	BHModel(int n_sites = 1, int n_particles = 1) {
		nSites = n_sites;
		nParticles = n_particles;
		firstVec = arma::zeros<basisVecType>(nSites);
		firstVec(0) = nParticles;
		lastVec = arma::zeros<basisVecType>(nSites);
		lastVec(nSites - 1) = nParticles;
	}

	// use algorithm to generate all vectors
	int getK(const basisVecType&) const;
	basisVecType nextVec(const basisVecType&) const;

	void createBasisMatrix();
	void printBasis();
	void getNBasis();
	

	static std::string h5name(int ns, int np) {
		std::string str = "nS=" + std::to_string(ns) + "_nP=" + std::to_string(np) + "_basis.h5";
		return str;
	}

	void mkBasisMatrix() {
		if (fs::exists(h5name(nSites, nParticles))) { 
			std::cout << h5name(nSites, nParticles) << " already exists!.\n"; return;
		}
		else {
			getNBasis();
			createBasisMatrix();
		}
	}

};


