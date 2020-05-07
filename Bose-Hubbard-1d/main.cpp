
#include "pch.h"
//#include "hdf5.h"
//#include "H5Cpp.h"
//
#include "BHModel.h"

const int numSites = 5;
const int numParticles = 3;
double intStrength; //the interaction strength U (with J set to 1).{


int main() {
	/*
	std::cout << "1d Bose-Hubbard model.\n\nPlease input the number of sites M:" << std::endl;
	std::cin >> numSites;
	std::cout << "\nPlease enter the number of total particles N:" << std::endl;
	std::cin >> numParticles;
	*/
	BHModel bh(numSites, numParticles);
	bh.generateHDF5DataFile();
	bh.readh5dims();

	int nBasis = bh.nBasis;
	int basisLen = bh.basisLen;
	//arma::Mat<int> data = bh.readData();
	//data.print("data:");
	arma::Mat<int> data;
	data.load(arma::hdf5_name("data.h5", "Dataset1"));
	data.print("data matrix is:");
}
