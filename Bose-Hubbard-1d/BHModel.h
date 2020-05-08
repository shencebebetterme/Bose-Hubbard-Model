#pragma once
#include "pch.h"

const std::string wl_file_name = "generateBasis.wl";
const std::string h5_file_name = "data.h5";
const std::string dim_file_name = "dims.txt";

using basisVecType = arma::Row<int>;

class BHModel {
public:
	int nSites;
	int nParticles;
	int nBasis = 1; // dim0 in .h5
	//int basisLen = 1; // equan to nSites
	basisVecType firstVec;
	basisVecType lastVec;
	basisVecType currentVec;

	//constructor of BHModel class
	BHModel(int n_sites = 1, int n_particles = 1) {
		nSites = n_sites;
		nParticles = n_particles;
		firstVec = arma::zeros<basisVecType>(nSites);
		firstVec(0) = nParticles;
		lastVec = arma::zeros<basisVecType>(nSites);
		lastVec(nSites - 1) = nParticles;
	}

#if 0
	//deprecated
	void generateWlFile() {
		//wolframscript -code 'Export[\"data.h5\",Flatten[Permutations/@(PadRight[#,3]&/@IntegerPartitions[5,3]),1]]'
		/*
		std::string cmd1 = "wolframscript -code 'Export[\"data.h5\",Flatten[Permutations/@(PadRight[#,";
		std::string cmd2 = "]&/@IntegerPartitions[";
		std::string cmd3 = "]),1]]'";
		std::string cmd_nSites = std::to_string(nSites);
		std::string cmd_nParticles = std::to_string(nParticles);
		std::string cmd = cmd1 + cmd_nSites + cmd2 + cmd_nParticles + "," + cmd_nSites + cmd3;
		*/
		//create .wl file
		std::string line1 = "nS=" + std::to_string(nSites) + ";\n";
		std::string line2 = "nP=" + std::to_string(nParticles) + ";\n";
		std::string line3 = "data=Flatten[Permutations/@(PadRight[#,nS]&/@IntegerPartitions[nP,nS]),1];\n";
		std::string line4 = "Export[\"" + h5_file_name + "\",data];\n";
		std::string line5 = "Print/@Dimensions[data];";
		//std::string cmd1 = "Export[\"data.h5\",Flatten[Permutations/@(PadRight[#,";
		//std::string cmd2 = "]&/@IntegerPartitions[";
		//std::string cmd3 = "]),1]]";
		//std::string cmd_nSites = std::to_string(nSites);
		//std::string cmd_nParticles = std::to_string(nParticles);
		//std::string cmd = cmd1 + cmd_nSites + cmd2 + cmd_nParticles + "," + cmd_nSites + cmd3;
		std::string wl_contents = line1 + line2 + line3 + line4 + line5;
		std::ofstream wlfile(wl_file_name, std::ios::out);
		wlfile << wl_contents;
		wlfile.close();
	}
#endif

	//deprecated, call wolframscript in shell
	void MMA_GenerateH5() {
		//run .wl file
		std::string cmd_nSites = " " + std::to_string(nSites);
		std::string cmd_nParticles = " " + std::to_string(nParticles);
		std::string cmd = "wolframscript -file " + wl_file_name + cmd_nSites + cmd_nParticles + " > " + dim_file_name;
		std::system(cmd.c_str());
	}

	//store the h5 matrix dimensions to nBasis and basisLen
	void readh5dims() {
		std::ifstream dimfile(dim_file_name, std::ios::app);
		std::string line;//store temp dim
		std::getline(dimfile, line); nBasis = stoi(line);
		std::getline(dimfile, line); nSites = stoi(line);
	}

	// use algorithm to generate all vectors
	int getK(const basisVecType&) const;
	basisVecType nextVec(const basisVecType&) const;

	void printBasis() {
		currentVec = firstVec;
		while (1) {
			currentVec.print();
			if (sum(currentVec != lastVec)) 
				currentVec = nextVec(currentVec);
			else break;
		}
	}

	const char* h5name() {
		std::string str = "nS=" + std::to_string(nSites) + "_nP=" + std::to_string(nParticles) + "basis.h5";
		return str.c_str();
	}

	void mkBasisMatrix() {
		currentVec = firstVec;
		int nrows = 1;
		arma::Mat<int> basisMat(firstVec);
		while (1) {
			if (sum(currentVec != lastVec)) {
				currentVec = nextVec(currentVec);
				nrows++;
				basisMat = arma::join_cols(basisMat, currentVec);
			}
			else {
				basisMat.save(arma::hdf5_name(h5name(), "dataset"));
				std::cout << "basis matrix has dimension " << nrows << " * " << nSites << "\n";
				std::cout << "basis matrix is successfully created and saved to basis_matrix.h5\n";
				break;
			}
		}
	}

#if 0
	//deprecated
	arma::Mat<int> readData() {
		auto file = H5Fopen(h5_file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		auto dset = H5Dopen2(file, "Dataset1", H5P_DEFAULT);
		//read the dimensions
		auto dspace = H5Dget_space(dset);
		hsize_t dims[2];
		H5Sget_simple_extent_dims(dspace, dims, NULL);
		nBasis = dims[0];
		basisLen = dims[1];
		//move .h5 data to armadillo int matrix
		//arma::umat basisData(nBasis, basisLen);
		//arma::u64* read_buffer = new arma::u64[nBasis*basisLen];
		int* read_buffer = new int[nBasis * basisLen];
		auto status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_buffer); //1D vector
		status = H5Dclose(dset);
		status = H5Fclose(file);
		arma::Mat<int> basisData(read_buffer, nBasis, basisLen, false);
		delete[] read_buffer;

		return basisData;
	}
	// 1. read hdf5 data
	// 2. create armadillo uword matrix from the data
	// 3. all later manipulation devolved to armadillo
#endif
};


int BHModel::getK(const basisVecType& vec) const {
	int last = vec.n_cols - 1;
	if (vec.n_cols == 1) return 1;
	for (int i = last - 1; i >= 0; ) {
		if (vec(i) == 0) i--;
		else { 
			//std::cout << "k=" << i << "\n"; 
			return i; 
		}
	}
}

basisVecType BHModel::nextVec(const basisVecType& vec) const {
	int last = vec.n_cols - 1;
	basisVecType res = vec;
	int k = getK(res);
	res(k) -= 1;
	res(k + 1) = vec(last) + 1;
	if (last > (k + 1)) { res(last) = 0; }
	return res;
}