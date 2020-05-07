#pragma once
#include "pch.h"

class BHModel {
public:
	int nSites;
	int nParticles;
	int nBasis = 1; // dim0 in .h5
	int basisLen = 1; // dim1 in .h5
	std::string wl_file_name = "test.wl";
	std::string h5_file_name = "data.h5";

	BHModel(int n_particles = 5, int n_sites = 3) {
		nSites = n_sites;
		nParticles = n_particles;
	}

	void generateHDF5DataFile() {
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
		std::string line3 = "Export[\"" + h5_file_name + "\",Flatten[Permutations/@(PadRight[#,nS]&/@IntegerPartitions[nP,nS]),1]]";
		//std::string cmd1 = "Export[\"data.h5\",Flatten[Permutations/@(PadRight[#,";
		//std::string cmd2 = "]&/@IntegerPartitions[";
		//std::string cmd3 = "]),1]]";
		//std::string cmd_nSites = std::to_string(nSites);
		//std::string cmd_nParticles = std::to_string(nParticles);
		//std::string cmd = cmd1 + cmd_nSites + cmd2 + cmd_nParticles + "," + cmd_nSites + cmd3;
		std::string wl_contents = line1 + line2 + line3;
		std::ofstream wlfile(wl_file_name, std::ios::out);
		wlfile << wl_contents;
		wlfile.close();
		//run .wl file
		std::string cmd = "wolframscript -file " + wl_file_name;
		std::system(cmd.c_str());
	}

	arma::umat readData() {
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
		arma::u64* read_buffer = new arma::u64[nBasis*basisLen];
		auto status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_buffer); //1D vector
		status = H5Dclose(dset);
		status = H5Fclose(file);
		arma::umat basisData(read_buffer, nBasis, basisLen, false);
		delete[] read_buffer;
		
		return basisData;
	}
	// 1. read hdf5 data
	// 2. create armadillo uword matrix from the data
	// 3. all later manipulation devolved to armadillo
};

