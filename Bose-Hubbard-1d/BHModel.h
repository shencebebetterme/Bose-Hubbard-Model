#pragma once
#include "pch.h"

class BHModel {
public:
	int nSites;
	int nParticles;
	int nBasis = 1;

	BHModel(int n_particles = 5, int n_sites = 3) {
		nSites = n_sites;
		nParticles = n_particles;
	}

	void generateHDF5DataFile() {
		//wolframscript -code 'Export[\"data.h5\",Flatten[Permutations/@(PadRight[#,3]&/@IntegerPartitions[5,3]),1]]'
		std::string cmd1 = "wolframscript -code 'Export[\"data.h5\",Flatten[Permutations/@(PadRight[#,";
		std::string cmd2 = "]&/@IntegerPartitions[";
		std::string cmd3 = "]),1]]'";
		std::string cmd_nSites = std::to_string(nSites);
		std::string cmd_nParticles = std::to_string(nParticles);
		std::string cmd = cmd1 + cmd_nSites + cmd2 + cmd_nParticles + "," + cmd_nSites + cmd3;
		std::system(cmd.c_str());
	}

	// 1. read hdf5 data
	// 2. create armadillo uword matrix from the data
	// 3. all later manipulation devolved to armadillo
};

