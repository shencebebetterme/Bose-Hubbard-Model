#include "pch.h"
#include "Hamiltonian.h"
#include "BHModel.h"


extern int numSites;
extern int numParticles;
extern int numBasis;
extern double intStrength; //the interaction strength U/2 (with J set to 1)
extern std::string basis_h5_file_name;
extern std::string hamil_bin_file_name;
extern const char* dataset_name;

auto exe_policy = std::execution::par;


Hamiltonian::Hamiltonian() {
	this->loadBasisMat();
	dim = basisMat.n_rows;
	numBasis = dim;
	H = arma::sp_mat(dim, dim);
	T = std::vector<float>(dim, 0.0);
	ind = std::vector<int>(dim, 0);//record the index of T elements
}

void Hamiltonian::loadBasisMat() {
	if (fs::exists(basis_h5_file_name)) {
		basisMat.load(arma::hdf5_name(basis_h5_file_name, dataset_name));
		std::cout << basis_h5_file_name << " ==> " << dataset_name << "\thas been loaded!\n\n";
	}
	else {
		std::cout << basis_h5_file_name << " doesn't exist!\n Exiting.\n\n";
		std::exit(EXIT_FAILURE);
	}
}



//set hamiltonian matrix name
//sparse matrix can't be saved to h5 format
std::string Hamiltonian::bin_name() {
	std::string str = "nS=" + std::to_string(numSites) + "_nP=" + std::to_string(numParticles) + "_hamiltonian.bin";
	return str;
}

void Hamiltonian::saveHamiltonianMatrix() {
	H.save(bin_name(), arma::arma_binary);//use the default dataset name = "dataset"
	//H.save(arma::hdf5_name(hamil_h5_file_name, dataset_name));
	std::cout << "\nHamiltonian matrix has dimension " << dim << " * " << dim << "\n";
	std::cout << "Hamiltonian matrix is successfully created and saved to " << bin_name() << "\n\n";
}


void Hamiltonian::getH0() {
	for (int i = 0; i < dim; i++) {
		//take the i-th row of basisMat
		int iM = 0;
		for (int j = 0; j < numSites; j++) {
			iM += basisMat(i, j) * (basisMat(i, j) - 1);
		}
		H(i, i) = iM * intStrength;
	}
}



//extract the i-th row from the basis vector
basisVecType extractRow(const basisMatType& bm, const int i) {
	basisVecType bv(numSites);
	for (int j = 0; j < numSites; j++) {
		bv[j] = bm(i, j);
	}
	return bv;
}

//calculate hash of a row vector
float calculateHash(const basisVecType& bv) {
	float s = 0;
	for (auto pr = bv.begin(); pr < bv.end(); pr++) {
		s += std::sqrtf(3 + 100 * (pr - bv.begin())) * (*pr);
	}
	return s;
}

#if 1
//calculate the hash of the i-th row of matrix bm
float calculateHash(const basisMatType& bm, int i) {
	basisVecType bv = extractRow(bm, i);
	return calculateHash(bv);
}
#else
//calculate the hash of the i-th row of matrix bm
float calculateHash(const basisMatType& bm, int i) {
	float s = 0;
	int j = 0;
	for (auto jpr = bm.begin_row(i); jpr != bm.end_row(i); jpr++) {
		s += std::sqrtf(3 + 100 * (j)) * (*jpr);
		j++;
	}
}
#endif



template<class ForwardIt, class Ty, class Compare = std::less<>>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const Ty& value, Compare comp = {})
{
	// Note: BOTH type T and the type after ForwardIt is dereferenced 
	// must be implicitly convertible to BOTH Type1 and Type2, used in Compare. 
	// This is stricter than lower_bound requirement (see above)

	first = std::lower_bound(first, last, value, comp);
	return first != last && !comp(value, *first) ? first : last;
}

//determine whether the hopping i->i+1 exist
inline bool existHopping(const basisVecType& bv, int i, int j) {
	//periodic boundary
	if (i >= numSites) { i = i % numSites; }
	if (j >= numSites) { j = j % numSites; }
	//can't hop if there's no particle at site j
	if (bv(j) == 0) return false;
	else return true;
}

//hash value after hopping of j->i
float hashAfterHop(const basisVecType& bv, int i, int j) {
	basisVecType bw = bv;
	//periodic boundary
	if (i >= numSites) { i = i % numSites; }
	if (j >= numSites) { j = j % numSites; }
	bw(i) += 1;
	bw(j) -= 1;
	return calculateHash(bw);
}


void Hamiltonian::getH1() {
	//populate vector T
	for (auto i = T.begin(); i < T.end(); i++) {
		*i = calculateHash(basisMat, i - T.begin());
	}
	//initialize vector ind
	std::iota(ind.begin(), ind.end(), 0);
	//create and release the local variable Tcopy
	{
		std::vector<float> Tcopy(T);
		//sort T
		std::sort(exe_policy, T.begin(), T.end());
		//sort ind
		auto comparator = [&Tcopy](int i1, int i2) {return Tcopy[i1] < Tcopy[i2]; };
		std::sort(exe_policy, ind.begin(), ind.end(), comparator);
	}
	//TODO: parallelize this process
	//find the hopping-connected basis vectors
	for (int iv = 0; iv < dim; iv++) {
		hopping(iv);
	}

}


void Hamiltonian::hopping(int iv) {
	basisVecType bvi = extractRow(basisMat, iv);
	for (int j = 0; j < bvi.n_elem; j++) {
		//if neighboring hopping j+1 -> j exist
		if (existHopping(bvi, j, j + 1)) {
			float hj = hashAfterHop(bvi, j, j + 1);
			//search the hash of the basis vector after hopping
			auto p = binary_find(T.begin(), T.end(), hj);
			//location in sorted vector sort(T)
			int loc_v = p - T.begin();
			//location in original vector T
			int loc_v_ori = ind[loc_v];

			//populate H1 in coordinate (iv, loc_v_ori) and (loc_v_ori, iv)
			//becuase H1 is symmetric
			if (j == bvi.n_elem - 1) {
				H(iv, loc_v_ori) = -std::sqrt((bvi(j) + 1) * bvi(0));
				H(loc_v_ori, iv) = -std::sqrt((bvi(j) + 1) * bvi(0));
			}
			else {
				H(iv, loc_v_ori) = -std::sqrt((bvi(j) + 1) * bvi(j + 1));
				H(loc_v_ori, iv) = -std::sqrt((bvi(j) + 1) * bvi(j + 1));
			}
		}
	}
}