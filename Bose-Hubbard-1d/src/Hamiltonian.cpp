#include "pch.h"
#include "Hamiltonian.h"
#include "BHModel.h"



//auto exe_policy = std::execution::par_unseq;
auto exe_policy = std::execution::par;



Hamiltonian::Hamiltonian(int ns, int np, double intstr) {
	nSites = ns;
	nParticles = np;
	intStr = intstr;
	this->loadBasisMat();
	dim = basisMat.n_rows;
	nBasis = dim;
	H = arma::sp_mat(dim, dim);
	T = std::vector<float>(dim, 0.0);
	ind = std::vector<int>(dim, 0);//record the index of T elements
}

void Hamiltonian::loadBasisMat() {
	std::string h5_name = BHModel::h5name(nSites, nParticles);
	if (fs::exists(h5_name)) {
		basisMat.load(arma::hdf5_name(h5_name, "dataset"));
		std::cout << h5_name << "\thas been loaded!\n\n";
	}
	else {
		std::cout << h5_name << " doesn't exist!\t try creating " << h5_name << "\n\n";
		BHModel bh(nSites, nParticles);
		bh.mkBasisMatrix();
		basisMat.load(arma::hdf5_name(h5_name, "dataset"));
		std::cout << h5_name << "\thas been loaded!\n\n";
		//std::exit(EXIT_FAILURE);
	}
}



//set hamiltonian matrix name
//sparse matrix can't be saved to h5 format
std::string Hamiltonian::bin_name(int ns, int np, double intstr) {
	std::string str = "nS=" + std::to_string(ns) + "_nP=" + std::to_string(np) + "_intstr=" + std::to_string(intstr) + "_hamiltonian.bin";
	return str;
}

void Hamiltonian::saveHamiltonianMatrix() {
	H.save(bin_name(nSites, nParticles, intStr), arma::arma_binary);//use the default dataset name = "dataset"
	//H.save(arma::hdf5_name(hamil_h5_file_name, dataset_name));
	std::cout << "\nHamiltonian matrix has dimension " << dim << " * " << dim << "\n";
	std::cout << "Hamiltonian matrix is successfully created and saved to " << bin_name(nSites, nParticles, intStr) << "\n\n";
}


void Hamiltonian::getH0() {
	for (int i = 0; i < dim; i++) {
		//take the i-th row of basisMat
		int iM = 0;
		for (int j = 0; j < nSites; j++) {
			iM += basisMat(i, j) * (basisMat(i, j) - 1);
		}
		H(i, i) = iM * intStr;
	}
}



//extract the i-th row from the basis vector
basisVecType Hamiltonian::extractRow(const basisMatType& bm, const int i) {
	basisVecType bv(nSites);
	for (int j = 0; j < nSites; j++) {
		bv[j] = bm(i, j);
	}
	return bv;
}

//calculate hash of a row vector
float Hamiltonian::calculateHash(const basisVecType& bv) {
	float s = 0;
	for (auto pr = bv.begin(); pr < bv.end(); pr++) {
		s += std::sqrtf(3 + 100 * (pr - bv.begin())) * (*pr);
	}
	return s;
}

#if 1
//calculate the hash of the i-th row of matrix bm
float Hamiltonian::calculateHash(const basisMatType& bm, int i) {
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

//determine whether the hopping j->i exist
inline bool Hamiltonian::existHopping(const basisVecType& bv, int i, int j) {
	//periodic boundary
	if (i >= nSites) { i = i % nSites; }
	if (j >= nSites) { j = j % nSites; }
	//can't hop if there's no particle on site j
	if (bv(j) == 0) return false;
	else return true;
}


//hash value after hopping of j->i
float Hamiltonian::hashAfterHop(const basisVecType& bv, int i, int j) {
	basisVecType bw = bv;
	//periodic boundary
	if (i >= nSites) { i = i % nSites; }
	if (j >= nSites) { j = j % nSites; }
	bw(i) += 1;
	bw(j) -= 1;
	return calculateHash(bw);
}


void Hamiltonian::getH1() {
	//store the hash values in T
	for (auto i = T.begin(); i < T.end(); i++) {
		*i = calculateHash(basisMat, i - T.begin());
	}
	//initialize vector by sequence 0, 1, 2, ...
	std::iota(ind.begin(), ind.end(), 0);

	std::vector<int> indc(ind);//create the local copy
	std::vector<float> Tcopy(T);//create the local variable Tcopy
	//sort T
	std::sort(exe_policy, T.begin(), T.end());
	//sort ind
	auto comparator = [&Tcopy](int i1, int i2) {return Tcopy[i1] < Tcopy[i2]; };
	std::sort(exe_policy, ind.begin(), ind.end(), comparator);
	
	
	//find the hopping-connected basis vectors
	for (int iv = 0; iv < dim; iv++) {
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
	//parallelize this process
	//for_each(exe_policy, indc.begin(), indc.end(), [&](int& n) {hopping(n); });

}
