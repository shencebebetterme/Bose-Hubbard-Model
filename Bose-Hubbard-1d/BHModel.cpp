#include "pch.h"
#include "BHModel.h"




void BHModel::printBasis() {
	currentVec = firstVec;
	while (1) {
		currentVec.print();
		if (sum(currentVec != lastVec))
			currentVec = nextVec(currentVec);
		else break;
	}
}



int BHModel::getK(const basisVecType& vec) const {
	int last = vec.n_elem - 1;
	if (vec.n_elem == 1) return 1;
	for (int i = last - 1; i >= 0; ) {
		if (vec(i) == 0) i--;
		else {
			//std::cout << "k=" << i << "\n"; 
			return i;
		}
	}
}

basisVecType BHModel::nextVec(const basisVecType& vec) const {
	int last = vec.n_elem - 1;
	basisVecType res = vec;
	int k = getK(res);
	res(k) -= 1;
	res(k + 1) = vec(last) + 1;
	if (last > (k + 1)) { res(last) = 0; }
	return res;
}



// get the total number of basis, save the value to nBasis
// store the basis vectors into a std::vector
void BHModel::getNBasis() {
	currentVec = firstVec;
	matVec.push_back(currentVec);
	int nrows = 1;
	while (1) {
		if (sum(currentVec != lastVec)) {
			currentVec = nextVec(currentVec);
			//basisMat.insert_rows(nrows, currentVec);

			matVec.push_back(currentVec);
			//show progress
			if (nrows % 10000 == 1) {
				std::cout << nrows << "\tvectors generated!\n";
			}
			nrows++;
		}
		else {
			nBasis = nrows;
			std::cout << "\nin total " << nBasis << " basis vectors generated\n";
			break;
		}
	}
}

void BHModel::createBasisMatrix() {
	arma::Mat<uint8_t> basisMat(nBasis, nSites);
	for (int i = 0; i < nBasis; i++) {
		for (int j = 0; j < nSites; j++) {
			basisMat[i, j] = matVec[i](j);
		}
	}
	basisMat.save(arma::hdf5_name(h5name(), "dataset"));
	std::cout << "\nbasis matrix has dimension " << nBasis << " * " << nSites << "\n";
	std::cout << "basis matrix is successfully created and saved to " << h5name() << "\n\n";
}