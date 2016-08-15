/*
 * FiniteElem.h
 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */
#include <armadillo>
#include <string>
#include <vector>
#include "Mesh.h"
#include "PdeFun.h"

#ifndef FINITEELEM_H_
#define FINITEELEM_H_

class FiniteElem {
public:
	FiniteElem();
	FiniteElem(std::vector<double> meshprops,std::string fun);
	virtual ~FiniteElem();
	Mesh mesh;
	arma::vec u;
private:
	arma::vec calcRHS();
	arma::vec accumArray(arma::uvec subs,arma::vec ar,arma::uword N);
	PdeFun * pde;
	PdeFun * bdfun;
};

#endif /* FINITEELEM_H_ */
