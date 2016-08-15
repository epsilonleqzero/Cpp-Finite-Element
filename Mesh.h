/*
 * Mesh.h
 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include <vector>

#ifndef MESH_H_
#define MESH_H_

class Mesh {
public:
	Mesh();
	Mesh(std::vector<double> meshprops);
	virtual ~Mesh();
	arma::sp_mat stiffness;
	arma::vec area;
	arma::mat node;
	arma::umat elem;
	arma::uvec bdNode;
	arma::uvec isbdNode;
	arma::uvec freeNode;
	arma::uword N;
private:

	arma::sp_mat assembleMatrix();
	void makeMesh(arma::vec xr,arma::vec yr);
	void findBoundary();
	double h;
	double n;
	arma::uword NT;
	};

#endif /* MESH_H_ */
