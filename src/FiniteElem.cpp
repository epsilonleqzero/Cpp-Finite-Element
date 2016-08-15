/*
 * FiniteElem.cpp
 *
 * This class runs the finite element method and obtains a
 * numerical solution to a PDE on a domain given into the constructor.
 *
 * Currently, only two functions have been implemented, and square domains.
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FiniteElem.h"
#include "FunZeros.h"
#include "FunOnes.h"
#include "FunSC.h"
#include "FunSCbD.h"

using namespace std;
using namespace arma;

FiniteElem::FiniteElem() {

}

/**
 * Constructor to create and run.
 *
 * This constructor creates and runs the finite element method using the properties
 * provided in meshprops and using the name provided in fun.
 *
 * @param meshprops - Vector containing the information on the square.
 * @param fun - string containing the name of the function to use.
 *
 */
FiniteElem::FiniteElem(vector<double> meshprops,string fun) {
	// Create mesh.
	Mesh msh(meshprops);

	// Choose and apply function using polymorphism.
	if(fun.find("Bump") != string::npos ||
			fun.find("bump") != string::npos ){
		pde=new FunOnes;
		bdfun=new FunZeros;
	} else{
		pde=new FunSC;
		bdfun=new FunSCbD;
	}
	// Set the mesh variable to the mesh created.
	mesh=msh;
	uword N=mesh.node.n_rows;
	// Initialize solution vector.
	u=zeros<vec>(N);
	// Calculate the right hand side of the equation.
	vec b=calcRHS();
	mat node=mesh.node;
	// Calculate the boundary condition.
	u.rows(mesh.bdNode)=bdfun->evalF(node.rows(mesh.bdNode));
	// Release the function classes created.
	delete pde;
	delete bdfun;
	// Setup system to solve.
	vec r=b-(mesh.stiffness*u);
	mat A(mesh.stiffness);
	mat Af=A.submat(mesh.freeNode,mesh.freeNode);
	// Solve the linear system Au=r
	u.rows(mesh.freeNode)=solve(Af,r.rows(mesh.freeNode));
}

/**
 * calcRHS() is an internal method to find the right hand side.
 *
 * This is used above to set b for the system Au=b
 */
vec FiniteElem::calcRHS(){
	// Get needed Matrices and vectors.
	vec area=mesh.area; mat node=mesh.node; umat elem=mesh.elem;
	// Calculate Midpoints.
	mat mid1=(node.rows(elem.col(1))+node.rows(elem.col(2)))/2.0;
	mat mid2=(node.rows(elem.col(2))+node.rows(elem.col(0)))/2.0;
	mat mid3=(node.rows(elem.col(0))+node.rows(elem.col(1)))/2.0;
	// Calculate vectors to accumulate
	vec bt1=(area%((pde->evalF(mid2))+(pde->evalF(mid3))))/6.0;
	vec bt2=(area%((pde->evalF(mid3))+(pde->evalF(mid1))))/6.0;
	vec bt3=(area%((pde->evalF(mid1))+(pde->evalF(mid2))))/6.0;
	// Vertically concatenate vectors.
	vec bt12=join_vert(bt1,bt2);
	vec bts=join_vert(bt12,bt3);
	uvec elemv=vectorise(elem);
	// Accumulate the array.
	vec ret=accumArray(elemv,bts,node.n_rows);
	return ret;
}

/**
 * accumArray implements the accumarray method from MATLAB.
 *
 * This is a direct implementation and does the exact same thing
 * as the MATLAB function.
 *
 * @param subs - vector of indices.
 * @param ar - vector of values.
 * @param N - size of the array to be made.
 * @return accumulated array (same as MATLAB would).
 */
vec FiniteElem::accumArray(uvec subs,vec ar,uword N){

	vec S=zeros<vec>(N);

	for(uword i=0;i<N;i++){
		// Get the subscripts.
		uvec q1=find(subs ==(i));
		vec spot;
		double thesum=0.0;
		if(!q1.is_empty()){
			// Find elements at indices q1
			spot=ar.elem(q1);
			// Sum elements.
			thesum=sum(spot);
			// Set at position i.
			S(i)=thesum;
		}
		// If it is empty, add 0.
		else{
			S(i)=0.0;
		}
	}

	return S;
}

FiniteElem::~FiniteElem() {
}

