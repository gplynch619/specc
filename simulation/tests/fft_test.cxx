#include "Distribution.hpp"
#include "GenericIO.h"
#include "poisson.h"

#include "sc.h"
#include "complex-type.h"


#include <iostream>
#include <string.h>
#include <string>
#include <iomanip>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
	MPI::Init();

	int rank;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int ng = 8;
	int Ng = ng*ng*ng; 

	complex_t zero(0.0,0.0);
	complex_t one(1.0, 0.0);

	std::vector<complex_t> fxi;
	std::vector<complex_t> fx;
	fx.resize(Ng);
	fxi.resize(Ng);
	int i=0;
	
	for(int x=0; x<ng; ++x){
		for(int y=0; y<ng; ++y){
			for(int z=0; z<ng; ++z){
				//fx.at(i).real(x*x + y*y + z*z);
				//fxi.at(i).real(x*x + y*y + z*z);
				if(x==0 && y==0 && z==0){
					fx.at(i) = one;
					fxi.at(i) = one;
				} else{
					fx.at(i) = zero;
					fxi.at(i) = zero;
				}
				i++;
			}
		}
	}

	hacc::Distribution *dist = new hacc::Distribution(MPI_COMM_WORLD, ng);
	MPI_Barrier(MPI_COMM_WORLD);

	SolverDiscrete *solver = new SolverDiscrete(*dist);
	
	MPI_Barrier(MPI_COMM_WORLD);
	solver->forward_solve(fx.data());
	MPI_Barrier(MPI_COMM_WORLD);

	solver->backward_solve_only(fx.data());
	MPI_Barrier(MPI_COMM_WORLD);
/*
	double scal = pow(ng, -3.0);
	for(int j=0; j<Ng; ++j){
		fx.at(j) *= scal;
	}
*/
	i=0;
	
	std::string fname;
	std::ofstream OutFile;
	fname="fft.out";
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(int x=0; x<ng; ++x){
		for(int y=0; y<ng; ++y){
			for(int z=0; z<ng; ++z){
				OutFile << x <<" "<<y<<" "<<z<<" "<<fxi.at(i)<<" "<<fx.at(i)<<std::endl;
				i++;
			}
		}
	}

	OutFile.close();
	return 1;
}
