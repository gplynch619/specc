#include "initread.h"
#include "Distribution.hpp"
#include "GenericIO.h"
#include "poisson.h"

#include "sc.h"

#include <iostream>
#include <string.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
	
	MPI::Init();

	int rank;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int argvi = optind;

	std::string sng;
	sng = argv[argvi];

	int ng = atoi(sng.c_str());
	
	hacc::Distribution *dist = new hacc::Distribution(MPI_COMM_WORLD, ng);

	SolverDiscrete *solver = new SolverDiscrete(*dist);

	std::vector<complex_t> inarr;
	std::vector<complex_t> phi;
	
	std::ostringstream ost;

	int nelem;

	ost << "density.debug.rspace.gio";
	{
		gio::GenericIO GIO(dist->parent_comm(), ost.str());
		GIO.openAndReadHeader();

		std::cout<<"test"<<std::endl;
		nelem = GIO.readNumElems();
		inarr.resize(nelem);
		std::vector<double> real;
		std::vector<double> imag;
		
		real.resize(nelem+GIO.requestedExtraSpace()/sizeof(double));
		imag.resize(nelem+GIO.requestedExtraSpace()/sizeof(double));

		GIO.addVariable("real", real.data(), gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("imag", imag.data(), gio::GenericIO::VarHasExtraSpace);

		GIO.readData();
		std::cout<<"test"<<std::endl;
		for(long i=0; i<nelem; ++i){
			inarr.at(i).real(real[i]);
			inarr.at(i).imag(imag[i]);
		}
	}

	ost.str("");

	phi.resize(inarr.size());

	solver->solve(inarr.data(), phi.data());

	MPI_Barrier(MPI_COMM_WORLD);

	ost<<"potential.specc.gio";

  	{
		gio::GenericIO GIO(dist->parent_comm(), ost.str());
		
		std::vector<double> real;
		std::vector<double> imag;

		for(long i=0; i<nelem; ++i){
			real.push_back(phi[i].real());
			imag.push_back(phi[i].imag());
		}

		GIO.setNumElems(nelem);

		real.resize(nelem+GIO.requestedExtraSpace()/sizeof(double));
		imag.resize(nelem+GIO.requestedExtraSpace()/sizeof(double));

		GIO.addVariable("real", real.data(), gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("imag", imag.data(), gio::GenericIO::VarHasExtraSpace);

		GIO.write();
  	}

	return 1;
}
