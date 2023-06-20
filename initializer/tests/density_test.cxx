#include "CosmoClass.h"
#include "initread.h"
#include "Initializer.h"
#include "Distribution.hpp"
#include "Wavefunction.h"
#include "GenericIO.h"

#include <iostream>
#include <string.h>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
	MPI::Init();

	int rank;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int argvi = optind;

	std::string pName;
	pName = argv[argvi];

	Parameters params(pName);

	std::string paramout = params.getParamsString();

	if(rank == 0){
		std::cout << paramout << std::endl;
	}

	CosmoClass CC(params.omega_matter(),
			params.omega_cdm(),
			params.omega_baryon(),
			params.omega_cb(),
			params.omega_nu(),
			params.f_nu_massless(),
			params.f_nu_massive(),
			params.omega_radiation(),
			params.hubble(),
			params.w_de(),
			params.wa_de(),
			params.cdm_mass(),
			params.scatter());
	
	CosmoClass::CosmoClassParams Cosmology = CC.getCosmoClassParams();

	if(rank==0){	
		std::cout<<"kappa1 coeff for this cosmology: "<<std::scientific<<Cosmology.k1coeff<<std::endl;
		std::cout<<"kappa2 coeff for this cosmology: "<<std::scientific<<Cosmology.k2coeff<<std::endl;
	}
	
	TransferClass TC(Cosmology, 
			params.ns(),
			params.ss8(),
			params.trans(),
			params.inTransfer(),
			rank==0);

	hacc::Distribution *dist = new hacc::Distribution(MPI_COMM_WORLD, params.ng());
	
	Wavefunction wf(params, dist);
	
	long Ng = params.ng()*params.ng()*params.ng();

	initialize_wavefunction(Ng, CC, TC, params.rL(), params.zin(), params.iseed(), *dist, wf);

	MPI_Barrier(MPI_COMM_WORLD);

	wf.writeWavefunction(params.outBase(),0);
	std::cout<<"here"<<std::endl;

	return 1;
}
