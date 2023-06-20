#include "CosmoClass.h"
#include "initread.h"
#include "sctypes.h"
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
	
	TransferClass TC(Cosmology, 
			params.ns(),
			params.ss8(),
			params.trans(),
			params.inTransfer(),
			rank==0);

	if(rank==0){	
		std::cout<<"kappa1 coeff for this cosmology: "<<std::scientific<<Cosmology.k1coeff<<std::endl;
		std::cout<<"kappa2 coeff for this cosmology: "<<std::scientific<<Cosmology.k2coeff<<std::endl;
	}
	REAL k[1000];
	REAL growth[1000];
	REAL deriv[1000];

	for(int i=0; i<1000; ++i){
		REAL kk =0.05+0.00994999*i;
		REAL g;
		REAL gdot;
		k[i] = kk;
		CC.GrowthFactor(200.0, kk, g, gdot);
		//CC.GrowthFactorLCDM(200.0, &g, &gdot);
		growth[i]=g;
		deriv[i]=gdot;
	}

	std::string outBase;
	outBase = params.outBase();

	std::string fname;
	FILE *pfile;
	std::ofstream OutFile;
	int i;
	fname=outBase+".gf_test.dat";
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<1000; ++i){
		OutFile<<k[i]<<" "<<growth[i]<<std::endl;
	}
	OutFile.close();

	fname=outBase+".gdot_test.dat";
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<1000; ++i){
		OutFile<<k[i]<<" "<<deriv[i]<<std::endl;
	}
	OutFile.close();
	
	return 1;
}
