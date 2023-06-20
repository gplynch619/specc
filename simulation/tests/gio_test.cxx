#include "CosmoClass.h"
#include "initread.h"
#include "Initializer.h"
#include "Distribution.hpp"
#include "Wavefunction.h"
#include "GenericIO.h"
#include "TimeStepper.h"
#include "poisson.h"
#include "Timer.h"

#include "sc.h"

#include <iostream>
#include <string.h>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#define DIMENSION 3
#define EPS 1e-8

int main(int argc, char *argv[]){
	int step0 = 0;

	Timer::TimerRef t_total = Timer::getTimer("total");
	Timer::TimerRef t_init = Timer::getTimer("init");
	Timer::TimerRef t_step = Timer::getTimer("step");
	
	MPI::Init();

	int rank;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int argvi = optind;

	std::string pName;
	pName = argv[argvi];

	Parameters params(pName);

	Timer::startTimer(t_total);

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
			rank==0,
			params.initTrans());

	hacc::Distribution *dist = new hacc::Distribution(MPI_COMM_WORLD, params.ng());
	
	Wavefunction wf(params, dist);

	long Ng = params.ng()*params.ng()*params.ng();

	std::vector<complex_p> before_save;
	std::vector<complex_p> after_save;
	long int my_Ng;


	initialize_wavefunction(Ng, CC, TC, params.rL(), params.zin(), params.iseed(), *dist, wf);

	MPI_Barrier(MPI_COMM_WORLD);

	my_Ng = wf.Ng();

	before_save.resize(my_Ng);
	after_save.resize(my_Ng);

	for(int i=0; i<my_Ng; i++){
		before_save.at(i) = wf.m_psi.at(i);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0) printf("\n");

	wf.writeWavefunction(params.outBase(), -1);
	MPI_Barrier(MPI_COMM_WORLD);

	std::ostringstream temp;
	temp<<params.outBase()<<".step.ini";
	
	wf.readWavefunction(temp.str());

	for(int i=0; i<my_Ng; i++){
		after_save.at(i) = wf.m_psi.at(i);
	}
	
	int neq_count=0;

	for(int i=0; i<my_Ng; i++){
		if((after_save.at(i).mag != before_save.at(i).mag)
				&&(after_save.at(i).phase != before_save.at(i).phase)){
			neq_count++;
		}
	}

	std::cout<<"Rank "<<rank<<": "<<neq_count<<std::endl;
		
	
	Timer::stopTimerStats(t_total);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
