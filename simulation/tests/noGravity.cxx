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

int main(int argc, char *argv[]){
	int step0 = 0;

	Timer::TimerRef t_total = Timer::getTimer("total");
	Timer::TimerRef t_init = Timer::getTimer("init");
	
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
	Timer::startTimer(t_init);

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

	Timer::stopTimerStats(t_init);
	if(rank==0) printf("\n");

	SolverDiscrete *solver = new SolverDiscrete(*dist);

	TimeStepper ts(params.alpha(),
			params.ain(),
			params.afin(),
			params.nsteps(),
			params.omega_matter(),
			params.omega_cdm(),
			params.omega_baryon(),
			params.omega_cb(),
			params.omega_nu(),
			params.omega_radiation(),
			params.f_nu_massless(),
			params.f_nu_massive(),
			params.w_de(),
			params.wa_de(),
			params.zpr(),
			params.nzpr());

	MPI_Barrier(MPI_COMM_WORLD);
	solver->forward_solve(wf.m_rho.data());
	MPI_Barrier(MPI_COMM_WORLD);

	writePk(solver, params, params.outBase(), "ini");
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	wf.writeWavefunction(params.outBase(), -1);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){ std::cout<<"tau: "<<ts.tau()<<std::endl;  }
	if(rank==0){ std::cout<<"tau2: "<<ts.tau2()<<std::endl;  }

	for(int step=0; step<params.nsteps(); step++){
		if(rank==0){ std::cout<<"TIMESTEP "<<step<<std::endl;  }

		/////////////////////////////////////////////////////	
		MPI_Barrier(MPI_COMM_WORLD);
		wf.forward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		/////////////////////////////////////////////////////	
		
	
		///////////////////////STREAM/////////////////////////	
		wf.stream(ts.tau(), ts.aa(), ts.adot());
		ts.advanceFullStep();
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.backward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		/////////////////////////////////////////////////////	
		wf.computeRho(ts.aa());
		
		MPI_Barrier(MPI_COMM_WORLD);
		solver->forward_solve(wf.m_rho.data());
		MPI_Barrier(MPI_COMM_WORLD);
	
		std::ostringstream ost;
		ost << step <<".z"<<ts.zz();
		
		if(step%10==0){ writePk(solver, params, params.outBase(), ost.str()); }
		if(step%10==0){ wf.writeWavefunction(params.outBase(), step); }
		
		if(rank==0){ std::cout<<"Step "<<step<<" done!"<<std::endl; }
	}

	Timer::stopTimer(t_total);

	return 1;
}
