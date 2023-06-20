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
			rank==0,
			params.initTrans());

	hacc::Distribution *dist = new hacc::Distribution(MPI_COMM_WORLD, params.ng());
	
	Wavefunction wf(params, dist);

	long Ng = params.ng()*params.ng()*params.ng();

	initialize_wavefunction(Ng, CC, TC, params.rL(), params.zin(), params.iseed(), *dist, wf);

	MPI_Barrier(MPI_COMM_WORLD);

	Timer::stopTimerStats(t_init);
	if(rank==0) printf("\n");

	SolverDiscrete *solver = new SolverDiscrete(*dist);

	std::ostringstream inter1;
	inter1.str("intermediate_phase");

	wf.writeWavefunction(params.outBase(), -1);
	MPI_Barrier(MPI_COMM_WORLD);
	
	std::ostringstream tempname;
	tempname<<inter1.str()<<".step.-1";
	
	wf.writePhase(tempname.str());
	//this is the initial wave function with zero velocity

	MPI_Barrier(MPI_COMM_WORLD);

	std::ostringstream temp;
	temp<<params.outBase()<<".step.ini";

	for(int step=0; step<params.nsteps(); step++){
			
		std::ostringstream filename;
		std::ostringstream outfile;
		filename<<inter1.str()<<".step."<<step-1;	
		outfile<<inter1.str()<<".step."<<step;	
		//////////////////////SETUP///////////////////////////
		wf.readWavefunction(temp.str()); //read the zero vel wave function 
		if(rank==0){ std::cout<<"The following should reflect 0s "<<std::endl; }
		wf.check_phase(step);
		MPI_Barrier(MPI_COMM_WORLD);
		wf.readPhase(filename.str()); //replace phase with phase from previous step

		wf.write1D(params.outBase(), step);

		Timer::startTimer(t_step);	
		
		if(rank==0){ std::cout<<"Iteration: "<<step<<std::endl; }
		wf.check_phase(step);

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
		
		ts.set_step(); //sets tau to the default (afin-ain)/nsteps
		
		if(rank==0){ std::cout<<"a: "<<ts.aa()<<std::endl; }
		MPI_Barrier(MPI_COMM_WORLD);
		wf.forward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		//////////////////////////////////////////////////////
		
		///////////////////////STREAM/////////////////////////	
		wf.stream(ts.tau2(), ts.aa(), ts.adot());
		ts.advanceHalfStep();
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.backward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		/////////////////////////////////////////////////////	
			
		////////////////////////KICK/////////////////////////	
		
		//compute physical rho
		//a needed to get proper density
		wf.computeRho(1.0);	

		MPI_Barrier(MPI_COMM_WORLD);
		solver->solve(wf.m_rho.data(), wf.m_phi.data());
		MPI_Barrier(MPI_COMM_WORLD);
		
		wf.kick(ts.tau(), ts.aa(), ts.adot(), ts.phiscal());
		ts.advanceHalfStep();
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.forward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		/////////////////////////////////////////////////////	
		

		///////////////////////STREAM/////////////////////////	
	
		wf.stream(ts.tau2(), ts.aa(), ts.adot());
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.backward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		
		wf.writePhase(outfile.str());
		MPI_Barrier(MPI_COMM_WORLD);
		

		/////////////////////////////////////////////////////	

		Timer::stopTimerStats(t_step);
	}

	Timer::stopTimerStats(t_total);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
