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
	
	hacc::Distribution *dist = new hacc::Distribution(MPI_COMM_WORLD, params.ng());
	
	Wavefunction wf(params, dist);

	initialize_planewave(params.kmode(),  params.rL(), params.iseed(), *dist, wf);

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
	
	wf.computeRho(1.0);
	
	//wf.writeWavefunction(params.outBase(), -1);
	wf.write1D(params.outBase(), -1);


	MPI_Barrier(MPI_COMM_WORLD);
	solver->forward_solve(wf.m_rho.data());
	MPI_Barrier(MPI_COMM_WORLD);

	writePk(solver, params, params.outBase(), "ini");
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){ std::cout<<"tau: "<<ts.tau()<<std::endl;  }
	if(rank==0){ std::cout<<"tau2: "<<ts.tau2()<<std::endl;  }
	if(rank==0){ std::cout<<"dt: "<<ts.dt()<<std::endl;  }
	if(rank==0){ std::cout<<"dt2: "<<ts.dt2()<<std::endl;  }
	if(rank==0){ std::cout<<"alpha: "<<ts.alpha()<<std::endl;  }
		
	MPI_Barrier(MPI_COMM_WORLD);
	wf.forward_psi(); //Loads psi_hat into wf.buf1
	MPI_Barrier(MPI_COMM_WORLD);

	REAL previous_step=ts.tau();
	bool just_printed = false;

	for(int step=0; step<params.nsteps(); step++){
		Timer::startTimer(t_step);	
		if(rank==0){ std::cout<<"TIMESTEP "<<step<<std::endl;  }
		if(rank==0){ std::cout<<"z: "<<ts.zz()<<" a: "<<ts.aa()<<std::endl; }

		previous_step =ts.tau();
		ts.set_step(); //sets tau to the default (afin-ain)/nsteps

		if(ts.m_apr.front() - ts.aa() < ts.tau()){
			ts.set_step(ts.m_apr.front() - ts.aa());
		} else if(just_printed) {//this compensates for the previous understep and puts us where we should be
			ts.set_step(2.0*ts.tau() - previous_step); 
			just_printed=false;
		}

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
		wf.check_phi();
		wf.kick(ts.tau(), ts.aa(), ts.adot(), ts.phiscal());
		ts.advanceHalfStep();
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.forward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		/////////////////////////////////////////////////////	
		

		///////////////////////STREAM/////////////////////////	
	
		wf.stream(ts.tau2(), ts.aa(), ts.adot());
		
		/////////////////////////////////////////////////////	
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(rank==0){ std::cout<<"Ending a: "<<ts.aa()<<std::endl;}

		if(ts.aa()==ts.m_apr.front() || step+1==params.nsteps() || params.debug()){
			MPI_Barrier(MPI_COMM_WORLD);
			wf.backward_psi();
			MPI_Barrier(MPI_COMM_WORLD);
			//wf.computeRho(ts.aa());
			wf.computeRho(1.0);
		
			MPI_Barrier(MPI_COMM_WORLD);
			solver->forward_solve(wf.m_rho.data());
			MPI_Barrier(MPI_COMM_WORLD);

			if(rank==0){ std::cout<<wf.m_rho[0]<<std::endl; }

			std::ostringstream ost;
			ost<<step<<".z"<<ts.zz();
			writePk(solver, params, params.outBase(), ost.str());
			MPI_Barrier(MPI_COMM_WORLD);
			
			//wf.writeWavefunction(params.outBase(), step);
			wf.write1D(params.outBase(), step);

			MPI_Barrier(MPI_COMM_WORLD);
			wf.forward_psi();
			MPI_Barrier(MPI_COMM_WORLD);
			if(ts.aa()==ts.m_apr.front()){
				ts.m_apr.erase(ts.m_apr.begin()); //updates next a to print at
			}
			just_printed=true;
		}
		if(rank==0){ std::cout<<"Step "<<step<<" done!"<<std::endl; }
		Timer::stopTimerStats(t_step);
	}

	Timer::stopTimerStats(t_total);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
