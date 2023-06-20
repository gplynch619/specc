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

//	Timer::startTimer(t_total);
//	Timer::startTimer(t_init);

	std::string paramout = params.getParamsString();

	if(rank == 0){ std::cout << paramout << std::endl; }
		
	hacc::Distribution d(MPI_COMM_WORLD, params.ng());

	std::vector<GRID_T> xx;
	std::vector<GRID_T> yy;
	std::vector<GRID_T> zz;
	std::vector<complex_t> psi;
	std::vector<complex_t> s;
	long My_Ng = d.local_ng_3d(0)*d.local_ng_3d(1)*d.local_ng_3d(2); 

	std::cout<<My_Ng<<std::endl;
	xx.resize(My_Ng);
	yy.resize(My_Ng);
	zz.resize(My_Ng);
	psi.resize(My_Ng);
	s.resize(My_Ng);

	if(rank==0){ std::cout<<"x dim: "<<d.local_ng_3d(0)<<std::endl; }
	if(rank==0){ std::cout<<"y dim: "<<d.local_ng_3d(1)<<std::endl; }
	if(rank==0){ std::cout<<"z dim: "<<d.local_ng_3d(2)<<std::endl; }

	double cen = params.ng()/2;
	
	long index =0;
	double k0 = 0.5*(8.0*atan(1.0))*params.ng()/params.rL();
	for(long i=0; i<d.local_ng_3d(0); ++i){
		for(long j=0; j<d.local_ng_3d(1); ++j){
			for(long k=0; k<d.local_ng_3d(2); ++k){
				long global_x, global_y, global_z;
				global_x = d.local_ng_3d(0)*d.self_3d(0) + i; 
				global_y = d.local_ng_3d(1)*d.self_3d(1) + j;
				global_z = d.local_ng_3d(2)*d.self_3d(2) + k;

				double rdotr = (global_x-cen)*(global_x - cen) + (global_y-cen)*(global_y-cen) 
					+ (global_z-cen)*(global_z-cen);
				xx.at(index)= global_x;
				yy.at(index) = global_y;
				zz.at(index) = global_z;
#ifdef SODGAUSS
				double mag = exp(-rdotr/(2*params.scatter()*params.scatter()));
				complex_t phase;
				phase.imag(k0*sqrt(rdotr));
				//std::cout<<"("<<global_x<<","<<global_y<<","<<global_z<<") "<<wp<<std::endl;
				psi.at(index) = mag*std::exp(phase);
#else
				double mag;
				if(rdotr<params.scatter()*params.scatter()){
					mag=1.0;
				} else {
					mag=0;
				}
				psi.at(index) = std::polar(mag, 0.0);
#endif
				index++;
			}	
		}
	}

	Wavefunction wf(params, &d);
	wf.setPsiFromCopy(xx, yy, zz, psi);

	MPI_Barrier(MPI_COMM_WORLD);

//	Timer::stopTimerStats(t_init);
	if(rank==0) printf("\n");

	SolverDiscrete *solver = new SolverDiscrete(d);

	MPI_Barrier(MPI_COMM_WORLD);

	solver->solve(wf.m_rho.data(), wf.m_phi.data());

	MPI_Barrier(MPI_COMM_WORLD);
	
	int nsteps = params.nsteps();
	float delt=params.zfin();
	
	wf.writeWavefunction(params.outBase(),-1);

	for(int step=0; step<params.nsteps(); step++){
		if(rank==0){ std::cout<<"TIMESTEP "<<step<<std::endl;  }
		
		wf.kick(delt, 0, 1.0, 1.0);
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.forward_psi();
		MPI_Barrier(MPI_COMM_WORLD);

		wf.stream(2.0*delt, 1.0, 1.0);

		MPI_Barrier(MPI_COMM_WORLD);
		wf.backward_psi();
		MPI_Barrier(MPI_COMM_WORLD);

		wf.kick(delt, 0, 1.0, 1.0);

		wf.computeRho(1.0);
		MPI_Barrier(MPI_COMM_WORLD);
		solver->solve(wf.m_rho.data(), wf.m_phi.data());	
		MPI_Barrier(MPI_COMM_WORLD);

		wf.writeWavefunction(params.outBase(), step);
		//wf.writeWavefunction(params.outBase(), step);
		
		std::ostringstream ost;
		ost << step;
		
		if(rank==0){ std::cout<<"Step "<<step<<" done!"<<std::endl; }
	}


	return 1;
}
