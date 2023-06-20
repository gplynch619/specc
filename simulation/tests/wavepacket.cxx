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

	if(rank==0){ std::cout<<"x dim: "<<d.local_ng_3d(0)<<std::endl; }
	if(rank==0){ std::cout<<"y dim: "<<d.local_ng_3d(1)<<std::endl; }
	if(rank==0){ std::cout<<"z dim: "<<d.local_ng_3d(2)<<std::endl; }

	int ng = params.ng();
	double cen = ng/2.; //center
	double tpi = 2.0*4.0*atan(1.0);
	double pi  = 4.0*atan(1.0);

	double cell = params.rL()/(1.0*ng);

	cen= cell*cen;

	double k0 = 0;
	
	double width = params.width();

	std::cout<<"width "<<width<<std::endl;
	std::cout<<"params.width "<<params.width()<<std::endl;

	s.resize(ng);
	long index =0;
	int debug_ind = 0;
	for(long i=0; i<d.local_ng_3d(0); ++i){
		for(long j=0; j<d.local_ng_3d(1); ++j){
			for(long k=0; k<d.local_ng_3d(2); ++k){
				double global_x, global_y, global_z;
				global_x = d.local_ng_3d(0)*d.self_3d(0) + i; 
				global_y = d.local_ng_3d(1)*d.self_3d(1) + j;
				global_z = d.local_ng_3d(2)*d.self_3d(2) + k;

				xx.at(index)= global_x;
				yy.at(index) = global_y;
				zz.at(index) = global_z;
				
				global_x = cell*global_x;
				global_y = cell*global_y;
				global_z = cell*global_z;
	
				//double shiftx = global_x-cen;
				double shiftx = ((global_x<=cen) ? global_x : params.rL()-global_x); 
				shiftx=shiftx-cen;								
				double coeff = pow(pi, -0.25)*pow(1.0, -0.5);
				double mag = coeff*exp(-(shiftx*shiftx)/(2.));
				psi.at(index) = mag;
				if(global_y==cen && global_z==cen){	
					s.at(debug_ind)=mag;
					debug_ind++;
				}
				index++;
			}
		}
	}

	
	std::string fname;
	std::ofstream OutFile;
	int it;
	fname="debug.out";
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(it=0; it<ng; ++it){
		OutFile << s[it].real() <<" "<<s[it].imag()<<std::endl;
	}

	Wavefunction wf(params, &d);
	wf.setPsiFromCopy(xx, yy, zz, psi);

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0) printf("\n");

	MPI_Barrier(MPI_COMM_WORLD);
	
	//we now have a gaussian we would like to shift

	wf.forward_psi();

	MPI_Barrier(MPI_COMM_WORLD);

	double shift = params.shift();
	
	int local_k[3];
	int self_k[3];
	int k[3];

	double nq=ng/2;
/*
	for(int i=0; i<3; i++){
		local_k[i] = d.local_ng_2d_z(i);
		self_k[i] = d.self_2d_z(i);
	}

	for(int local_k0=0; local_k0 < local_k[0]; ++local_k0){
		k[0] = local_k0 + self_k[0]*local_k[0];
		if(k[0]>nq){ k[0] = -MOD(params.ng()-k[0], params.ng());}
		for(int local_k1=0; local_k1 < local_k[1]; ++local_k1){
			k[1] = local_k1 + self_k[1]*local_k[1];
			if(k[1]>=nq){ k[1] = -MOD(params.ng()-k[1], params.ng());}
			for(int local_k2=0; local_k2 < local_k[2]; ++local_k2){
				k[2] = local_k2 + self_k[2]*local_k[2];
				if(k[2]>=nq){ k[2] = -MOD(params.ng()-k[2], params.ng());}
				long id = (local_k0*local_k[1]+local_k1)*local_k[2]+local_k2;
			
				double kk = (tpi/(1.0*params.rL()))*k[0];

				double exponent = -cen*kk;
				
				wf.m_psi.at(id).phase += exponent;
			}
		}
	}

*/
	MPI_Barrier(MPI_COMM_WORLD);
	wf.backward_psi();
	MPI_Barrier(MPI_COMM_WORLD);

	std::cout<<"cen "<<cen<<std::endl;

	wf.writeWavefunction(params.outBase(),-1);
	
	wf.write1D(params.outBase(), -1);

	return 1;
	//Init done//

	MPI_Barrier(MPI_COMM_WORLD);
	wf.forward_psi();
	MPI_Barrier(MPI_COMM_WORLD);

	REAL omega = params.freq();

	//compute step stuff here//
	REAL nsteps = params.nsteps();
	REAL dt = params.duration()/(1.0*nsteps);
	REAL dt2 = dt/2.0;

	///////////////////////////

	if(rank==0){ std::cout<<"T: "<<params.duration()<<" dt: "<<dt<<std::endl; }


	for(int step=0; step<nsteps; step++){
		if(rank==0){ std::cout<<"TIMESTEP "<<step<<std::endl;  }
		
		wf.coherent_stream(dt2, omega);
		
		MPI_Barrier(MPI_COMM_WORLD);	
		wf.backward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		
		wf.coherent_kick(0.0, 0.0); //no potential
		
		MPI_Barrier(MPI_COMM_WORLD);
		wf.forward_psi();
		MPI_Barrier(MPI_COMM_WORLD);
		
		wf.coherent_stream(dt2, omega);

		MPI_Barrier(MPI_COMM_WORLD);
		
		//if(step%50==0){
		if(step%25==0){
			wf.backward_psi();
			MPI_Barrier(MPI_COMM_WORLD);
			wf.write1D(params.outBase(), step);
			MPI_Barrier(MPI_COMM_WORLD);
			wf.forward_psi();
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//wf.writeWavefunction(params.outBase(), step);
		MPI_Barrier(MPI_COMM_WORLD);
		
		std::ostringstream ost;
		ost << step;
		
		if(rank==0){ std::cout<<"Step "<<step<<" done!"<<std::endl; }
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	wf.backward_psi();
	MPI_Barrier(MPI_COMM_WORLD);
	
	wf.write1D(params.outBase(), nsteps);

	return 1;
}
