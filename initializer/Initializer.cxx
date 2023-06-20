#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <exception> 

//Include RNG
#include "CBRNG_Random.h"

#include "sctypes.h"
#include "BasicDefinitions.h"
#include "CosmoClass.h"

//Include Wavefunction
#include "Wavefunction.h"
//Include SWFFT headers

#include "complex-type.h"
#include "Distribution.hpp"
#include "Dfft.hpp"

inline int MOD(int x, int y) { return (x - y*(x/y)); }

class Initializer {
	
	public:
		void initPsi(int32_t NG_Max, CosmoClass &CC, TransferClass &TC, double rL, REAL z_in, 
				unsigned long seed);

		Initializer(hacc::Distribution &dist);
		~Initializer();

		void initPowerSpec(TransferClass &TC, CosmoClass &CC, double rL, REAL z_in);
		void initDensity(unsigned long seed);
		
		void initVelocity(CosmoClass &CC, TransferClass &TC, unsigned long seed, double rL, REAL z_in);
		
		void initPlaneWave(int *mode, double rL, unsigned long seed);
		void handoff(Wavefunction &wf);
		void writePk(std::vector<REAL> &tf, std::vector<REAL> &pk, std::vector<REAL> &kks, 
				TransferClass &TC, double rL);
	
		struct NegativeDensException : public std::exception
		{
			const char * what () const throw ()
    		{
    		return "Negative density in initializer : delta>1";
    		}
		};
		
	protected:
		hacc::Distribution &d;

		int ng;
		size_t My_Ng;

		MPI_Comm comm;
		int MasterProc;
		int MyProc;
		int NumProcs;

		hacc::Dfft *dfft;

		//k-space variables
		int nx1, ny1, nz1;
		int ngx, ngy, ngz;

		//r-space variables
		int tnx1, tny1, tnz1;
		int tngx, tngy, tngz;

		std::vector<complex_t> psit;

		std::vector<double> vel;
		std::vector<double> density;

		std::vector<complex_t> buf1;
		std::vector<complex_t> buf2;
		
		std::vector<complex_t> scratch;
		
		std::vector<double> Pk;

		std::vector<GRID_T> xx;
		std::vector<GRID_T> yy;
		std::vector<GRID_T> zz;

		std::ostringstream m_ost;
};

Initializer::Initializer(hacc::Distribution &dist) :
	d(dist),
	ng(d.global_ng(0)),
	My_Ng(d.local_size()),
	comm(d.parent_comm()),
	MasterProc(0),
	MyProc(d.self()),
	NumProcs(d.nproc()),
	dfft(NULL),
	nx1(-1),
	ny1(-1),
	nz1(-1),
	ngx(-1),
	ngy(-1),
	ngz(-1),
	tnx1(-1),
	tny1(-1),
	tnz1(-1),
	tngx(-1),
	tngy(-1),
	tngz(-1)
{
	Pk.resize(My_Ng);
	vel.resize(My_Ng);
	density.resize(My_Ng);
	scratch.resize(My_Ng);
	psit.resize(My_Ng);
	buf1.resize(My_Ng);
	buf2.resize(My_Ng);
	
	dfft = new hacc::Dfft(d,
			(complex_t *)&buf1[0], // forward_output in -> fs -> fo
			(complex_t *)&buf2[0], // forward_psit
			(complex_t *)&buf1[0], // backward_input
			(complex_t *)&buf2[0]);// backward_psit
	//dfft->forward(in) goes in->fs->...->fo
	//dfft->backward(out) goes bin->bs...->out
	//fo = buf1, bin = buf2

	ngx = dfft->local_ng_kspace(0);
	ngy = dfft->local_ng_kspace(1);
	ngz = dfft->local_ng_kspace(2);
	nx1 = ngx*dfft->self_kspace(0);
	ny1 = ngy*dfft->self_kspace(1);
	nz1 = ngz*dfft->self_kspace(2);
	
	tngx = dfft->local_ng_rspace(0);
	tngy = dfft->local_ng_rspace(1);
	tngz = dfft->local_ng_rspace(2);
	tnx1 = tngx*dfft->self_rspace(0);
	tny1 = tngy*dfft->self_rspace(1);
	tnz1 = tngz*dfft->self_rspace(2);
}

Initializer::~Initializer() {
	delete dfft;
}

void Initializer::initPowerSpec(TransferClass &TC, CosmoClass &CC, double rL, REAL z_in) {
	const double tpi=2.0*4.0*atan(1.0);
	const double tpiL = tpi/rL;
	const int nq = ng/2; //nyquist

	double file_kmax = 0.5*tpi*ng/rL;

	if(MyProc==MasterProc){
		std::cout<<"Generating Pk"<<std::endl;
	}

	if(TC.kmax() < file_kmax){
		std::ostringstream ost;
		ost << "Initializer::initPowerSpec CMBFast file only goes to k = " << TC.kmax() << std::endl;
		throw std::runtime_error(ost.str());
	}

	double Dplus, Ddot;
	CC.GrowthFactorLCDM(z_in, &Dplus, &Ddot);
	
	for(long i=0; i<ngx; ++i){
		long k_i = i+nx1;
		if(k_i>nq){k_i = -MOD(ng-k_i, ng);}
		for(long j=0; j<ngy; ++j){
			long k_j = j+ny1;
			if(k_j >= nq){k_j = -MOD(ng-k_j, ng);}
			for(long k=0; k<ngz; ++k){
				long k_k = k+nz1;
				if(k_k >= nq){k_k = -MOD(ng-k_k, ng);}
				long index = (i*ngy+j)*ngz+k;
				double kk = tpiL*sqrt((double) (k_i*k_i+k_j*k_j+k_k*k_k));
				double tmp = TC.Pk_cb(kk);
				if(TC.TF()==2){
					Pk[index] = tmp;
				} else if (TC.TF()==0) {
					Pk[index] = Dplus*Dplus*tmp;
				}
			}
		}
	}

	if(MyProc==MasterProc){
		std::cout<<"Ending Pk"<<std::endl;
	}

	return;
}

void Initializer::writePk(std::vector<REAL> &tf, std::vector<REAL> &pk, std::vector<REAL> &kks, 
		TransferClass &TC, double rL) {
	const double tpi=2.0*4.0*atan(1.0);
	const double tpiL = tpi/rL;
	const int nq = ng/2; //nyquist

	double file_kmax = 0.5*tpi*ng/rL;

	if(TC.kmax() < file_kmax){
		std::ostringstream ost;
		ost << "Initializer::initPowerSpec CMBFast file only goes to k = " << TC.kmax() << std::endl;
		throw std::runtime_error(ost.str());
	}

	for(long i=0; i<ngx; ++i){
		long k_i = i+nx1;
		if(k_i>nq){k_i = -MOD(ng-k_i, ng);}
		for(long j=0; j<ngy; ++j){
			long k_j = j+ny1;
			if(k_j >= nq){k_j = -MOD(ng-k_j, ng);}
			for(long k=0; k<ngz; ++k){
				long k_k = k+nz1;
				if(k_k >= nq){k_k = -MOD(ng-k_k, ng);}
				long index = (i*ngy+j)*ngz+k;
				double kk = tpiL*sqrt((double) (k_i*k_i+k_j*k_j+k_k*k_k));
				tf.at(index) = TC.Tf_cb(kk);
				pk.at(index) = Pk[index];
				//pk.at(index) = TC.Pk_cb(kk);
				kks.at(index) = kk;
			}
		}
	}

	return;
}

void Initializer::initDensity(unsigned long seed){
	
	if(MyProc == MasterProc){
		std::cout<< "Generating density field" << std::endl;
	}
	
	std::vector<unsigned long> slab_keys;
	slab_keys.resize(tngx);
	GetSlabKeys(&slab_keys[0], tnx1, tngx, seed);

	for(long i=0; i<tngx; ++i) {
		unsigned long ikey = slab_keys[i];
		for(long j=0; j<tngy; ++j){
			long jg = tny1+j;
			for(long k=0; k<tngz; ++k){
				long kg = tnz1+k;
				unsigned long index=jg*ng + kg;
				double rn1, rn2;
				GetRandomDoublesWhiteNoise(rn1, rn2, ikey, index);
				index = (i*tngy+j)*tngz+k;
				double rn = rn1*rn1 + rn2*rn2;
				buf1[index].real(rn2 * sqrt(-2.0*log(rn)/rn));
				buf1[index].imag(0.0);
			}
		}
	}

	MPI_Barrier(comm);
	dfft->forward((complex_t *)&buf1[0]);
	MPI_Barrier(comm);

	if(MyProc==MasterProc){
		std::cout<< "FFT to k-space done!" << std::endl;
	}

	double scal = pow(1.0*ng, -1.5);

	for(long i=0; i<My_Ng; ++i){
		complex_t tmp = buf1[i]; 
		buf1[i].real(sqrt(Pk[i])*scal*tmp.real());
		buf1[i].imag(sqrt(Pk[i])*scal*tmp.imag());
	}

	if(nx1 == 0 && ny1 == 0 && nz1 ==0){
		buf1[0].real(0.0);
		buf1[0].imag(0.0);
	}

	//Back to r-space
	MPI_Barrier(comm);
	dfft->backward((complex_t *)&scratch[0]);
	MPI_Barrier(comm);
	
	if(MyProc==MasterProc){
		std::cout<< "FFT to r-space done!" << std::endl;
	}
	
} 
//At end of the above: 
// delta at z=z_in cached in r-space in psit
// delta at z=0 is caches in r-space in rho

void Initializer::initVelocity(CosmoClass &CC, TransferClass &TC, unsigned long seed, double rL, REAL z_in){
	const REAL tpi=2.0*4.0*atan(1.0);
	const REAL tpiL = tpi/rL;
	const REAL tping = tpi/ng;
	const int nq = ng/2; //nyquist

	double file_kmax = 0.5*tpi*ng/rL;

	double ain=1.0/(1.0+z_in);

	if(MyProc==MasterProc){
		std::cout<<"Generating Pk"<<std::endl;
	}

	if(TC.kmax() < file_kmax){
		std::ostringstream ost;
		ost << "Initializer::initPowerSpec CMBFast file only goes to k = " << TC.kmax() << std::endl;
		throw std::runtime_error(ost.str());
	}
/*	
	for(long i=0; i<ngx; ++i){
		long k_i = i+nx1;
		if(k_i>nq){k_i = -MOD(ng-k_i, ng);}
		for(long j=0; j<ngy; ++j){
			long k_j = j+ny1;
			if(k_j >= nq){k_j = -MOD(ng-k_j, ng);}
			for(long k=0; k<ngz; ++k){
				long k_k = k+nz1;
				if(k_k >= nq){k_k = -MOD(ng-k_k, ng);}
				long index = (i*ngy+j)*ngz+k;
				double kk = tpiL*sqrt((double) (k_i*k_i+k_j*k_j+k_k*k_k));
				double tmp = TC.Pk_total(kk);
				Pk[index] = tmp;
			}
		}
	}

	if(MyProc==MasterProc){
		std::cout<<"Ending Pk"<< std::endl;
	}
	
	std::vector<unsigned long> slab_keys;
	slab_keys.resize(tngx);
	GetSlabKeys(&slab_keys[0], tnx1, tngx, seed);

	for(long i=0; i<tngx; ++i) {
		unsigned long ikey = slab_keys[i];
		for(long j=0; j<tngy; ++j){
			long jg = tny1+j;
			for(long k=0; k<tngz; ++k){
				long kg = tnz1+k;
				unsigned long index=jg*ng + kg;
				double rn1, rn2;
				GetRandomDoublesWhiteNoise(rn1, rn2, ikey, index);
				index = (i*tngy+j)*tngz+k;
				double rn = rn1*rn1 + rn2*rn2;
				buf1[index].real(rn2 * sqrt(-2.0*log(rn)/rn));
				buf1[index].imag(0.0);
			}
		}
	}

	MPI_Barrier(comm);
	dfft->forward((complex_t *)&buf1[0]);
	MPI_Barrier(comm);

	if(MyProc==MasterProc){
		std::cout<< "FFT to k-space done!" << std::endl;
	}

	double scal = pow(1.0*ng, -1.5);
	
	for(long i=0; i<My_Ng; ++i){
		complex_t tmp = buf1[i]; 
		buf1[i].real(sqrt(Pk[i])*scal*tmp.real());
		buf1[i].imag(sqrt(Pk[i])*scal*tmp.imag());
	}

	if(nx1 == 0 && ny1 == 0 && nz1 ==0){
		buf1[0].real(0.0);
		buf1[0].imag(0.0);
	}

	//Back to r-space
	MPI_Barrier(comm);
	dfft->backward((complex_t *)&buf1[0]);
	MPI_Barrier(comm);
	
   	scal = pow(1.0*rL, -1.5);
	for(long i=0; i<My_Ng; ++i){
		buf1[i] *= scal;
	}

	scal = pow(1.0*ng, -1.5);
	
	MPI_Barrier(comm);
	dfft->forward((complex_t *)&buf1[0]);
	MPI_Barrier(comm);
	
	for(long i=0; i<My_Ng; ++i){
		buf1[i] *= scal;
	}
*/

	// density currently holds real double zin density
	
	double Dplus, Ddot, Lk;
	//returns Ddot/H0 (dimensionless)	
	CC.GrowthFactorLCDM(z_in, &Dplus, &Ddot);
	
	double H0 = CC.getCosmoClassParams().h*100.0;

	buf1.resize(density.size());

	for(int i=0; i<buf1.size(); i++){
		double tmp = density.at(i)-1;
		//tmp = tmp/Dplus; //only works for scale invariant Dplus
		buf1.at(i).real(tmp);
		buf1.at(i).imag(0.0);
	}

	double scal = pow(1.0*ng, -1.5);
	double alpha = NU*CC.getCosmoClassParams().h/CC.getCosmoClassParams().cdm_mass;

	MPI_Barrier(comm);
	dfft->forward((complex_t *)&buf1[0]);
	MPI_Barrier(comm);
	
	for(int i=0; i<buf1.size(); i++){
		buf1.at(i) *= scal;
	}

	for(long i=0; i<ngx; ++i){
		long k_i = i+nx1;
		if(k_i>nq){k_i = -MOD(ng-k_i, ng);}
		for(long j=0; j<ngy; ++j){
			long k_j = j+ny1;
			if(k_j >= nq){k_j = -MOD(ng-k_j, ng);}
			for(long k=0; k<ngz; ++k){
				long k_k = k+nz1;
				if(k_k >= nq){k_k = -MOD(ng-k_k, ng);}
				long index = (i*ngy+j)*ngz+k;
				double kk = tpiL*sqrt((double) (k_i*k_i+k_j*k_j+k_k*k_k));
				if(kk==0){
					kk=0.0000001;
				}
				complex_t tmp = buf1[index];
				//CC.MeanGrowth(z_in, kk, &Lk);
				//buf1[index] = ain*Lk*Ddot*tmp/(kk*kk);
				buf1[index] = ain*Ddot*tmp/(alpha*kk*kk);
			}
		}
	}

	if(nx1 == 0 && ny1 == 0 && nz1 ==0){
		buf1[0].real(0.0);
		buf1[0].imag(0.0);
	}

	MPI_Barrier(comm);
	dfft->backward((complex_t *)&scratch[0]);
	MPI_Barrier(comm);

}

void Initializer::initPsi(int32_t NG_Max, CosmoClass &CC, TransferClass &TC, double rL, 
		REAL z_in, unsigned long seed){
	
	if(NG_Max < My_Ng){
		std::ostringstream ost;
		ost << "Initializer::initPsi NG_Max="<<NG_Max<<" My_Ng="<<My_Ng<<std::endl;
		throw std::runtime_error(ost.str());
	}

	if(MyProc==MasterProc){
		std::cout<<"Initializing Power Spectrum..."<<std::endl;
	}
	
	initPowerSpec(TC, CC, rL, z_in);
	
	if(MyProc==MasterProc){
		std::cout<<"Initializing density field..."<<std::endl;
	}
	
	initDensity(seed);
	//"should" be rl	
	double scal = pow(1.0*rL, -1.5);
	for(long i=0; i<My_Ng; ++i){
		scratch[i] *= scal;
		density[i] = scratch[i].real() + 1.0;
	}
/*
	std::string fname;
	std::ofstream OutFile;
	fname="sorc.dens.debug";
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(int i=0; i<ng; ++i){
		OutFile << psit[i].real() <<" "<<psit[i].imag()<<std::endl;
	}
	OutFile.close();
*/
	//scales after last FFT in initDensity
	//psit holds delta at z=z_in, FFT scaled and in r-space
	if(MyProc==MasterProc){
		std::cout<<"Initializing velocity potential..."<<std::endl;
	}
	
	initVelocity(CC, TC, seed, rL, z_in);

	scal = pow(1.0*ng, -1.5);

	for(long i=0; i<My_Ng; i++){
		scratch[i] *= scal;
		vel[i] = scratch[i].real();	
	}

	/////////ZERO VEL START////////
	/*
	for(long i=0; i<My_Ng; i++){
		vel[i] = 0.0;
	}*/
	///////////////////////////////

	//REAL pcoeff = NU*CC.getCosmoClassParams().h/CC.getCosmoClassParams().cdm_mass;
	//REAL pcoeff = .095;
	REAL pcoeff = 1.0;

	std::cout<<"PCOEFFFF "<<pcoeff<<std::endl;

	for(long i=0; i<My_Ng;++i){
		if(density.at(i) < 0){
			throw NegativeDensException();		
		}
		double re, im;
		re = sqrt(density.at(i))*cos(vel.at(i));
		im = sqrt(density.at(i))*sin(vel.at(i));
		std::complex<double> c(re, im);
		psit[i] = c;
	}

	//set x,y,z arrays	
	for(long i=0; i<dfft->local_ng_rspace(0); ++i) {	
		for(long j=0; j<dfft->local_ng_rspace(1); ++j){
			for(long k=0; k<dfft->local_ng_rspace(2); ++k){
				long global_x, global_y, global_z;
				global_x = dfft->local_ng_rspace(0)*dfft->self_rspace(0) + i;	
				global_y = dfft->local_ng_rspace(1)*dfft->self_rspace(1) + j;	
				global_z = dfft->local_ng_rspace(2)*dfft->self_rspace(2) + k;	
				xx.push_back(global_x);
				yy.push_back(global_y);
				zz.push_back(global_z);
			}
		}
	}
}

void Initializer::initPlaneWave(int *mode, double rL, unsigned long seed){
	
	if(MyProc==MasterProc){
		std::cout<<"Initializing plane wave..."<<std::endl;
	}

	const double tpi=2.0*4.0*atan(1.0);
	const double tpiL = tpi/rL;
	const int nq = ng/2; //nyquist

	double file_kmax = 0.5*tpi*ng/rL;

	double phys_mode[3];
	for(int i=0; i<3; i++){
		phys_mode[i]=mode[i];
	}

	for(long i=0; i<My_Ng; i++){
		buf1[i] = 0.0;
	}

	for(long i=0; i<ngx; ++i){
		long k_i = i+nx1;
		if(k_i>nq){k_i = -MOD(ng-k_i, ng);}
		for(long j=0; j<ngy; ++j){
			long k_j = j+ny1;
			if(k_j >= nq){k_j = -MOD(ng-k_j, ng);}
			for(long k=0; k<ngz; ++k){
				long k_k = k+nz1;
				if(k_k >= nq){k_k = -MOD(ng-k_k, ng);}
				long index = (i*ngy+j)*ngz+k;
				double kk = tpiL*sqrt((double) (k_i*k_i+k_j*k_j+k_k*k_k));
				if(fabs(k_i)==phys_mode[0] && fabs(k_j)==phys_mode[1] && fabs(k_k)==phys_mode[2]){
					//fabs takes care of conjugate symmetry so that the iFFT is real
					buf1[index].real(1.0);
					std::cout<<"Setting mode"<<std::endl;
					std::cout<<"("<<k_i<<","<<k_j<<","<<k_k<<")"<<std::endl;
				}
			}
		}
	}

	if(nx1 == 0 && ny1 == 0 && nz1 ==0){
		buf1[0].real(0.0);
		buf1[0].imag(0.0);
	}
	
	MPI_Barrier(comm);
	dfft->backward((complex_t *)&scratch[0]);
	MPI_Barrier(comm);

	if(MyProc==MasterProc){
		std::cout<< "FFT to k-space done!" << std::endl;
	}

	double scal = pow(1.0*rL, -1.5);

	double max_im = -1e30;
	double min_im = 1e30;
	for(long i=0; i<My_Ng; ++i){
		scratch[i] *= scal;
		if(scratch[i].imag() > max_im){ max_im = scratch[i].imag(); }	
		if(scratch[i].imag() < min_im){ min_im = scratch[i].imag(); }	
		density[i] = scratch[i].real() + 1.0;
	}

	//std::cout<<"Rank "<<MyProc<<" has (max, min) of imag: ("<<max_im<<","<<min_im<<")"<<std::endl;

	for(long i=0; i<My_Ng;++i){
		if(density.at(i) < 0){
			throw NegativeDensException();		
		}
		double re, im;
		re = sqrt(density.at(i))*cos(0.0);
		im = sqrt(density.at(i))*sin(0.0);
		std::complex<double> c(re, im);
		psit[i] = c;
	}

	//set x,y,z arrays	
	for(long i=0; i<dfft->local_ng_rspace(0); ++i) {	
		for(long j=0; j<dfft->local_ng_rspace(1); ++j){
			for(long k=0; k<dfft->local_ng_rspace(2); ++k){
				long global_x, global_y, global_z;
				global_x = dfft->local_ng_rspace(0)*dfft->self_rspace(0) + i;	
				global_y = dfft->local_ng_rspace(1)*dfft->self_rspace(1) + j;	
				global_z = dfft->local_ng_rspace(2)*dfft->self_rspace(2) + k;	
				xx.push_back(global_x);
				yy.push_back(global_y);
				zz.push_back(global_z);
			}
		}
	}
}

void Initializer::handoff(Wavefunction &wf){
	wf.setPsiFromCopy(xx,yy,zz,psit);
}

void initialize_wavefunction(int32_t NG_Max, CosmoClass &CC, TransferClass &TC, double rL, 
		double z_in, unsigned long seed, hacc::Distribution &d, Wavefunction &wf){
	
	Initializer init(d);
	init.initPsi(NG_Max, CC, TC, rL, z_in, seed);

	init.handoff(wf);

}

void initialize_planewave(int* mode, double rL, unsigned long seed, 
		hacc::Distribution &d, Wavefunction &wf){
	
	Initializer init(d);
	init.initPlaneWave(mode, rL, seed);

	init.handoff(wf);
}

void initialize_powerspec(int32_t NG_Max, CosmoClass &CC, TransferClass &TC, double rL, 
		double z_in, unsigned long seed, hacc::Distribution &d, Wavefunction &wf){

	std::vector<REAL> ltf;
	std::vector<REAL> lpk;
	std::vector<REAL> lkks;


	long Ng,lNg;
	
	Ng = d.global_size();
	lNg = d.local_size();
	
	if(d.self() == 0){
		std::cout<<"global ng "<<d.global_ng(0)<<std::endl;
		std::cout<<"local Ng "<<lNg<<std::endl;
		std::cout<<"Ng "<<Ng<<std::endl;
	}
	
	ltf.resize(lNg);
	lpk.resize(lNg);
	lkks.resize(lNg);
	
	Initializer init(d);
	try{
	init.initPsi(NG_Max, CC, TC, rL, z_in, seed);
	}
	catch (Initializer::NegativeDensException& e)
	{
		std::cout<<e.what()<<std::endl;
		std::cout<<"This usually means you have to start at a higher redshift."<<std::endl;
	}

	if(d.self() == 0){
		std::cout<<"done with init"<<std::endl;
	}
	init.writePk(ltf, lpk, lkks, TC, rL);

	MPI_Barrier(d.parent_comm());

#ifdef TESTING
	std::vector<REAL> gtf;
	std::vector<REAL> gpk;
	std::vector<REAL> gkks;
	
	if(d.self() == 0){
		gtf.resize(Ng);
		gpk.resize(Ng);
		gkks.resize(Ng);
	}

	if(d.self() == 0){
		std::cout<<"Before gather "<<std::endl;
	}
	MPI_Gather(lkks.data(), lNg, MPI_DOUBLE, gkks.data(), lNg, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(d.self() == 0){
		std::cout<<"k gathered "<<std::endl;
	}
	MPI_Gather(ltf.data(), lNg, MPI_DOUBLE, gtf.data(), lNg, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(d.self() == 0){
		std::cout<<"tf gathered "<<std::endl;
	}
	
	MPI_Gather(lpk.data(), lNg, MPI_DOUBLE, gpk.data(), lNg, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(d.self() == 0){
		std::cout<<"pk gathered "<<std::endl;
	}
	
	if(d.self() == 0){
		std::string fname;
		std::ofstream OutFile;
		int i;
		fname="debug.pk.out";
		OutFile.open(fname, std::ios::out | std::ios::app);
		for(i=0; i<Ng; ++i){
			OutFile << gkks.at(i) <<" "<< gtf.at(i)<<" "<<gpk.at(i)<<std::endl;
		}
	}
#endif 
	MPI_Barrier(MPI_COMM_WORLD);

	init.handoff(wf);

}
