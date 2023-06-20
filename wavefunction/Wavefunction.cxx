#include "Wavefunction.h"
#include "BasicDefinitions.h" 

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#define DIMENSION 3

inline int MOD(int x, int y) { return (x-y*(x/y)); }

Wavefunction::Wavefunction(const Parameters & params, hacc::Distribution *dist) :
	m_ng(dist->global_ng(0)),
	m_Ng(dist->local_size()),
	rank(dist->self()),
	m_d(*dist)
{
	m_psi.resize(m_Ng);
	m_xx.resize(m_Ng);
	m_yy.resize(m_Ng);
	m_zz.resize(m_Ng);

	m_rho.resize(m_Ng);
	m_phi.resize(m_Ng);
	
	m_buf1.resize(m_Ng);
	m_buf2.resize(m_Ng);

	
	m_dfft = new hacc::Dfft(m_d,
			(complex_t *)&m_buf1[0],
			(complex_t *)&m_buf2[0],
			(complex_t *)&m_buf1[0],
			(complex_t *)&m_buf2[0]);
	
	for(int i=0; i<DIMENSION; i++){
		m_local_k_dim[i] = m_dfft->local_ng_kspace(i);
		m_self_k_coord[i] = m_dfft->self_kspace(i);
	}
	m_rL = params.rL();
	m_hubble = params.hubble();
	m_mass = params.cdm_mass();
	//std::cout<<m_rL<<" "<<NU<<" "<<m_hubble<<" "<<m_mass<<" "<<m_ng<<std::endl;
	//REAL x0 = m_rL/m_ng;
	coeff = NU*m_hubble/(m_mass*(m_rL/(m_ng-1))*(m_rL/(m_ng-1)));
	//coeff = NU*m_hubble/(m_mass*x0*x0);
	//coeff = 1.0; 
	if(rank==0){ std::cout<<"WF COEFF "<<coeff<<std::endl; }
}

Wavefunction::~Wavefunction(){
	delete m_dfft;
};


void Wavefunction::computeRho(double a){
	for(long i=0; i<m_Ng; ++i){
		m_rho[i].real(pow(a, -3.0)*pow(m_psi[i].mag, 2.0));
		m_rho[i].imag(0.0);
	}
}

//might need a unit change before doing this?
//assumes psi in kspace
void Wavefunction::stream(double timestep, double a, double adot){ 
	int k[3]; //kx, ky, kz for a given point
	double tpi = 2.0*atan(1.0)*4.0;	
	double tpiL = tpi/(1.0*m_ng);
	const int nq = m_ng/2;

	double max_k = 0.0;
	double mk[3];
	long index=0;
	for(int local_k0 =0; local_k0 < m_local_k_dim[0]; ++local_k0){ 
		k[0] = local_k0 + m_self_k_coord[0]*m_local_k_dim[0];
		if(k[0]>nq){ k[0] = -MOD(m_ng-k[0], m_ng); }
		for(int local_k1 = 0; local_k1 < m_local_k_dim[1]; ++local_k1){
			k[1] = local_k1 + m_self_k_coord[1]*m_local_k_dim[1];
			if(k[1]>=nq){ k[1] = -MOD(m_ng-k[1], m_ng); }
			for(int local_k2=0; local_k2 < m_local_k_dim[2]; ++local_k2){
				k[2] = local_k2 + m_self_k_coord[2]*m_local_k_dim[2];
				if(k[2]>=nq){ k[2] = -MOD(m_ng-k[2], m_ng); }
				long index = (local_k0*m_local_k_dim[1]+local_k1)*m_local_k_dim[2]+local_k2;
				double kk = tpiL*sqrt((double) (k[0]*k[0]+k[1]*k[1]+k[2]*k[2]));

				REAL stream_term;	
				stream_term = -kk*kk*timestep*coeff/(2.0*a*a*adot);
				
				if((rank==0) && (index==1)){ std::cout<<"Stream coeff: "<<stream_term/(kk*kk)<<std::endl; }
				if((rank==0) && (index==1)){ std::cout<<"Stream term: (1.0,"<<stream_term<<")"<<std::endl; }
				
				m_psi.at(index).phase += stream_term;

				index++;
			}
		}
	}
	 //std::cout<<"RANK "<<rank<<std::endl<<"MAX K CALCULATED: "<<max_k<<" k[0] = "
	//	<<mk[0]<<" k[1] = "<<mk[1]<<" k[2] = "<<mk[2]<<std::endl; 
}

void Wavefunction::kick(double timestep, double a_in, double adot, double ps){

	for(int i = 0; i < m_Ng; ++i){
		if((rank==0) && (i==1)){ std::cout<<"psi[1]: ("<<m_psi.at(1).mag<<","
			<<m_psi.at(1).phase<<")"<<std::endl; }
		REAL kick_term;
		kick_term = -timestep*ps*m_phi[i].real()/(adot*coeff); 
		
		if((rank==0) && (i==1)){ std::cout<<"kick term: (1.0,"<<kick_term<<")"<<std::endl; }
		
		m_psi.at(i).phase += kick_term;
	}
}

void Wavefunction::coherent_kick(double timestep, double omega){

	REAL v0 = m_ng/2.0; //position of the potential minimum (centered)
	REAL cell = m_rL/(1.0*m_ng);
	REAL tpi = 2.0*4.0*atan(1.0); 

	long index=0;
	for(long i=0; i<m_d.local_ng_3d(0); ++i){
		REAL global_x = m_d.local_ng_3d(0)*m_d.self_3d(0) + i;
		for(long j=0; j<m_d.local_ng_3d(1); ++j){
			for(long k=0; k<m_d.local_ng_3d(2); ++k){
				REAL kick;
				kick = -timestep*0.5*tpi*omega*pow(cell, 2.0)*pow(global_x-v0, 2.0);
				m_psi.at(index).phase += kick;
				index++;
			}
		}
	}

}

void Wavefunction::coherent_stream(double timestep, double omega){
	
	int k[3]; //kx, ky, kz for a given point
	double tpi = 2.0*atan(1.0)*4.0;	
	double tpiNg = tpi/(1.0*m_ng);
	const int nq = m_ng/2;
	
	REAL cell = m_rL/(1.0*m_ng);

	long index=0;
	for(int local_k0 =0; local_k0 < m_local_k_dim[0]; ++local_k0){ 
		k[0] = local_k0 + m_self_k_coord[0]*m_local_k_dim[0];
		if(k[0]>nq){ k[0] = -MOD(m_ng-k[0], m_ng); }
		for(int local_k1 = 0; local_k1 < m_local_k_dim[1]; ++local_k1){
			k[1] = local_k1 + m_self_k_coord[1]*m_local_k_dim[1];
			if(k[1]>=nq){ k[1] = -MOD(m_ng-k[1], m_ng); }
			for(int local_k2=0; local_k2 < m_local_k_dim[2]; ++local_k2){
				k[2] = local_k2 + m_self_k_coord[2]*m_local_k_dim[2];
				if(k[2]>=nq){ k[2] = -MOD(m_ng-k[2], m_ng); }
				long index = (local_k0*m_local_k_dim[1]+local_k1)*m_local_k_dim[2]+local_k2;
				double kk = tpiNg*k[0];

				REAL stream_term;	
				stream_term = -0.5*timestep*tpi*omega*pow(1.0/cell, 2.0)*kk*kk;
				
				m_psi.at(index).phase += stream_term;

				index++;
			}
		}
	}

}

void Wavefunction::forward_psi(){
	double pi = 4.0*atan(1.0);
	for(long i=0; i<m_psi.size(); i++){
		/*
		double re, im;
		re = m_psi.at(i).mag*cos(m_psi.at(i).phase);
		im = m_psi.at(i).mag*sin(m_psi.at(i).phase);
		std::complex<double> c(re, im);
		*/
		m_buf1.at(i) = std::polar(m_psi.at(i).mag, m_psi.at(i).phase);
	}
	m_dfft->forward(&m_buf1[0]);
	for(long i=0; i<m_psi.size(); i++){
		m_psi.at(i).mag = std::abs(m_buf1.at(i)); 
		m_psi.at(i).phase = std::arg(m_buf1.at(i)); 
	}
	fft_scale();
}

//apply scalings here if needed
void Wavefunction::backward_psi(){
	double pi = 4.0*atan(1.0);
	for(long i=0; i<m_psi.size(); i++){
		/*
		double re, im;
		re = m_psi.at(i).mag*cos(m_psi.at(i).phase);
		im = m_psi.at(i).mag*sin(m_psi.at(i).phase);
		std::complex<double> c(re, im);
		*/
		m_buf1.at(i) = std::polar(m_psi.at(i).mag, m_psi.at(i).phase);
	}
	m_dfft->backward(&m_buf1[0]);
	for(long i=0; i<m_psi.size(); i++){
		m_psi.at(i).mag = std::abs(m_buf1.at(i));
		m_psi.at(i).phase = std::arg(m_buf1.at(i));
	}
	fft_scale();
}

void Wavefunction::fft_scale(){
	float scal=pow(1.0*m_ng, -1.5);
	for(int i=0; i<m_Ng; ++i){
		m_psi.at(i).mag*=scal;
	}
}

void Wavefunction::naiveWrite(){
	std::string fname;
	std::ofstream OutFile;
	int i;
	fname="xvec2.out";
	OutFile.open(fname, std::ios::out | std::ios::app);
	for(i=0; i<m_xx.size(); ++i){
		OutFile << m_xx[i] <<std::endl;
	}
}

void Wavefunction::write1D(std::string outFileBase, int step){
	
	std::ostringstream ost;
	int ycoord, zcoord;
	ycoord = m_ng/2;
	zcoord = m_ng/2;
	int local_ngx, local_ngy, local_ngz;
	int self_x, self_y, self_z;
	self_x=m_d.self_3d(0);
	self_y=m_d.self_3d(1);
	self_z=m_d.self_3d(2);
	local_ngx=m_d.local_ng_3d(0);
	local_ngy=m_d.local_ng_3d(1);
	local_ngz=m_d.local_ng_3d(2);

	std::vector<GRID_T> local_x;
	std::vector<double> mag;
	std::vector<double> arg;


	int count=0;
	for(int i=0; i<local_ngx; i++){
		int global_i = i+self_x*local_ngx;
		for(int j=0; j<local_ngy; j++){
			int global_j = j+self_y*local_ngy;
			for(int k=0; k<local_ngz; k++){
				int global_k = k+self_z*local_ngz;
				if(global_j==ycoord & global_k==zcoord){
					long index = (i*local_ngy+j)*local_ngz+k;
					local_x.push_back(m_xx.at(index));
					mag.push_back(m_psi.at(index).mag);
					arg.push_back(m_psi.at(index).phase);
					count+=1;
				}	
			}
		}
	}
	
	int color=MPI_UNDEFINED;
	if(count>0){
		color = 0;
	}

	if(step==-1){
		ost<<outFileBase<<".1D.ini.gio";
	} else {
		ost<<outFileBase<<".1D."<<step<<".gio";
	}
	
	//We will now split MPI_COMM_WORLD and initialize a GIO writer using only the ranks that actually have
	//elements

	MPI_Comm subcomm;
	MPI_Comm_split(m_d.parent_comm(), color, rank, &subcomm);

	if(subcomm != MPI_COMM_NULL){
		gio::GenericIO GIO(subcomm, ost.str());
		
		unsigned CoordFlagsX = gio::GenericIO::VarIsPhysCoordX;

		int elem_per_rank = m_ng/m_d.nproc_3d(0);

		GIO.setNumElems(elem_per_rank);
		GIO.setPhysOrigin(0.0);
		local_x.resize(elem_per_rank+GIO.requestedExtraSpace()/sizeof(GRID_T));
		mag.resize(elem_per_rank+GIO.requestedExtraSpace()/sizeof(double));
		arg.resize(elem_per_rank+GIO.requestedExtraSpace()/sizeof(double));

		GIO.addVariable("x", local_x.data(), gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("mag", mag.data(),  gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("arg", arg.data(),  gio::GenericIO::VarHasExtraSpace);

		GIO.write();
		MPI_Barrier(subcomm);
		MPI_Comm_free(&subcomm);
	}
}

void Wavefunction::writeWavefunction(std::string outFileBase, int step){
	std::ostringstream ost;
	
	if(step==-1){
		ost<<outFileBase<<".step.ini.gio";
	} else {
		ost<<outFileBase<<".step."<<step<<".gio";
	}
	
	gio::GenericIO GIO(m_d.parent_comm(), ost.str());
	
	std::vector<double> mag;
	std::vector<double> arg;

	for(long i=0; i<m_Ng; ++i){
		mag.push_back(m_psi.at(i).mag);
		arg.push_back(m_psi.at(i).phase);
	}
	
	unsigned CoordFlagsX = gio::GenericIO::VarIsPhysCoordX;
	unsigned CoordFlagsY = gio::GenericIO::VarIsPhysCoordY;
	unsigned CoordFlagsZ = gio::GenericIO::VarIsPhysCoordZ;
	
	
	GIO.setNumElems(m_Ng);
	GIO.setPhysOrigin(0.0);
	
	m_xx.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	m_yy.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	m_zz.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	mag.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(double));
	arg.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(double));

	GIO.addVariable("x", m_xx.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("y", m_yy.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("z", m_zz.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("mag", mag.data(),  gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("arg", arg.data(),  gio::GenericIO::VarHasExtraSpace);

	GIO.write();

	m_xx.resize(m_Ng);
	m_yy.resize(m_Ng);
	m_zz.resize(m_Ng);
}

void Wavefunction::readWavefunction(std::string inFileBase){
	
	std::ostringstream ost;
	ost<<inFileBase<<".gio";

	gio::GenericIO GIO(m_d.parent_comm(), ost.str());
	GIO.openAndReadHeader();

	m_Ng = GIO.readNumElems();

	std::vector<double> mag;
	std::vector<double> arg;
	
	m_xx.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	m_yy.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	m_zz.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	mag.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(double));
	arg.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(double));


	GIO.addVariable("x", m_xx.data(), gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("y", m_yy.data(), gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("z", m_zz.data(), gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("mag", mag.data(), gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("arg", arg.data(), gio::GenericIO::VarHasExtraSpace);

	GIO.readData();

	for(long i=0; i<m_Ng; ++i){
		m_psi.at(i).mag = mag[i];
		m_psi.at(i).phase = arg[i];
	}

}

void Wavefunction::readPhase(std::string inFileBase){
	
	std::ostringstream ost;
	ost<<inFileBase<<".gio";

	gio::GenericIO GIO(m_d.parent_comm(), ost.str());
	GIO.openAndReadHeader();

	m_Ng = GIO.readNumElems();

	std::vector<double> arg;
	arg.resize(m_Ng);

	GIO.addVariable("arg", arg.data(), gio::GenericIO::VarHasExtraSpace);

	GIO.readData();

	for(long i=0; i<m_Ng; ++i){
		m_psi.at(i).phase = arg[i];
	}

}

void Wavefunction::writePhase(std::string outFileBase){
	std::ostringstream ost;
	
	ost<<outFileBase<<".gio";
	
	gio::GenericIO GIO(m_d.parent_comm(), ost.str());
	
	std::vector<double> arg;

	for(long i=0; i<m_Ng; ++i){
		arg.push_back(m_psi.at(i).phase);
	}
	
	GIO.setNumElems(m_Ng);
	GIO.setPhysOrigin(0.0);
	
	arg.resize(m_Ng+GIO.requestedExtraSpace()/sizeof(double));

	GIO.addVariable("arg", arg.data(),  gio::GenericIO::VarHasExtraSpace);

	GIO.write();
}

void Wavefunction::check_phase(int step){
	double min_psi, max_psi, ave_psi;

	min_psi=1.0e30;
	max_psi=-1.0e30;
	ave_psi=0.0;

	for(long i=0; i<m_Ng; ++i){
		ave_psi += m_psi.at(i).phase;
		if(m_psi.at(i).phase > max_psi){ max_psi = m_psi.at(i).phase; }
		if(m_psi.at(i).phase < min_psi){ min_psi = m_psi.at(i).phase; }
	}
	ave_psi/=m_Ng;

	MPI_Allreduce(MPI_IN_PLACE, &max_psi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &min_psi, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &ave_psi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(rank==0){
		std::cout<<std::endl << "Max and min value of potential in code units: "
			<< max_psi <<" "<<min_psi<<std::endl;
		ave_psi/=m_d.nproc();
		std::cout<< "Average value of phase in code units: " << ave_psi <<std::endl;
	}
	/*	
	if(rank==0){

		std::string outName = "phase_values.txt";
		
		FILE *outFile = fopen(outName.c_str(),"a");

		fprintf(outFile, "%i\t%f\t%f\t%f\n", 
			step, 
			ave_psi,	
			min_psi,
			max_psi);
	  
		fclose(outFile);
	}
	*/
	MPI_Barrier(MPI_COMM_WORLD);
}

void Wavefunction::check_phi(){
	double min_phi, max_phi, ave_phi, max_im, min_im;

	min_phi=1.0e30;
	max_phi=-1.0e30;
	min_im=1.0e30;
	max_im=-1.0e30;
	ave_phi=0.0;

	for(long i=0; i<m_Ng; ++i){
		ave_phi += m_phi.at(i).real();
		if(m_phi.at(i).real() > max_phi){ max_phi = m_phi.at(i).real(); }
		if(m_phi.at(i).imag() > max_im){ max_im = m_phi.at(i).imag(); }
		if(m_phi.at(i).real() < min_phi){ min_phi = m_phi.at(i).real(); }
		if(m_phi.at(i).imag() < min_im){ min_im = m_phi.at(i).imag(); }
	}
	ave_phi/=m_Ng;

	MPI_Allreduce(MPI_IN_PLACE, &max_phi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &max_im, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &min_phi, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &min_im, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &ave_phi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(rank==0){
		std::cout<<std::endl << "Max and min value of potential in code units: "
			<< max_phi <<" "<<min_phi<<std::endl;
		std::cout<<"Max and min imaginary part of potential in code units: "
			<< max_im <<" "<<min_im<<std::endl;
		ave_phi/=m_d.nproc();
		std::cout<< "Average value of phi in code units: " << ave_phi <<std::endl;
	}
}


void Wavefunction::setPsiFromCopy(std::vector<GRID_T> xx, std::vector<GRID_T> yy,
		std::vector<GRID_T> zz, std::vector<complex_t> original){ 
	
	m_xx.swap(xx);
	m_yy.swap(yy);
	m_zz.swap(zz);
	double pi = 4.0*atan(1.0);
	for(long i=0; i<original.size(); i++){
		m_psi.at(i).mag = std::abs(original.at(i));		
		//m_psi.at(i).phase = std::arg(original.at(i)) + pi;		
		m_psi.at(i).phase = std::arg(original.at(i));		
	}
}
