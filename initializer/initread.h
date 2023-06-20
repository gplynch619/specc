////////////////////////////////////
//
// initread.h
// This is the basic class that reads in the raw parameter file. Based essentially on Basedata.h 
//
/////////////////////////////////////

#ifndef INITREAD_H
#define INITREAD_H

#include <unistd.h>
#include <string>
#include <stdlib.h>
#include <map>
#include <mpi.h> //needed for rank 0 readin
#include <math.h> //needed for pow

#include <iostream>
#include <fstream>
#include <sstream>


#include "sctypes.h"

class Parameters {
	
	public:
		
		Parameters(std::string pName);
		~Parameters() {}

		
		void readNewParams(std::string &fn);
		std::string getParamsString();
		
		REAL omega_cdm() const { return m_omega_cdm; }
		REAL deut() const { return m_deut; }
		REAL omega_nu() const { return m_omega_nu; }
		REAL hubble() const{ return m_hubble; }
		REAL ss8() const { return m_ss8; }
		REAL ns() const { return m_ns; }
		REAL w_de() const { return m_w_de; }
		REAL wa_de() const { return m_wa_de; }
		REAL Tcmb() const { return m_Tcmb; }
		REAL Zdec() const { return m_Zdec; }
		REAL neff_massless() const { return m_neff_massless; }
		REAL neff_massive() const { return m_neff_massive; }
		double cdm_mass() const { return m_cdm_mass; }
		double scatter() const {return m_scatter; }

		REAL zin() const { return m_zin; }
		bool new_ics() const { return m_new_ics; }
		bool debug() const { return m_debug; }
		int trans() const { return m_trans; }

		int iseed() const { return m_iseed; }
		REAL alpha() const { return m_alpha; }
		int ng() const { return m_ng; }
		REAL rL() const { return m_rL; }
		REAL zfin() const { return m_zfin; }
		std::string inTransfer() const { return m_inTransfer; }
		std::string initTrans() const { return m_initTrans; }
		std::string outBase() const { return m_outBase; }

  		float omega_baryon() const { return m_omega_baryon; }
  		float omega_cb() const { return m_omega_cb; }
  		float omega_matter() const { return m_omega_matter; } 
		float omega_radiation() const { return m_omega_radiation; }
		float f_nu_massless() const { return m_f_nu_massless; }
		float f_nu_massive() const { return m_f_nu_massive; }
		
		size_t nzpr() { return zpr_size; }
		
		REAL* zpr() { return m_zpr; }
		REAL zpr(int i) const { return m_zpr[i]; }
		REAL ain() const { return m_ain; }
		REAL afin() const { return m_afin; }
		REAL* apr() { return m_apr; }
		REAL apr(int i) const { return m_apr[i]; }
		REAL pp() const { return m_pp; }
		REAL pfin() const { return m_pfin; }
		REAL ppr() const { return m_ppr; }
		REAL adot() const { return m_adot; }
		REAL prefactor() const { return m_prefactor; }

		REAL freq() const {return m_qho_freq; } 
		REAL width() const { return m_width; }
		REAL shift() const { return m_shift; }
		REAL duration() const {return m_duration; }

		int pknbin1d() const { return m_pknbin1d; }
		int pknmode3d() const { return m_pknmode3d; }
		bool pklogbins() const { return m_pklogbins; }

		int nsteps() const { return m_nsteps; }

		int* kmode() { return m_kmode; }
		int kmode(int i) const { return m_kmode[i]; }

	private:

		std::map<std::string, std::string> m_params; //map that holds parameter:value pairs from read in

		//Cosmological parameters
		REAL m_omega_cdm;
		REAL m_deut;
		REAL m_omega_nu;
		REAL m_hubble;
		REAL m_ss8;
		REAL m_ns;
		REAL m_w_de;
		REAL m_wa_de;
		REAL m_Tcmb;
		REAL m_Zdec;
		REAL m_neff_massless;
		REAL m_neff_massive;
		double m_cdm_mass;
		double m_scatter;

  		float m_omega_baryon;
  		float m_omega_cb;
  		float m_omega_matter;
		float m_omega_radiation;
		float m_f_nu_massless;
		float m_f_nu_massive;

		//Initializer set up
		REAL m_zin;
		REAL m_zfin;
		REAL* m_zpr;
		size_t zpr_size;
		bool m_new_ics;
		int m_trans;	

		//Code Parameters
		int m_iseed;
		REAL m_alpha;
		int m_ng;
		REAL m_rL;
		
		REAL m_ain;
		REAL m_afin;
		REAL* m_apr;
		REAL m_pp;
		REAL m_pfin;
		REAL m_ppr;
		REAL m_adot;
		REAL m_prefactor;

		int m_nsteps;
		bool m_debug;

		std::string m_inTransfer;
		std::string m_initTrans;
		std::string m_outBase;

		//power spec params
		int m_pknbin1d;
		int m_pknmode3d;
		bool m_pklogbins;

		//plane wave params
		int m_kmode[3];

		//Gaussian params
		REAL m_qho_freq;
		REAL m_width;
		REAL m_shift;
		REAL m_duration;
};

inline bool getRank0Stream(const char *filename, std::stringstream &ss){
	MPI_Comm comm;

	int rank;

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	ss.clear();
	if(rank==0){
		std::ifstream fs(filename);
		if(!fs.is_open()){
			fprintf(stderr, "ERROR: failed to open file \"%s\"\n", filename);
			MPI_Abort(comm,1);
			return false;
		}

		ss.str(std::string(""));
		ss << fs.rdbuf();
		fs.close();

		std::string s = ss.str();
		int sz = (int) s.size()+1;
		MPI_Bcast(&sz, 1, MPI_INT, 0, comm);
		MPI_Bcast(const_cast<char*>(s.c_str()), sz, MPI_CHAR, 0, comm);
	} else {
		int sz;
		MPI_Bcast(&sz, 1, MPI_INT, 0, comm);

		char *buffer = new char[sz];
		MPI_Bcast(buffer, sz, MPI_CHAR, 0, comm);

		ss.str(std::string(buffer));
		delete [] buffer;
	}

	return true;
}
#endif
