#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <string>
#include <vector>

//SWFFT
#include "complex-type.h"
#include "Distribution.hpp"
#include "Dfft.hpp"

#include "initread.h"

//common
#include "GenericIO.h" 
#include "sctypes.h"

class Wavefunction {
	
	public:
		Wavefunction(const Parameters &, hacc::Distribution *d);
		~Wavefunction();

		void setPsiFromCopy(std::vector<GRID_T> xx, std::vector<GRID_T> yy, 
			   	std::vector<GRID_T> zz,  std::vector<complex_t> original);
		
		void computeRho(double a);
		void kick(double step, double a, double prefactor, double ps);
		void stream(double step, double a, double prefactor);
		void forward_psi();
		void backward_psi();
		void fft_scale();

		void coherent_kick(double timestep, double omega);
		void coherent_stream(double timestep, double omega);

		void check_phi();
		void check_phase(int step);

		void readPhase(std::string inFileBase);
		void writePhase(std::string outFileBase);
		
		void writeWavefunction(std::string outFileBase, int step); 
		void write1D(std::string outFileBase, int step); 
		void readWavefunction(std::string inFileBase);
		void naiveWrite();

		int Ng() { return m_Ng; } 
		float mass() { return m_mass; }
		hacc::Distribution d() { return m_d; }
		
		//primary arrays	
		std::vector<complex_p> m_psi;
		
		std::vector<complex_t> m_rho; //complex_t to faciliate FFTs
		std::vector<complex_t> m_phi;

		std::vector<complex_t> m_buf1;
		std::vector<complex_t> m_buf2;

		std::vector<GRID_T> m_xx;	
		std::vector<GRID_T> m_yy;	
		std::vector<GRID_T> m_zz;	
	private:
		int m_ng;
		long m_Ng;	
		int rank;

		int m_local_k_dim[3];
		int m_self_k_coord[3];

		double m_rL;
		double m_mass;
		float m_hubble;

		double coeff;

		hacc::Distribution &m_d;
		hacc::Dfft *m_dfft;
};

#endif
