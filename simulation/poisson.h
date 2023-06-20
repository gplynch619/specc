#ifndef POISSON_H
#define POISSON_H

#include <algorithm>
#include <vector>
#include <cassert>

#include "complex-type.h"
#include "Distribution.hpp"
#include "power.h"
#include "Dfft.hpp"

//////////////////////////////////////////////////
//
//Base class for Poisson solver. Specific implementations
//must be defined as derived classes.
//
//////////////////////////////////////////////////

inline int MOD(int x, int y) { return (x-y*(x/y)); }

class SolverBase {

	public:

		SolverBase(hacc::Distribution &dist)
			: d(dist)
		{
			initialize();
		}

		virtual ~SolverBase() {
			delete m_dfft;
		}

		//Forward
		
		void forward_solve(complex_t const *rho){
			m_dfft->forward(rho);
		}

		void forward_solve(std::vector<complex_t> const & rho){
			forward_solve(&rho[0]);
		}

		//Backward

		void backward_solve_only(complex_t *phi){
			int index = 0;
			for(int local_k0 = 0; local_k0 < m_local_dim[0]; ++local_k0){
				for(int local_k1 = 0; local_k1<m_local_dim[1]; ++local_k1){
					for(int local_k2=0; local_k2 < m_local_dim[2]; ++local_k2){
						m_buf1[index] = m_buf2[index];
						index++;
					}
					//index += d.m_d.padding[2];
				}
				//index += d.m_d.padding[1];
			}
			m_dfft->backward(phi);
		}

		void backward_solve(complex_t *phi){
			//buf2->buf1
			kspace_solve(&m_buf2[0], &m_buf1[0]);
			
			//buf1 --> phi
			m_dfft->backward(phi);
		}

		void backward_solve(std::vector<complex_t> & phi){
			backward_solve(&phi[0]);
		}

		//Both
		void solve(const complex_t *rho, complex_t *phi){
			forward_solve(rho);
			backward_solve(phi);
		}

		void solve(std::vector<complex_t> const & rho, std::vector<complex_t> & phi){
			solve(&rho[0], &phi[0]);
		}

		void shift(complex_t *rho){
			forward_solve(rho);
			if(d.self() == 0){
				m_buf2[0] = 0.0;
			}
			backward_solve_only(rho);
		}

		void power_spectrum(Bins1D& pkbins, Modes3D & pkmodes, int los = -1){
			PowerSpectrum pk(*m_dfft);
			pk.power_spectrum(&m_buf2[0], pkbins, pkmodes, los);
		}

		//Initialization
		
		void initialize(){
			m_greens_functions_initialized = false;

			m_buf1.resize(d.local_size());
			m_buf2.resize(d.local_size());
			m_buf3.resize(d.local_size());

			m_dfft = new hacc::Dfft(d,
					&m_buf2[0], //forward_output = m_buf2
					&m_buf1[0], //forward_scratch = m_buf1
					&m_buf1[0], //backward_input = m_buf1
					&m_buf3[0]);  //backward_scratch = m_buf3
			for(int i=0; i<3; i++){
				m_local_dim[i] = m_dfft->local_ng_kspace(i);
				m_self_coord[i] = m_dfft->self_kspace(i);
			}
		}
	
		//TESTING
		void set_buf(std::vector<complex_t> buf, int index){
			if(index==1){
				m_buf1.swap(buf);
			} else if(index==2) {
				m_buf2.swap(buf);
			} else if(index==3){
				m_buf3.swap(buf);
			}
		}

		//////////////////////////////////////////////////
		//
		//Solve for potential by applying G function to
		//kspace density
		//	rho ---------------- density (input)
		//	phi ---------------- potential (output)
		//
		//////////////////////////////////////////////////

		void kspace_solve(const complex_t *rho, complex_t *phi){
		
			initialize_greens_function();
			int index = 0;
			for(int local_k0 = 0; local_k0 < m_local_dim[0]; ++local_k0){
				
				for(int local_k1 = 0; local_k1<m_local_dim[1]; ++local_k1){
					
					for(int local_k2=0; local_k2 < m_local_dim[2]; ++local_k2){
						phi[index] = m_green[index]*rho[index];				
						index++;
					}
					//index += d.m_d.padding[2];
				}
				//index += d.m_d.padding[1];
			}
		}

		virtual void initialize_greens_function() = 0;

		hacc::Distribution & get_d() {return d;}

	protected:
			
		double max(double a, double b) {return a>b ? a : b;}

		std::vector<double> m_green;
		std::vector<double> m_gradient;

		std::vector<complex_t> m_buf1;
		std::vector<complex_t> m_buf2;
		std::vector<complex_t> m_buf3;

		bool m_greens_functions_initialized;

		hacc::Distribution &d;
		hacc::Dfft *m_dfft;

		int m_local_dim[3];
		int m_self_coord[3];
};

class SolverDiscrete : public SolverBase {

	public:

		SolverDiscrete(hacc::Distribution &dist) :
			SolverBase(dist)
		{
		}

		void initialize_greens_function(){
			double kstep;
			int index;
			int k[3];
			std::vector<double> cosine;
	
			int ng = d.m_d.n[0];
			const int nq = ng/2;

			if(m_greens_functions_initialized){
				return;
			}

			m_greens_functions_initialized = true;

			m_green.resize(d.local_size());
			m_gradient.resize(d.m_d.n[0]);
			cosine.resize(d.m_d.n[0]);

			double tpi = 2.0*atan(1.0)*4.0;
			double pi = tpi/2.0;
			kstep = tpi / (double) d.m_d.n[0];
			double kcoeff = tpi / (1.0*ng);
		/*
			for(int kk=0; kk < d.m_d.n[0]; ++kk) {
				cosine[kk] = cos(kk*kstep);
				//m_gradient[kk] = sin(kk*kstep);
			}
		*/	
			index=0;
			//double coeff = 0.5/double(d.global_size());
			double coeff = 1.0/double(d.global_size()); //NOTE: This is to account for forward+backward FFT
			for(int local_k0 = 0; local_k0 < m_local_dim[0]; ++local_k0){
				k[0] = local_k0 + m_self_coord[0]*m_local_dim[0];
				if(k[0]>=nq){ k[0] = -MOD(ng-k[0], ng); }
				for(int local_k1 = 0; local_k1<m_local_dim[1]; ++local_k1){
				
					k[1] = local_k1 + m_self_coord[1]*m_local_dim[1];	
					if(k[1]>=nq){ k[1] = -MOD(ng-k[1], ng); }
					for(int local_k2=0; local_k2 < m_local_dim[2]; ++local_k2){
						k[2] = local_k2 + m_self_coord[2]*m_local_dim[2];
						if(k[2]>=nq){ k[2] = -MOD(ng-k[2], ng); }
						double kk = kcoeff*sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
						//m_green[index] = 0.5*coeff / (cosine[k[0]] + cosine[k[1]] + cosine[k[2]] - 3.0);
						m_green[index] = -coeff/(kk*kk); 
						index++;
					}
				}
			}
			if(d.self() == 0){
				m_green[0] = 0.0;
			}
		}

};

#endif
