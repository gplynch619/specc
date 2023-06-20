#ifndef COSMOCLASS_H
#define COSMOCLASS_H

#include <vector>
#include <string>
#include "sctypes.h"

enum TF_COLUMNS{
	TF_K,
	TF_CDM,
	TF_BAR,
	TF_GAMMA,
	TF_NU_MASSLESS,
	TF_NU_MASSIVE,
	TF_TOTAL,
	TF_NUM_COLUMNS,
};

enum PK_COLUMNS{
	PK_K,
	PK_AX,
	PK_NUM_COLUMNS,
};

class CosmoClass {

	public:
		
		CosmoClass(REAL Omega_m,
				REAL Omega_cdm,
				REAL Omega_bar,
				REAL Omega_cb,
				REAL Omega_nu,
				REAL f_nu_massless,
				REAL f_nu_massive,
				REAL Omega_r,
				REAL h,
				REAL w_de,
				REAL wa_de,
				REAL m,
				REAL scatter);
		~CosmoClass(){}

		void GrowthFactor(REAL z, REAL kk, REAL &gf, REAL &g_dot);
		void MatterGrowthFactor(REAL z, REAL kk, REAL &gf);
		void MatterGrowthDeriv(REAL z, REAL kk, REAL &gdot);
		void GrowthFactorLCDM(REAL z, REAL *gf, REAL *g_dot);
		REAL Hubble(REAL z_in);
		REAL EvolveDelta(REAL zin, REAL kk, REAL pk); 

		void gdot(REAL z, REAL kk,  REAL &g_dot);
		void MeanGrowth(REAL z, REAL kk,  REAL *Lk);

		struct CosmoClassParams {
			REAL Omega_m;
			REAL Omega_cdm;
			REAL Omega_bar;
			REAL Omega_cb;
			REAL Omega_nu;
			REAL f_nu_massless;
			REAL f_nu_massive;
			REAL Omega_r;
			REAL h;
			REAL w_de;
			REAL wa_de;
			REAL cdm_mass;
			REAL k1coeff;
			REAL k2coeff;
			REAL kappa1;
			REAL kappa2;
			REAL setKappa1(REAL kk);
			REAL setKappa2(REAL kk);
			REAL Omega_nu_massive(REAL a);
		private:
			REAL setKappa1coeff(REAL mass);
			REAL setKappa2coeff(REAL mass, REAL scatter);
			friend class CosmoClass;
		};

		const CosmoClassParams &getCosmoClassParams() { return CCP; }

	private:

		CosmoClassParams CCP;

		// The following are from Numerical Recipes
		void growths(REAL a, REAL y[], REAL dydx[]); 
		void LCDMgrowths(REAL a, REAL y[], REAL dydx[]); 

		void odesolve(REAL ystart[], int nvar, REAL x1, REAL x2, REAL eps, REAL h1,
				void (CosmoClass::*derivs)(REAL, REAL [], REAL[]),
				bool print_stat);

		void rkqs(REAL y[], REAL dydx[], int n, REAL *x, REAL htry, REAL eps, REAL yscal[],
				REAL *hdid, REAL *hnext, int *feval,
				void (CosmoClass::*derivs)(REAL, REAL [], REAL []));

  		void rkck(REAL y[], REAL dydx[], int n, REAL x, REAL h,
				REAL yout[], REAL yerr[],
				void (CosmoClass::*derivs)(REAL, REAL [], REAL []));
};

//Transfer class currently only handles CAMB transfer function but this is easily (?) extendable

class TransferClass{
	public:
		TransferClass(const CosmoClass::CosmoClassParams &CCP,
				REAL n_s_,
				REAL Sigma_8_,
				int TFFLag,
				std::string tfName,
				bool printWarningsFromThisRank=true,
				std::string tfInName="");
		~TransferClass(){}

		void SetSigma8(REAL sigma8);
		REAL Pk_cb(REAL k);
		REAL Tf_cb(REAL k);
		REAL Pk_total(REAL k);
		REAL kmax() { return last_k; }

		int TF() { return TFFlag; }
	protected:
		REAL n_s;
		REAL Sigma_8;
		REAL R_M;
		REAL Pk_norm;
		const REAL cobe_temp;
		const REAL tt;
		REAL kh_tmp;
		REAL last_k;
		REAL hubble;

		// CAMB number of rows
		unsigned long table_size;
		unsigned long table_size_vel;

//////////////////////////////////////////
//Density Transfer Function Tables
//////////////////////////////////////////
		std::vector<REAL> table_kk;
		std::vector<REAL> table_tf_cdm; 
		std::vector<REAL> table_tf_bar;
		std::vector<REAL> table_tf_gamma;
		std::vector<REAL> table_tf_nu_massless;
		std::vector<REAL> table_tf_nu_massive;
		std::vector<REAL> table_tf_total;
		//constructed columns
		std::vector<REAL> table_tf_my_total;
		std::vector<REAL> table_tf_cb;
//////////////////////////////////////////

		REAL (TransferClass::*tan_f_tot)(REAL k);
		REAL (TransferClass::*tan_f_cb)(REAL k);

		REAL camb_cb(REAL k);
		REAL camb_tot(REAL k);
		REAL axion_camb_cb(REAL k);
		REAL klypin_holtzmann(REAL k);

		REAL sigma2(REAL k);

		//Numerical Recipes
		REAL integrate(REAL (TransferClass::*func)(REAL), REAL a, REAL b);
		REAL midpoint(REAL (TransferClass::*func)(REAL), REAL a, REAL b, long n);
		REAL interpolate(REAL xx[], REAL yy[], unsigned long n, REAL x);
		void locate(REAL xx[], unsigned long n, REAL x, unsigned long *j);
	private:

		int TFFlag;
};

#endif
