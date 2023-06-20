#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <cmath>
#include <vector>
#include <iostream>

typedef double REAL;

class TimeStepper {
	public:

		TimeStepper(REAL alpha_,
				REAL ain_,
				REAL afin_,
				int nsteps_,
				REAL omega_matter_,
				REAL omega_cdm_,
				REAL omega_baryon_,
				REAL omega_cb_,
				REAL omega_nu_,
				REAL omega_radiation_,
				REAL f_nu_massless_,
				REAL f_nu_massive_,
				REAL w_,
				REAL wa_,
				REAL* zpr_,
				size_t nzpr_);
		~TimeStepper();

		void advanceHalfStep();
		void advanceFullStep();

		REAL aa() { return m_aa; }
		REAL pp() { return m_pp; }
		REAL zz() { return m_zz; }
		REAL tt() { return m_tt; }
		REAL alpha() { return m_alpha; }
		REAL tau() { return m_tau; }
		REAL tau2() { return m_tau2; }
		REAL dt() { return m_dt; }
		REAL dt2() { return m_dt2; }
		REAL adot() { return m_adot; }
		REAL prefactor() { return m_prefactor; }
		REAL H_ratio();
		REAL omega_matter() { return m_omega_matter; }
		REAL omega_cdm() { return m_omega_cdm; }
		REAL omega_baryon() { return m_omega_baryon; }
		REAL omega_cb() { return m_omega_cb; }
		REAL omega_nu() { return m_omega_nu; }
		REAL omega_radiation() { return m_omega_radiation; }
		REAL f_nu_massless() { return m_f_nu_massless; }
		REAL f_nu_massive() { return m_f_nu_massive;}
		REAL ain() { return m_ain; }
		REAL afin() { return m_afin; }
		REAL tin() { return m_tin; }
		REAL tfin() { return m_tfin; }
		REAL pin() { return m_pin; }
		REAL pfin() { return m_pfin; }
		REAL zin() { return m_zin; }
		REAL zfin() {return m_zfin; }
		int nsteps() { return m_nsteps; }
		REAL phiscal() { return m_phiscal; }
		REAL fscal() { return m_fscal; }
		REAL w() { return m_w; }
		REAL wa() { return m_wa; }
		int currentHalfStep() {return m_currentHalfStep; }
		
		std::vector<REAL> m_zpr;
		std::vector<REAL> m_apr;

		void set_step(REAL da=0);

	private:

		TimeStepper();
		TimeStepper( const TimeStepper& );
		TimeStepper& operator = (const TimeStepper& );

		REAL omega_nu_massive(REAL a);
		void set_adot();
		void set_scal();

		REAL m_aa;
		REAL m_pp;
		REAL m_zz;
		REAL m_tt;
		REAL m_alpha;
		REAL m_tau;
		REAL m_tau2;
		REAL m_dt;
		REAL m_dt2;
		REAL m_adot;
		REAL m_prefactor;

		REAL m_omega_matter;
		REAL m_omega_cdm;
		REAL m_omega_baryon;
		REAL m_omega_cb;
		REAL m_omega_nu;
		REAL m_omega_radiation;
		REAL m_f_nu_massless;
		REAL m_f_nu_massive;

		REAL m_ain;
		REAL m_afin;
		REAL m_pin;
		REAL m_pfin;
		REAL m_tin;
		REAL m_tfin;
		REAL m_zin;
		REAL m_zfin;
		int m_nsteps;


		REAL m_phiscal;
		REAL m_fscal;
		REAL m_w;
		REAL m_wa;

		int m_currentHalfStep;
};

#endif
