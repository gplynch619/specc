#include "TimeStepper.h"

TimeStepper::TimeStepper(REAL alpha_,
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
		size_t nzpr_) :
	m_aa(-1.0),
	m_tt(-1.0),
	m_pp(-1.0),
	m_zz(-1.0),
	m_alpha(-1.0),
	m_tau(-1.0),
	m_tau2(-1.0),
	m_dt(-1.0),
	m_dt2(-1.0),
	m_adot(-1.0),
	m_prefactor(-1.0),
	m_omega_matter(-1.0),
	m_omega_cdm(-1.0),
	m_omega_baryon(-1.0),
	m_omega_cb(-1.0),
	m_omega_nu(-1.0),
	m_omega_radiation(-1.0),
	m_f_nu_massless(-1.0),
	m_f_nu_massive(-1.0),
	m_ain(-1.0),
	m_afin(-1.0),
	m_pin(-1.0),
	m_pfin(-1.0),
	m_zin(-1.0),
	m_zfin(-1.0),
	m_nsteps(-1),
	m_phiscal(-1.0),
	m_fscal(-1.0),
	m_w(-1.0),
	m_wa(-1.0),
	m_currentHalfStep(0)
{
	m_alpha = alpha_;
	m_ain = ain_;
	m_afin = afin_;
	m_nsteps = nsteps_;
	m_omega_matter = omega_matter_;
	m_omega_cdm = omega_cdm_;
	m_omega_baryon = omega_baryon_;
	m_omega_cb = omega_cb_;
	m_omega_nu = omega_nu_;
	m_omega_radiation = omega_radiation_;
	m_f_nu_massless = f_nu_massless_;
	m_f_nu_massive = f_nu_massive_;
	m_w = w_;
	m_wa = wa_;

	m_pin = pow(m_ain, m_alpha);
	m_pfin = pow(m_afin, m_alpha);
	m_tin = pow(m_ain, 3.0/2.0);
	m_tfin = pow(m_afin, 3.0/2.0);
	m_zin = 1.0/m_ain - 1.0;
	m_zfin = 1.0/m_afin - 1.0;
	m_tau = (m_afin - m_ain)/(1.0*m_nsteps);
	m_tau2 = 0.5*m_tau;
	m_dt = (m_tfin - m_tin)/(1.0*m_nsteps);
	m_dt2 = 0.5*m_dt;
	m_pp = m_pin;
	m_aa = m_ain;
	m_zz = m_zin;
	m_tt = m_tin;

	m_zpr.assign(zpr_, zpr_ + nzpr_);

	for(int i=0; i<m_zpr.size(); i++){
		m_apr.push_back(1.0/(1.0+m_zpr.at(i)));
	}

	set_adot();
	set_scal();	
}

TimeStepper::~TimeStepper() {}

REAL TimeStepper::omega_nu_massive(REAL a){
	REAL mat = m_omega_nu/pow(a,3.0);
	REAL rad = m_f_nu_massive*m_omega_radiation/pow(a,4.0);
	return (mat>=rad)*mat + (rad>mat)*rad;
}

REAL TimeStepper::H_ratio(){
	return sqrt(m_omega_cb/pow(m_aa,3.0f)
			+ (1.0 + m_f_nu_massless)*m_omega_radiation/pow(m_aa, 4.0f)
			+ omega_nu_massive(m_aa)
			+ (1.0 - m_omega_matter - (1.0+m_f_nu_massless)*m_omega_radiation)
			*pow(m_aa, (float)(-3.0*(1.0+m_w+m_wa)))*exp(-3.0*m_wa*(1.0-m_aa))
			);
}

void TimeStepper::set_adot() {
	REAL pp1 = pow(m_aa, -3.0*(m_w+m_wa))*exp(-3.0*m_wa*(1.0-m_aa));
	REAL tmp = m_omega_cb
		+ (1.0+m_f_nu_massless)*m_omega_radiation/m_aa
		+ omega_nu_massive(m_aa)*pow(m_aa,3.0)
		+ (1.0-m_omega_matter-(1.0+m_f_nu_massless)*m_omega_radiation)*pp1;
	tmp /= m_aa;
	m_adot = sqrt(tmp);
	return;
}

void TimeStepper::set_scal() {
	set_adot();
	m_phiscal = 1.5*m_omega_cb/m_aa;
	//m_phiscal = 1.5*m_omega_cb;
	m_fscal = m_phiscal*m_aa/(m_alpha*m_adot*m_pp);
	return;
}

void TimeStepper::set_step(REAL da){
	if(da==0){
		m_tau = (m_afin - m_ain)/(1.0*m_nsteps);
		m_tau2 = 0.5*m_tau;
	} else {
		m_tau = da;
		m_tau2 = 0.5*m_tau;
	}
}

void TimeStepper::advanceHalfStep() {
	m_aa += m_tau2; 
	m_zz = 1.0/m_aa - 1.0;
	set_adot();
	set_scal();
	m_currentHalfStep++;
	return;
}

/*
void TimeStepper::reverseHalfStep() {
	m_pp -= m_tau2;
	m_aa = pow(m_pp, t.0/m_alpha);
	m_zz = 1.0/m_aa - 1.0;
	set_adot();
	set_scal();
	m_currentHalfStep--;
	return;
}
*/
void TimeStepper::advanceFullStep() {
	advanceHalfStep();
	advanceHalfStep();
	return;
}
/*
void TimeStepper::reverseFullStep() {
	reverseHalfStep();
	reverseHalfStep();
	return;
}*/
