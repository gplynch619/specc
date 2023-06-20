/*
	This file contains routines needed to specify the cosmology for initial conditions.
	This includes routines such as calculating the growth function, storing the transfer function
	from CAMB, and also contains a struct with easily accessible cosmo parameters.


	See $(HACC_TRUNK)/nbody/initializer/InitCosmology.cxx for comparison.

	Gabe Lynch, Aug 2019
	lynchg@anl.gov
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

#ifdef GSL
#include <gsl/gsl_sf_hyperg.h> // Only needed when computing closed form, Lamba=0 growth fxn. 
#endif

#include "CosmoClass.h"
#include "BasicDefinitions.h"
#include "initread.h" //for reading transfer func

#define TF_RELATIVE_TOLERANCE 0.00001

CosmoClass::CosmoClass(REAL Omega_m,
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
				REAL scatter) {
	CCP.Omega_m = Omega_m;
	CCP.Omega_cdm = Omega_cdm;
	CCP.Omega_bar = Omega_bar;
	CCP.Omega_cb = Omega_cb;
	CCP.Omega_nu = Omega_nu;
	CCP.f_nu_massless = f_nu_massless;
	CCP.f_nu_massive = f_nu_massive;
	CCP.Omega_r = Omega_r;
	CCP.h = h;
	CCP.w_de = w_de;
	CCP.wa_de = wa_de;
	CCP.cdm_mass = m;
	CCP.k1coeff = CCP.setKappa1coeff(m);
	CCP.k2coeff = CCP.setKappa2coeff(m, scatter);
}

//kk = k = sqrt(k dot k)
void CosmoClass::GrowthFactor(REAL z, REAL kk,  REAL &gf, REAL &g_dot){
	REAL x1, x2, dplus, ddot;
	const REAL zinfinity = 100000.0;

	CCP.setKappa1(kk);
	CCP.setKappa2(kk);

	x1 = 1.0/(1.0+zinfinity);
	x2 = 1.0/(1.0+z);
	REAL ystart[2];
	ystart[0] = x1;
	ystart[1] = 0.0;

	odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::growths, false);

	dplus = ystart[0];
	ddot = ystart[1];
	x1 = 1.0/(1.0+zinfinity);
	x2 = 1.0;
	ystart[0] = x1;
	ystart[1] = 0.0;

	odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::growths, false);

	gf = dplus/ystart[0];
	g_dot = ddot/ystart[0];
}

///////////////////////////////////////////////////////
// This function calculates the derivative of the 
// growth factor for use in velocity initialization.
// it is normalized such that gdot(1)=1.0
///////////////////////////////////////////////////////

//kk should be dimensionful

void CosmoClass::gdot(REAL z, REAL kk,  REAL &g_dot){
	REAL norm, gd;

	REAL l2 = CCP.k1coeff;

	REAL pi = 4.0*atan(1.0);

	REAL a = 1.0/(1+z);
	REAL arg = l2*kk*kk/sqrt(a);

	REAL H = Hubble(z)*100.0*CCP.h;
	REAL adot = H*a;

	//Bessel functions J_{-3/2} and J_{-5/2}
	REAL J3 = sqrt(2.0/(pi*arg))*(-cos(arg)/arg - sin(arg));
	REAL J5 = sqrt(2.0/(pi*arg))*(3.0*cos(arg)/(arg*arg)+3.0*sin(arg)/arg 
			- cos(arg));

	REAL term1 = 0.5*l2*kk*kk*adot*powf(a,-7.0/4.0)*J3;
	REAL term2 = adot*powf(a, -5.0/4.0)*J5;

	gd=term1+term2;

	//Now we must normalize it by D(k,1)

	arg = l2*kk*kk;
	norm = sqrt(2.0/(pi*arg))*(3.0*cos(arg)/(arg*arg)+3.0*sin(arg)/arg);

	g_dot = gd/norm;

}

void CosmoClass::MeanGrowth(REAL z, REAL kk, REAL *Lk){

	REAL a = 1.0/(1+z);

	REAL mratio = CCP.cdm_mass/1e-24;
	REAL omega_ratio = 1.0; //Omega_a / Omega_cdm
	REAL jeans = 66.5*powf(a,0.25)*powf(CCP.cdm_mass/1e-22, 0.5)
		*powf(CCP.Omega_cdm*CCP.h*CCP.h/0.12, 0.25)/CCP.h;

	REAL alpha = 0.194*powf(mratio, -0.501)*powf(omega_ratio, 0.0829)*CCP.h; 
	// alpha in units pf Mpc/h
	
	REAL k0 = 0.0334*powf(mratio, -0.00485)*powf(omega_ratio, 0.527)*jeans;

	REAL result = (1.0 - powf(1.0+exp(-2*alpha*(kk-k0)), -8));

	*Lk = result;
}

REAL CosmoClass::EvolveDelta(REAL z, REAL kk, REAL Pk){
	REAL x1, x2, pkf; 
	const REAL zinfinity = 100000.0;

	CCP.setKappa1(kk);
	CCP.setKappa2(kk);

	x1 = 1.0/(1.0+zinfinity);
	x2 = 1.0/(1.0+z);
	
	REAL ystart[2];
	ystart[0] = sqrt(Pk)*x1;
	ystart[1] = 0.0;

	odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::growths, false);

	pkf=ystart[0]*ystart[0];
	
	return pkf;
}

#ifdef GSL
void CosmoClass::MatterGrowthFactor(REAL z, REAL kk, REAL &gf){
	double mu,j,xi;
	double alpha,b,arg;
	REAL a = 1.0/(1.0+z);
	
	CCP.setKappa1(kk);
	CCP.setKappa2(kk);
	//REAL a=z;
	xi = sqrt(6.0*CCP.kappa2);

	j=1.5*CCP.kappa1/xi;
	mu=-5.0/4.0;
	arg=xi/a;
	
	alpha=mu-j+0.5;
	b=2.0*mu+1;
	
	double norm = pow(1.0, 0.25)*exp(-xi/2.0)*pow(xi, mu+0.5)*gsl_sf_hyperg_U(alpha, b, xi); //growth today
	//double norm = 1.0;
	gf = pow(a,0.25)*exp(-arg/2.0)*pow(arg, mu+0.5)*gsl_sf_hyperg_U(alpha,b,arg)/norm;
}

void CosmoClass::MatterGrowthDeriv(REAL z, REAL kk, REAL &gdot){
	double mu, j, xi;
	double alpha, b, arg;
	double term1, term2;
	REAL a = 1.0/(1.0+z);
	
	CCP.setKappa1(kk);
	CCP.setKappa2(kk);
	
	xi = sqrt(6.0*CCP.kappa2); //sqrt(6kappa2)
	j = 1.5*CCP.kappa1/xi;
	mu = 5.0/4.0;
	arg = xi/a;

	alpha = mu-j+0.5;
	b=2.0*mu+1.0;

	double norm = pow(1, 0.25)*exp(-xi/2.0)*pow(xi, mu+0.5)*gsl_sf_hyperg_U(alpha, b, xi); //D(1) growth today

	term1 = 0.25*pow(a, -0.75)*exp(-arg/2.0)*pow(arg, mu+0.5)*gsl_sf_hyperg_U(alpha, b, arg);
	term2 = (-0.5*exp(-arg/2.0)*pow(arg, mu+0.5) + 
			exp(-arg/2.0)*(mu+0.5)*pow(arg, mu-0.5))*gsl_sf_hyperg_U(alpha, b, arg);
	term2 += exp(-arg/2.0)*pow(arg, mu+0.5)*(-alpha)*gsl_sf_hyperg_U(alpha+1, b+1, arg);
	term2 *= -pow(a, 0.25)*xi/(a*a);

	gdot = (term1+term2)/norm;

}
#endif

void CosmoClass::GrowthFactorLCDM(REAL z, REAL *gf, REAL *g_dot){
	REAL x1, x2, dplus, ddot;
	const REAL zinfinity = 100000.0;

	x1 = 1.0/(1.0+zinfinity);
	x2 = 1.0/(1.0+z);
	REAL ystart[2];
	ystart[0] = x1;
	ystart[1] = 0.0;

	odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::LCDMgrowths, false);

	dplus = ystart[0];
	ddot = ystart[1];
	x1 = 1.0/(1.0+zinfinity);
	x2 = 1.0;
	ystart[0] = x1;
	ystart[1] = 0.0;

	odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, &CosmoClass::LCDMgrowths, false);

	*gf = dplus/ystart[0];
	*g_dot = ddot/ystart[0];
}

REAL CosmoClass::CosmoClassParams::Omega_nu_massive(REAL a){
	REAL mat = Omega_nu/pow(a,3.0f);
	REAL rad = f_nu_massive*Omega_r/pow(a,4.0f);
	return (mat>=rad)*mat + (rad>mat)*rad;
}

//Computes the first small coefficient kappa1

double CosmoClass::CosmoClassParams::setKappa1coeff(REAL mass){
	REAL alpha=NU*h; //units of ev*Mpc**2
	//REAL coeff = (1/6.0)*(alpha/mass)*(alpha/mass);
	REAL coeff = alpha/mass;
	return coeff; //coeff = l_n^2
}

//Computes second growth function coefficient kappa 2. Returns result in units of Mpc**2
double CosmoClass::CosmoClassParams::setKappa2coeff(REAL mass, REAL scatter){
	if(scatter==0){ return 0; }
	REAL c2 = VLIGHT*VLIGHT;
	double coeff = c2*HBAR;
	coeff*=HBAR*c2;
	coeff*=c2*scatter/(mass*mass*mass*GRAV_C);
	coeff*=M2MPC*M2MPC;
	return coeff;
}

double CosmoClass::CosmoClassParams::setKappa1(REAL kk){
	kappa1 = k1coeff*kk*kk*kk*kk;
}

double CosmoClass::CosmoClassParams::setKappa2(REAL kk){
	kappa2 =  k2coeff*kk*kk;
}

void CosmoClass::LCDMgrowths(REAL a, REAL y[], REAL dydx[]){
	REAL H;
	H = sqrt(CCP.Omega_cb/pow(a, 3.0f) + (1.0 + CCP.f_nu_massless)*CCP.Omega_r/pow(a,4.0f)
			+ CCP.Omega_nu_massive(a) 
			+ (1.0 - CCP.Omega_m - (1.0+CCP.f_nu_massless)*CCP.Omega_r)
			*pow(a, (-3.0*(1.0+CCP.w_de+CCP.wa_de)))*exp(-3.0*CCP.wa_de*(1.0-a))
			);
	dydx[0] = y[1]/(a*H);
	dydx[1] = -2.0*y[1]/a + 1.5*CCP.Omega_cb*y[0]/(H*pow(a, 4.0f));
}

void CosmoClass::growths(REAL a, REAL y[], REAL dydx[]){
	REAL H;
	H = sqrt(CCP.Omega_cb/pow(a, 3.0f) + (1.0 + CCP.f_nu_massless)*CCP.Omega_r/pow(a,4.0f)
			+ CCP.Omega_nu_massive(a) 
			+ (1.0 - CCP.Omega_m - (1.0+CCP.f_nu_massless)*CCP.Omega_r)
			*pow(a, (-3.0*(1.0+CCP.w_de+CCP.wa_de)))*exp(-3.0*CCP.wa_de*(1.0-a))
			);
	dydx[0] = y[1]/(a*H);
	//dydx[1] = -2.0*y[1]/a + (1.5*CCP.Omega_cb*y[0]/(H*pow(a, 4.0f)))*(1.0 - CCP.kappa1/a - CCP.kappa2/(a*a));
	dydx[1] = -2.0*y[1]/a + (1.5*CCP.Omega_cb*y[0]/(H*pow(a, 4.0f)))*(1.0 - CCP.kappa1/a); 
}

REAL CosmoClass::Hubble(REAL z_in){
	REAL H;
	REAL a = 1.0/(1.0+z_in);
	H = sqrt(CCP.Omega_cb/pow(a, 3.0f) + (1.0 + CCP.f_nu_massless)*CCP.Omega_r/pow(a,4.0f)
			+ CCP.Omega_nu_massive(a) 
			+ (1.0 - CCP.Omega_m - (1.0+CCP.f_nu_massless)*CCP.Omega_r)
			*pow(a, (-3.0*(1.0*CCP.w_de+CCP.wa_de)))*exp(-3.0*CCP.wa_de*(1.0-a))
			);
	return H;
}

TransferClass::TransferClass(const CosmoClass::CosmoClassParams &CCP,
		REAL n_s_,
		REAL Sigma_8_,
		int TFFlag_,
		std::string tfName,
		bool printWarningsFromThisRank,
		std::string tfInName) :
	n_s(n_s_),
	Sigma_8(Sigma_8_),
	R_M(8.0),
	TFFlag(TFFlag_),
	Pk_norm(1.0),
	cobe_temp(2.728),
	tt(cobe_temp/2.7*cobe_temp/2.7),
	kh_tmp(-1.0),
	last_k(-1.0),
	table_size(-1),
	hubble(CCP.h)
{
	REAL Omega_m = CCP.Omega_m;
	REAL Omega_cdm = CCP.Omega_cdm;
	REAL Omega_bar = CCP.Omega_bar;
	REAL Omega_nu = CCP.Omega_nu;
	REAL h = CCP.h;

	if(TFFlag == 0){ //this will be true by default, but this is added now for future extensibility
		std::stringstream inputFile;
		if(!getRank0Stream(tfName.c_str(), inputFile)){
			std::ostringstream ost;
			ost << "TransferClass: cannot open '" << tfName <<"'";
			throw std::runtime_error(ost.str());
		}

		long ln = 0;	
		if(printWarningsFromThisRank == true){ 
			std::cout<<"IN NAME: "<<tfName.c_str()<<std::endl; 
		}
		REAL inputLine[TF_NUM_COLUMNS];
		while(! inputFile.eof()){
			for(int j=0; j<TF_NUM_COLUMNS; j++){
				inputFile >> inputLine[j]; 
				//Takes first line and fills input line column by column
				//So inputline holds the current line indexed by column
			}
			
			table_kk.push_back(inputLine[TF_K]);
			table_tf_cdm.push_back(inputLine[TF_CDM]);
			table_tf_bar.push_back(inputLine[TF_BAR]);
			table_tf_gamma.push_back(inputLine[TF_GAMMA]);
			table_tf_nu_massless.push_back(inputLine[TF_NU_MASSLESS]);
			table_tf_nu_massive.push_back(inputLine[TF_NU_MASSIVE]);
			table_tf_total.push_back(fabs(inputLine[TF_TOTAL]));

			// tfbar*Omega_bar/Omega_m + tfcfm(Omega_m -Omega_bar-Omega_nu)/Omega_m + tfnu*Omega_nu/Omega_m
			table_tf_cb.push_back(fabs(table_tf_bar[ln]*Omega_bar/Omega_m 
					+ table_tf_cdm[ln]*Omega_cdm/Omega_m));

			table_tf_my_total.push_back(table_tf_cb[ln] + table_tf_nu_massive[ln]*Omega_nu/Omega_m);

			++ln;
		}
		table_size = ln-1;
		last_k = table_kk[table_size-1];

		//Check that it matches CAMB's total transfer function
		bool tfmismatch = false;
		for(int i=0; i<table_size; i++){
			if( fabs(table_tf_my_total[i] - table_tf_total[i])/table_tf_total[i] >= TF_RELATIVE_TOLERANCE){
				tfmismatch = true;
			}
		}

		if(tfmismatch != false && printWarningsFromThisRank ==true){
			std::cout<< "WARNING: Mismatch between read and computed total TF" <<std::endl;
		}

		tan_f_cb = &TransferClass::camb_cb;
		tan_f_tot = &TransferClass::camb_tot;
	} else if (TFFlag==1){
		// Klypin-Holtzmann tf

		if(printWarningsFromThisRank == true){ std::cout<<"Using Klypin-Holtzmann TF"<<std::endl; }
    	
		REAL akh1=pow(46.9*Omega_m*h*h, 0.670)*(1.0+pow(32.1*Omega_m*h*h, -0.532));
   		REAL akh2=pow(12.0*Omega_m*h*h, 0.424)*(1.0+pow(45.0*Omega_m*h*h, -0.582));
    	REAL alpha=pow((double) akh1, -1.0*Omega_bar/Omega_m)*pow((double) akh2, 
				pow(-1.0*Omega_bar/Omega_m, 3.0));
    	kh_tmp = Omega_m*h*sqrt(alpha)*pow(1.0-Omega_bar/Omega_m, 0.6);

		//h was removed due to KH expecting physical k, and I pass k/h.

    	tan_f_cb = &TransferClass::klypin_holtzmann;
    	tan_f_tot = &TransferClass::klypin_holtzmann;

    	last_k = 200.0;
	} else if (TFFlag==2){
		if(tfInName.empty()){
			std::ostringstream ost;
			ost<<"AxionCAMB ICs specified but no z=zin Transfer file supplied!";
			throw std::runtime_error(ost.str());
		}
		std::stringstream inputFile;
		if(!getRank0Stream(tfInName.c_str(), inputFile)){
			std::ostringstream ost;
			ost << "TransferClass: cannot open '" << tfInName <<"'";
			throw std::runtime_error(ost.str());
		}

		long ln = 0;
		REAL inputLine[TF_NUM_COLUMNS]; //z=zin
		while(! inputFile.eof()){
			for(int j=0; j<TF_NUM_COLUMNS; j++){
				inputFile >> inputLine[j];
				//Takes first line and fills input line column by column
				//So inputline holds the current line indexed by column
			}
			//kk will hold the values of k the Power spectrum is sampled at
			//tf_cdm, in this case, is really "pk_fdm"
			table_kk.push_back(inputLine[TF_K]);
			table_tf_cdm.push_back(fabs(inputLine[TF_TOTAL])); 
			//This is TF_TOTAL because we want the total transfer function at z=zin. We are storing it in
			//table_tf_cdm for convenience.
			table_tf_bar.push_back(inputLine[TF_BAR]);
			table_tf_gamma.push_back(inputLine[TF_GAMMA]);
			table_tf_nu_massless.push_back(inputLine[TF_NU_MASSLESS]);
			table_tf_nu_massive.push_back(inputLine[TF_NU_MASSIVE]);
			++ln;
		}
		table_size = ln-1;
		last_k = table_kk[table_size-1];
		
		inputFile.str("");
		if(!getRank0Stream(tfName.c_str(), inputFile)){
			std::ostringstream ost;
			ost << "TransferClass: cannot open '" << tfName <<"'";
			throw std::runtime_error(ost.str());
		}

		
		REAL inputLine1[TF_NUM_COLUMNS];
		while(! inputFile.eof()){
			for(int j=0; j<TF_NUM_COLUMNS; j++){
				inputFile >> inputLine1[j]; 
			}
			table_tf_total.push_back(fabs(inputLine1[TF_TOTAL]));
		}

		tan_f_cb = &TransferClass::axion_camb_cb;
		tan_f_tot = &TransferClass::camb_tot;
		
		if(printWarningsFromThisRank ==true){
			std::cout<< "Initial Power Spectrum read in from AxionCAMB" <<std::endl;
		}
	
	}


	SetSigma8(Sigma_8);
	if(printWarningsFromThisRank == true){ std::cout<<"Pk_norm = "<<Pk_norm<<std::endl; }
	
	return;
}

//Transfer Function Routines

REAL TransferClass::camb_cb(REAL k){
	REAL t_f = interpolate(&table_kk[0], &table_tf_cb[0], table_size, k);
	return(t_f);
}

REAL TransferClass::axion_camb_cb(REAL k){
	//reusing tables, table_tf_cdm is really table_pk_fdm
	REAL t_f = interpolate(&table_kk[0], &table_tf_cdm[0], table_size, k);
	return(t_f);
}

REAL TransferClass::camb_tot(REAL k){
	REAL t_f = interpolate(&table_kk[0], &table_tf_total[0], table_size, k);
	return(t_f);
}

REAL TransferClass::klypin_holtzmann(REAL k){
  if (k == 0.0) return(0.0);
  REAL qkh = k*tt/kh_tmp;
  
  // NOTE: the following line has 0/0 for k=0.
  // This was taken care of at the beginning of the routine.
  REAL t_f = log(1.0 + 2.34*qkh)/(2.34*qkh) *
    pow(1.0
	+ 13.0*qkh
	+ pow(10.5*qkh, 2.0)
	+ pow(10.4*qkh, 3.0)
	+ pow(6.51*qkh, 4.0),
	-0.25);
  return(t_f);
}

REAL TransferClass::Tf_cb(REAL k){
	REAL tf = ((this->*tan_f_cb)(k));
	return (tf);
}

REAL TransferClass::Pk_cb(REAL k){
	REAL tf = ((this->*tan_f_cb)(k));
	return Pk_norm*pow(k, n_s)*tf*tf; //normalized to today
}

REAL TransferClass::Pk_total(REAL k){
	REAL tf = ((this->*tan_f_tot)(k));
	return Pk_norm*pow(k, n_s)*tf*tf; //This assumes growth at 1/k = 8 h^-1 Mpc is the same as CDM
}

//Sigma normalization routines

REAL TransferClass::sigma2(REAL k) {
	REAL w_f = 3.0*(sin(k*R_M) - k*R_M*cos(k*R_M))/pow((double) k*R_M, 3.0);
	REAL s2 = k*k*w_f*w_f*Pk_total(k);
	return(s2);
}

void TransferClass::SetSigma8(REAL s8){
	R_M = 8.0;
	REAL pi = 4.0*atan(1.0);
	const REAL k_min=0.0, k_max=last_k;
	REAL s2 = 1.0/(2.0*pi*pi) * integrate(&TransferClass::sigma2, k_min, k_max);
	Pk_norm = s8*s8/s2;
	return;
}

//////////////////////////////////////////
//
//Numerical Recipes
//
//////////////////////////////////////////

#define MAXSTP 10000
#define TINY 1.0e-30
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void CosmoClass::odesolve(REAL ystart[], int nvar, REAL x1, REAL x2,
			  REAL eps, REAL h1,
                          void (CosmoClass::*derivs)(REAL, REAL [], REAL []),
			  bool print_stat)
{
  int i, nstp, nok, nbad, feval;
  REAL x,hnext,hdid,h;
  REAL *yscal,*y,*dydx;
  const REAL hmin=0.0;

  feval = 0;
  yscal= (REAL *)malloc(nvar*sizeof(REAL));
  y= (REAL *)malloc(nvar*sizeof(REAL));
  dydx= (REAL *)malloc(nvar*sizeof(REAL));
  x=x1;
  h=SIGN(h1,x2-x1);
  nok = nbad = 0;
  for (i=0; i<nvar; ++i) {y[i]=ystart[i];}

  for (nstp=0; nstp<MAXSTP; ++nstp) {
    (this->*derivs)(x, y, dydx);
    ++feval;
    for (i=0; i<nvar; ++i)
    {yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;}
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,&feval,derivs);
    if (hdid == h) ++nok; else ++nbad;
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=0; i<nvar; ++i) {ystart[i]=y[i];}
      free(dydx);
      free(y);
      free(yscal);
      if (print_stat){
	printf("ODEsolve:\n");
	printf(" Evolved from x = %f to x = %f\n", x1, x2);
	printf(" successful steps: %d\n", nok);
	printf(" bad steps: %d\n", nbad);
	printf(" function evaluations: %d\n", feval);
      }
      return;
    }
    if (fabs(hnext) <= hmin) {
      printf("Step size too small in ODEsolve");
      exit(1);
    }
    h=hnext;
  }
  printf("Too many steps in ODEsolve");
  exit(1);
}
#undef MAXSTP
#undef TINY
#undef SIGN

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
static REAL maxarg1,maxarg2, minarg1, minarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
void CosmoClass::rkqs(REAL y[], REAL dydx[], int n, REAL *x, REAL htry,
		      REAL eps,
                      REAL yscal[], REAL *hdid, REAL *hnext, int *feval,
                      void (CosmoClass::*derivs)(REAL, REAL [], REAL []))
{
  int i;
  REAL errmax,h,htemp,xnew,*yerr,*ytemp;

  yerr= (REAL *)malloc(n*sizeof(REAL));
  ytemp= (REAL *)malloc(n*sizeof(REAL));
  h=htry;

  for (;;) {
    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    *feval += 5;
    errmax=0.0;
    for (i=0; i<n; ++i) {errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));}
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp=SAFETY*h*pow((double) errmax,PSHRNK);
    h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
    xnew=(*x)+h;
    if (xnew == *x) {
      printf("Stepsize underflow in ODEsolve rkqs");
      exit(1);
    }
  }
  if (errmax > ERRCON) *hnext=SAFETY*h*pow((double) errmax,PGROW);
  else *hnext=5.0*h;
  *x += (*hdid=h);
  for (i=0; i<n; ++i) {y[i]=ytemp[i];}
  free(ytemp);
  free(yerr);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef FMAX
#undef FMIN

void CosmoClass::rkck(REAL y[], REAL dydx[], int n, REAL x, REAL h,
		      REAL yout[], REAL yerr[],
		      void (CosmoClass::*derivs)(REAL, REAL [], REAL []))
{
  int i;
  static REAL a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  REAL dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  REAL *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

  ak2= (REAL *)malloc(n*sizeof(REAL));
  ak3= (REAL *)malloc(n*sizeof(REAL));
  ak4= (REAL *)malloc(n*sizeof(REAL));
  ak5= (REAL *)malloc(n*sizeof(REAL));
  ak6= (REAL *)malloc(n*sizeof(REAL));
  ytemp= (REAL *)malloc(n*sizeof(REAL));

  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (this->*derivs)(x+a2*h,ytemp,ak2);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (this->*derivs)(x+a3*h,ytemp,ak3);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (this->*derivs)(x+a4*h,ytemp,ak4);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (this->*derivs)(x+a5*h,ytemp,ak5);
  for (i=0; i<n; ++i)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (this->*derivs)(x+a6*h,ytemp,ak6);
  for (i=0; i<n; ++i)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0; i<n; ++i)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

  free(ytemp);
  free(ak6);
  free(ak5);
  free(ak4);
  free(ak3);
  free(ak2);
}

// Numerical integration routines from Numerical Recipes
// calling:
// x = integrate(&CosmoClass::func, a, b)
//
// where
//    func -- is the (REAL) function whose integral is calculated,
//                               and has to be member of the CosmoClass
//    a    -- is the (REAL) lower boundary for integration
//    b    -- is the (REAL) upper boundary for integration

#define FUNC(x) ((this->*func)(x))
REAL TransferClass::midpoint(REAL (TransferClass::*func)(REAL),
			     REAL a, REAL b, long n)
{
  REAL x,tnm,sum,del,ddel;
  static REAL s;
  long it,j;

  if (n == 1) {
    return (s=(b-a)*FUNC(0.5*(a+b)));
  }
  else {
    it = 1;
    for(j=1; j<n-1; ++j) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1; j<=it; ++j){
      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return (s);
  }
}
#undef FUNC

#define EPS 1.0e-4
#define JMAX 20
#define JMIN 5
REAL TransferClass::integrate(REAL (TransferClass::*func)(REAL), REAL a, REAL b){
  REAL sol, old;
  long j;

  sol = 0.0;
  old = -1.0e26;
  for (j=1; j<=JMAX; ++j) {
    sol=midpoint(func, a, b, j);
    if (j > JMIN){
      if (fabs(sol-old) < EPS*fabs(old) ||
	  (sol==0.0 && old==0.0)) return(sol);
    }
    old = sol;
  }
  printf("integrate: no convergence!");
  abort();
  return 0;
}
#undef EPS
#undef JMAX
#undef JMIN

// Linear interpolation routine
// calling:
// y = interpolate(xx, yy, n, x)
//
// where
// xx -- is the (REAL) x array of ordered data
// yy -- is the (REAL) y array
// n  -- size of above arrays
// x  -- is the (REAL) point whose y value should be interpolated

REAL TransferClass::interpolate(REAL xx[], REAL yy[], unsigned long n, REAL x){
  REAL y, dx, dy;
  unsigned long jlo;

  locate(xx, n, x, &jlo);
  // Linear interpolation:
  dx = xx[jlo] - xx[jlo+1];
  dy = log10(yy[jlo]) - log10(yy[jlo+1]);
  y = dy/dx*(x-xx[jlo]) + log10(yy[jlo]);
  y=pow(10.0,y);
/*
  locate(xx,n,x,&jlo);

  dx = xx[jlo] - xx[jlo+1];
  dy = yy[jlo]-yy[jlo+1];
  y = dy/dx*(x-xx[jlo]) + yy[jlo];
 */ 
  return(y);
}

void TransferClass::locate(REAL xx[], unsigned long n, REAL x, unsigned long *j){
  unsigned long ju,jm,jl;
  int ascnd;

  jl=0;
  ju=n-1;
  ascnd=(xx[n-1] >= xx[0]);
  while (ju-jl > 1) {
    jm=(ju+jl)/2;
    if ((x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x == xx[0]) *j=0;
  else if(x == xx[n-1]) *j=n-2;
  else *j=jl;
  return;
}
