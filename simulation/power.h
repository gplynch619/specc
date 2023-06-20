#ifndef POWER_HPP
#define POWER_HPP

// for accumulating 1D bins and 3D modes
#include <mpi.h>

#include <vector>

#include "complex-type.h"
#include "Dfft.hpp"
#include "Distribution.hpp"


class Bins1D {

public:

  Bins1D()
    : comm(MPI_COMM_WORLD),
      nbin(0),
      xmin(0.0),
      xmax(0.0),
      logbins(false),
      dx(0.0),
      initializedq(false)
  {}

  Bins1D(MPI_Comm comm_,
	 int nbin_,
	 double xmin_,
	 double xmax_,
	 bool logbins_)
  {
    initialize(comm_, nbin_, xmin_, xmax_, logbins_);
  }

  ~Bins1D(){}

  void initialize(MPI_Comm comm_,
		  int nbin_,
		  double xmin_,
		  double xmax_,
		  bool logbins_ = false)
  {
    comm = comm_;
    nbin = nbin_;
    xmin = xmin_;
    xmax = xmax_;
    logbins = logbins_;

    if(logbins == true) {
      xmin = log(xmin);
      xmax = log(xmax);
    }
    dx = (xmax-xmin)/((double)nbin);

    xBins.resize(nbin);
    xBins.assign(nbin, 0.0);
    wBins.resize(nbin);
    wBins.assign(nbin, 0.0);
    monoBins.resize(nbin);
    monoBins.assign(nbin, 0.0);
    quadBins.resize(nbin);
    quadBins.assign(nbin, 0.0);

    initializedq = true;
  }

  inline
  void addVal(double x,
	      double mono,
	      double quad = 0.0)
  {
    if(logbins == true)
      x = log(x);
    int bin = (int)((x-xmin)/dx);
    if(bin >= 0 && bin < nbin) {
      xBins[bin] += x;
      wBins[bin] += 1.0;
      monoBins[bin] += mono;
      quadBins[bin] += quad;
    }
  }

  void accumulate()
  {
    MPI_Allreduce(MPI_IN_PLACE, &wBins[0],    nbin, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &xBins[0],    nbin, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &monoBins[0], nbin, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &quadBins[0], nbin, MPI_DOUBLE, MPI_SUM, comm);

    for(int i=0; i<nbin; i++)
      if(wBins[i] > 0.0) {
	xBins[i] /= wBins[i];
	monoBins[i] /= wBins[i];
	quadBins[i] /= wBins[i];
      }
  }

  MPI_Comm comm;
  int nbin;
  double xmin;
  double xmax;
  bool logbins;
  double dx;
  bool initializedq;

  std::vector<double> wBins;
  std::vector<double> xBins;
  std::vector<double> monoBins;
  std::vector<double> quadBins;
};



class Modes3D {

public:

  Modes3D()
    : comm(MPI_COMM_WORLD),
      nmode3d(0),
      nmodetotal(0),
      initializedq(false)
  {}

  Modes3D(MPI_Comm comm_,
	  int nmode3d_)
  {
    initialize(comm_, nmode3d_);
  }

  ~Modes3D(){}

  void initialize(MPI_Comm comm_,
		  int nmode3d_)
  {
    comm = comm_;
    nmode3d = nmode3d_;

    nmodetotal = 8*nmode3d*nmode3d*nmode3d;
    val3d.resize(nmodetotal);
    val3d.assign(nmodetotal, 0.0);
    legendre3d.resize(nmodetotal);
    legendre3d.assign(nmodetotal, 0.0);
 
    initializedq = true;
 }

  inline
  void addVal(int tk0,
	      int tk1,
	      int tk2,
	      double val,
	      double legendre=0.0)
  {
    if( nmode3d > 0 && 
	tk0 > -1*nmode3d && tk0 <= nmode3d &&
	tk1 > -1*nmode3d && tk1 <= nmode3d &&
	tk2 > -1*nmode3d && tk2 <= nmode3d ) {
      int indx = 0;
      indx += tk0 + nmode3d - 1;
      indx *= (2*nmode3d);
      indx += tk1 + nmode3d - 1;
      indx *= (2*nmode3d);
      indx += tk2 + nmode3d - 1;

      val3d[indx] = val;
      legendre3d[indx] = legendre;
    }
  }

  void accumulate()
  {
    if( nmode3d > 0 ) {
      MPI_Allreduce(MPI_IN_PLACE, &val3d[0],      nmodetotal, 
		    MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(MPI_IN_PLACE, &legendre3d[0], nmodetotal, 
		    MPI_DOUBLE, MPI_SUM, comm);
    }
  }

  MPI_Comm comm;
  int nmode3d;
  int nmodetotal;
  bool initializedq;

  std::vector<double> val3d;
  std::vector<double> legendre3d;
};



class PowerSpectrum {

public:

  /*
!-----3-D anti-CIC filter for deconvolution

      forall(ii=1:ng)kk(ii)=(ii-1)*tpi/(1.0*ng)

      do ii=1,ng
      if(ii.ge.ng/2+1)kk(ii)=(ii-ng-1)*tpi/(1.0*ng)
      enddo

      mult(1)=1.0
      
      forall(ii=2:ng)mult(ii)=
      #      1.0/(sin(kk(ii)/2.0)/(kk(ii)/2.0))**2
      
      forall(ii=1:ng,jj=1:ng,mm=1:ng)erhotr(ii,jj,mm)=
      #      mult(ii)*mult(jj)*mult(mm)*erhotr(ii,jj,mm)
  */

  PowerSpectrum(hacc::Dfft &dfft)
    : m_dfft(dfft),
      d(m_dfft.get_d()),
      ng(m_dfft.get_d().m_d.n[0]),
      m_local_dim(m_dfft.local_ng_kspace()),
      m_self_coord(m_dfft.self_kspace())
  {
    m_pk_ksq.resize(ng);
    m_pk_ksq.assign(ng, 0.0);
    m_pk_cic.resize(ng);
    m_pk_cic.assign(ng, 0.0);

    // cache periodic ksq
    double tpi = 2.0*atan(1.0)*4.0;
    for (int k = 0; k < ng / 2; ++k) {
      m_pk_ksq[k] = k * k;

      m_pk_ksq[k + ng / 2] = (k - ng / 2) * (k - ng / 2);

      double kk = tpi*k/ng;
      m_pk_cic[k] = pow(sin(0.5*kk)/(0.5*kk),-4.0);

      kk = tpi*(k-ng/2)/ng;
      m_pk_cic[k + ng/2] = pow(sin(0.5*kk)/(0.5*kk),-4.0);
    }
    m_pk_cic[0] = 1.0;    
  } 

  
  virtual ~PowerSpectrum(){}


  ///
  // calculate the k-space power spectrum
  //   P(modk) = Sum { |rho(k)|^2 : |k| = modk, k <- [0, ng / 2)^3, periodically extended }
  ///


  // auto power
  void power_spectrum(complex_t *rho,
		      Bins1D & pkbins,
		      Modes3D & pkmodes,
		      int los = -1)
  {
    power_spectrum(rho, rho, pkbins, pkmodes, los);
  }


  // cross power
  void power_spectrum(complex_t *rho1,
		      complex_t *rho2,
		      Bins1D & pkbins,
		      Modes3D & pkmodes,
		      int los = -1)
  {
    double inversevolume = pow(1.0*ng,-3.0); 
	double pi = atan(1.0)*4.0;
    
	int index = 0;
    for (int local_k0 = 0; local_k0 < m_local_dim[0]; ++local_k0) {
		int k0 = local_k0 + m_self_coord[0] * m_local_dim[0];
		double ksq0 = m_pk_ksq[k0];
		
		for (int local_k1 = 0; local_k1 < m_local_dim[1]; ++local_k1) {
			int k1 = local_k1 + m_self_coord[1] * m_local_dim[1];
			double ksq1 = m_pk_ksq[k1];
			
			for (int local_k2 = 0; local_k2 < m_local_dim[2]; ++local_k2) {
				int k2 = local_k2 + m_self_coord[2] * m_local_dim[2];
				double ksq2 = m_pk_ksq[k2];
				
				double k = sqrt(ksq0 + ksq1 + ksq2);
				
				double a = real(rho1[index]);
				double b = imag(rho1[index]);
				double c = real(rho2[index]);
				double d = imag(rho2[index]);

	  			// manifestly symmetry definition of cross power
	  	 		double pk = sqrt(pow(a*c-b*d,2) + pow(a*d+b*c,2));

	  			// 3-D anti-CIC filter for deconvolution
				// pk *= m_pk_cic[k0]*m_pk_cic[k1]*m_pk_cic[k2];

	  			// correct for dimensionless fft volume
	  			pk *= inversevolume;

	  			int tk0 = k0*(k0 <= ng/2) + (k0-ng)*(k0 > ng/2);
	  			int tk1 = k1*(k1 <= ng/2) + (k1-ng)*(k1 > ng/2);
	  			int tk2 = k2*(k2 <= ng/2) + (k2-ng)*(k2 > ng/2);
	  			double ii[3] = { 1.0*tk0, 1.0*tk1, 1.0*tk2 };
	  			double ksqr = (ii[0]*ii[0] + ii[1]*ii[1] + ii[2]*ii[2]);
	  			double legendre2 = 0.0;
	  			if(los >= 0 && los < 3){
	    			legendre2 = 2.5*(3.0*ii[los]*ii[los]/ksqr - 1.0);
				}
	  			if(!(k0==0 && k1==0 && k2==0)){
	    			pkbins.addVal(k, pk, pk*legendre2);
				}
	  			if( pkmodes.nmode3d > 0){
	    			pkmodes.addVal(tk0, tk1, tk2, pk, legendre2);
				}
	  			
				index++;
			}
			//index += d.m_d.padding[2];
      	}
      	//index += d.m_d.padding[1];
    }

    pkbins.accumulate();
    // remember to correct weight by half for real symmetry after here

    if (pkmodes.nmode3d > 0){
      pkmodes.accumulate();
  	}
  }


protected:

  hacc::Dfft &m_dfft;
  hacc::Distribution &d;
  int ng;

  const int *m_local_dim;
  const int *m_self_coord;

  std::vector<double> m_pk_cic;
  std::vector<double> m_pk_weight;
  std::vector<int>    m_pk_ksq;
};

#endif
