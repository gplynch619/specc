#include "sc.h"
#include "Distribution.hpp"
#include "power.h"

#define PK_SUFFIX "pk"
#define PK3D_SUFFIX "3dpk"

void writePk(SolverBase *solver, Parameters & params, std::string outBase, 
		std::string suffix, int los){
	hacc::Distribution & d = solver->get_d();
	int ng = d.m_d.n[0];

	double kmin1d = 1.0;
	double kmax1d = sqrt(3.0*(ng/2.0)*(ng/2.0));
	int nbin1d = params.pknbin1d();
	if(nbin1d == 0)
	  nbin1d = (int)ceilf(kmax1d - kmin1d);
	bool logbins = params.pklogbins();
	Bins1D pkbins(d.cart_1d(), nbin1d, kmin1d, kmax1d, logbins);

	int nmode3d = params.pknmode3d();
	Modes3D pkmodes(d.cart_1d(), nmode3d);

	solver->power_spectrum(pkbins, pkmodes, los);
	MPI_Barrier(MPI_COMM_WORLD);

	int rank = d.self(); 
	if(rank==0) {
		float rL = params.rL(); 
		float ng = 1.0*d.global_ng(0);
		float pi = 4.0*atan(1.0);
		float kcoeff = 2.0*pi/rL;
		float pkcoeff = powf(rL/ng,3.0);

	  std::string outName1d = outBase + "." + PK_SUFFIX;
	  
	  if(suffix.length() > 0)
		outName1d += "." + suffix;

	  FILE *outFile1d = fopen(outName1d.c_str(),"w");
	  fprintf(outFile1d,"# [<k> (h/Mpc)]\t[<P_0(k)> (Mpc/h)^3]\t[ErrorBar (Mpc/h)^3]\tnModes\t[<P_2(k)> (Mpc/h)^3]\n");

	  for(int i=0; i<pkbins.xBins.size(); i++) {

		double errorbar = pkcoeff*pkbins.monoBins[i];
		if(pkbins.wBins[i] > 0.0){
	 		 errorbar *= 1.0/sqrt(0.5*pkbins.wBins[i]);
		}

		fprintf(outFile1d, "%e\t%e\t%e\t%f\t%e\n", 
			kcoeff*pkbins.xBins[i], 
			pkcoeff*pkbins.monoBins[i],
			errorbar,
			0.5*pkbins.wBins[i],
			pkcoeff*pkbins.quadBins[i]);
	  }
	  fclose(outFile1d);

	  if(nmode3d > 0) {
		std::string outName3d = outBase + "." + PK3D_SUFFIX;
		if(suffix.length() > 0)
	  outName3d += "." + suffix;
		FILE *outFile3d = fopen(outName3d.c_str(),"w");
		fprintf(outFile3d, "# i\tj\tk\t[P(i,j,k) (Mpc/h)^3]\tL2(i,j,k)\n");
	  
		for(int i=-1*nmode3d+1; i<=nmode3d; i++)
	  for(int j=-1*nmode3d+1; j<=nmode3d; j++)
		for(int k=-1*nmode3d+1; k<=nmode3d; k++) {
		  int indx = (k+nmode3d-1) + (2*nmode3d)*((j+nmode3d-1) + (2*nmode3d)*(i+nmode3d-1));
		  fprintf(outFile3d, 
			"%d\t%d\t%d\t%e\t%f\n", 
			i, j, k, 
			pkcoeff*pkmodes.val3d[indx],
			pkmodes.legendre3d[indx]);
		}
		fclose(outFile3d);
    }
  }
}

template <>
void writeVector(std::vector<complex_t> &vector, Wavefunction &wf, std::string outFileBase){
		
	std::ostringstream ost;
	ost<<outFileBase<<".gio";

	int Ng = wf.Ng();
	gio::GenericIO GIO(wf.d().parent_comm(), ost.str());
	
	std::vector<double> mag;
	std::vector<double> arg;

	for(long i=0; i<Ng; ++i){
		mag.push_back(std::abs(vector[i]));
		arg.push_back(std::arg(vector[i]));
	}

	unsigned CoordFlagsX = gio::GenericIO::VarIsPhysCoordX;
	unsigned CoordFlagsY = gio::GenericIO::VarIsPhysCoordY;
	unsigned CoordFlagsZ = gio::GenericIO::VarIsPhysCoordZ;
	

	GIO.setNumElems(Ng);
	GIO.setPhysOrigin(0.0);
	
	wf.m_xx.resize(Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	wf.m_yy.resize(Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	wf.m_zz.resize(Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	mag.resize(Ng+GIO.requestedExtraSpace()/sizeof(double));
	arg.resize(Ng+GIO.requestedExtraSpace()/sizeof(double));

	GIO.addVariable("x", wf.m_xx.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("y", wf.m_yy.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("z", wf.m_zz.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("mag", mag.data(),  gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("arg", arg.data(),  gio::GenericIO::VarHasExtraSpace);

	GIO.write();

	wf.m_xx.resize(Ng);
	wf.m_yy.resize(Ng);
	wf.m_zz.resize(Ng);
	std::cout<<"test 3"<<std::endl;
}
