#ifndef SC_H
#define SC_H

#include <iostream>
#include <string>
#include <cassert>
#include <algorithm>
#include <vector>

#include "poisson.h"
#include "initread.h"
#include "GenericIO.h" 
#include "sctypes.h"
#include "Wavefunction.h"

void writePk(SolverBase *solver,
		Parameters & params,
		std::string outBase,
		std::string suffix,
		int los=-1);


template <class T> 
void writeVector(std::vector<T> &vector, Wavefunction &wf, std::string outFileBase){
		
	std::ostringstream ost;
	ost<<outFileBase<<".gio";

	gio::GenericIO GIO(wf.d().parent_comm(), ost.str());
	
	unsigned CoordFlagsX = gio::GenericIO::VarIsPhysCoordX;
	unsigned CoordFlagsY = gio::GenericIO::VarIsPhysCoordY;
	unsigned CoordFlagsZ = gio::GenericIO::VarIsPhysCoordZ;
	
	int Ng = wf.Ng();

	GIO.setNumElems(Ng);
	GIO.setPhysOrigin(0.0);
	
	wf.m_xx.resize(Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	wf.m_yy.resize(Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	wf.m_zz.resize(Ng+GIO.requestedExtraSpace()/sizeof(GRID_T));
	vector.resize(Ng+GIO.requestedExtraSpace()/sizeof(T));

	GIO.addVariable("x", wf.m_xx.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("y", wf.m_yy.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("z", wf.m_zz.data(), CoordFlagsX | gio::GenericIO::VarHasExtraSpace);
	GIO.addVariable("vec", vector.data(),  gio::GenericIO::VarHasExtraSpace);

	GIO.write();

	wf.m_xx.resize(Ng);
	wf.m_yy.resize(Ng);
	wf.m_zz.resize(Ng);
}
#endif
