#ifndef INIT_H
#define INIT_H

#include "sctypes.h"
#include "Distribution.hpp"
#include "CosmoClass.h"

#include "Wavefunction.h"

int initialize_wavefunction(int32_t NG_Max, CosmoClass &CC, TransferClass &TC, double rL, 
		double z_in, unsigned long seed, hacc::Distribution &d, Wavefunction &wf);

void initialize_powerspec(int32_t NG_Max, CosmoClass &CC, TransferClass &TC, double rL, 
		double z_in, unsigned long seed, hacc::Distribution &d, Wavefunction &wf);

void initialize_planewave(int* mode, double rL, unsigned long seed, 
		hacc::Distribution &d, Wavefunction &wf);

#endif
