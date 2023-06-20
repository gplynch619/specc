#include <string>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>

#include "GenericIO.h"
#include "sctypes.h"
#include "complex-type.h"
#include "Distribution.hpp"

#define DIMENSION 3
#define EPS 1e-10
#define EPS1 0.0

struct globalArgs {
	std::string filename;
	int nbins;
	int ng;
	int rank;
	std::vector<double> mag;
	std::vector<GRID_T> xx;
	std::vector<GRID_T> yy;
	std::vector<GRID_T> zz;
};

void readFile(MPI_Comm comm, std::string filename, std::vector<double> &outarr){
	
		long int Ng;

		gio::GenericIO GIO(comm, filename);
		GIO.openAndReadHeader(gio::GenericIO::MismatchRedistribute);

		Ng = GIO.readNumElems();

		outarr.resize(Ng);

		GIO.addVariable("arg", outarr.data(), gio::GenericIO::VarHasExtraSpace);

		GIO.readData();
		GIO.close();

}

int main(int argc, char *argv[]){

	MPI::Init();
	int rank;
	int nranks;
	MPI_Comm comm;

	comm=MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	std::string filename;
	std::string dirname;

    int opt;
	int nbins;

	bool dirflag;
    // put ':' in the starting of the 
    // string so that program can  
    //distinguish between '?' and ':'  

	int endstep;	
	int argvi = optind;

	endstep = atoi(argv[argvi]);

/*
	while((opt = getopt(argc, argv, "f:b:d:")) != -1){
		switch(opt){
			case 'n':
			{
				endstep = atoi(optarg);
			}
				continue;
			case -1:
				break;
		}		
	}  
	if(dirflag){
		DIR *dir;
		struct dirent *ent;
		if((dir = opendir(dirname.c_str()))!=NULL){
			while((ent=readdir(dir))!=NULL){
				std::stringstream ost;
				ost<<dirname;
				printf ("%s\n", ent->d_name);
				ost<<ent->d_name;
				size_t pos=0;
				pos=ost.str().find_last_of(".");
				if(ost.str().substr(pos+1)=="gio"){ files.push_back(ost.str());	}
			}
			closedir(dir);
		} else {
			perror("COULD NOT OPEN");
			return EXIT_FAILURE;
		}
	} else {
		files.push_back(filename);
	}
*/

	double averages[endstep+1];
	double maxes[endstep+1];
	double mins[endstep+1];

	long int Ng;
	
	std::cout<<"Here"<<std::endl;

	for(int step=1; step<endstep;step++){
		std::string filename1;
		std::string filename2;
		
		std::ostringstream ost;
		ost<<"intermediate_phase.step."<<step<<".gio";

		filename1=ost.str();
		ost.str("");
		ost<<"intermediate_phase.step."<<step+1<<".gio";
		filename2=ost.str();
		ost.str("");

		std::vector<double> phase1;
		std::vector<double> phase2;

		readFile(comm, filename1, phase1);
		MPI_Barrier(comm);
		readFile(comm, filename2, phase2);
		MPI_Barrier(comm);

		std::vector<double> ratio;

		if(phase1.size()==phase2.size()){
			std::cout<<"true"<<std::endl;
		}

		Ng = phase1.size();
	
		int zerocount = 0;

		for(int i=0; i<phase1.size(); i++){
			double num, den, result;
			if(phase2.at(i)==0){
				num=phase2.at(i)+EPS;
				zerocount++;
			} else {
				num=phase2.at(i);
			}
			
			if(phase1.at(i)==0){
				den=phase1.at(i)+EPS;
			} else {
				den=phase1.at(i);
			}
			result = num/(den+EPS1);
			ratio.push_back(result);
		}
		
		std::cout<<"Number of zeroes "<<zerocount<<std::endl;

		MPI_Barrier(MPI_COMM_WORLD);
	
		
		double min_psi, max_psi, ave_psi;

		min_psi=1.0e30;
		max_psi=-1.0e30;
		ave_psi=0.0;

		for(long i=0; i<Ng; ++i){
			ave_psi += ratio.at(i);
			if(ratio.at(i) > max_psi){ max_psi = ratio.at(i); }
			if(ratio.at(i) < min_psi){ min_psi = ratio.at(i); }
		}
		
		ave_psi/=Ng;

		MPI_Allreduce(MPI_IN_PLACE, &max_psi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &min_psi, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &ave_psi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	
		ave_psi/=nranks;

		averages[step-1] = ave_psi;
		mins[step-1] = min_psi;
		maxes[step-1] = max_psi;
	
	}	
	
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){

		std::string outName = "phase_values.txt";
	
		FILE *outFile = fopen(outName.c_str(),"a");
		for(int i =0; i<endstep;i++){

			fprintf(outFile, "%i\t%f\t%f\t%f\n", 
				i, 
				averages[i],	
				mins[i],
				maxes[i]);
			}  
		fclose(outFile);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
