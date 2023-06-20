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

//void calculateProfile(globalArgs &params, std::vector<float> &rbins, std::vector<float> &mbins);

void calculateProfile(globalArgs &params, std::vector<float> &radiusBins, std::vector<float> &countBins,
		std::vector<float> &rhoBins){
	
	hacc::Distribution d(MPI_COMM_WORLD, params.ng);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int ng[3];
	int self[3];

	for(int dim=0; dim<DIMENSION; dim++){
		ng[dim] = d.local_ng_3d(dim);
		self[dim] = d.self_3d(dim);
	}
	int nbins = params.nbins;
	float maxRadius = params.ng/2;
	float minRadius = 1.0;
	float deltaRadius = (maxRadius - minRadius)/(nbins-1);

	radiusBins.at(0)=minRadius;
	for(int bin=1; bin<nbins; bin++){
		radiusBins.at(bin)=deltaRadius*bin+minRadius;
	}

	for(int bin=0; bin<nbins; bin++){
		countBins.at(bin) = 0;
		rhoBins.at(bin) = 0.0;	
	}

	int centerIndex[DIMENSION];
	for(int i=0; i<DIMENSION;i++){
		centerIndex[i]=params.ng/2;
	}


	float maxR = 0;	
	float minR = 128;
	float zero_count=0;
	float one_count=0;
	float total_count=0;
	long location[DIMENSION];
	long index=0;
	
	for(long i=0; i<params.mag.size();i++){

		location[0]=(int)params.xx.at(i);
		location[1]=(int)params.yy.at(i);
		location[2]=(int)params.zz.at(i);

		float diff[DIMENSION];
		for(int dim=0; dim<DIMENSION; dim++){
			diff[dim] = location[dim] - centerIndex[dim];
		}

		float r = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
		
		if(r>maxR){maxR=r;}
		if(r<minR){minR=r;}
		
		if(r<maxRadius){
			int bin;
			if(r>minRadius){
				bin = (int) (floor((r-minRadius)/deltaRadius)) + 1;
			} else {
				bin=0;
			}
			if(bin>nbins){
				bin=nbins-1;
			}
			countBins.at(bin)++;
			rhoBins.at(bin) += params.mag.at(index);
			total_count++;	
		}
		index++;
	}

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
    while((opt = getopt(argc, argv, "f:b:d:")) != -1){
		switch(opt){
			case 'f':
			{
				filename = optarg;
				dirflag=false;	
			}
				continue;
			case 'b':
			{
				std::stringstream s_nbins(optarg);
				s_nbins >> nbins;
			}
				continue;
			case 'd':
			{
				dirname= optarg;
				dirflag=true;
			}
				continue;
			case -1:
				break;
		}		
	}  

	std::vector<std::string> files;

	std::cout<<dirname<<std::endl;

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

	std::cout<<"nbins="<<nbins<<std::endl;

	for(int file=0; file<files.size();file++){
		filename=files.at(file);
		std::cout<<"Computing profile for "<<filename<<std::endl;
		globalArgs params;
		params.nbins=nbins;
		params.ng=0;
		
		long total_Ng; 
		long Ng;
		std::vector<double> mag;
		std::vector<GRID_T> xx;
		std::vector<GRID_T> yy;
		std::vector<GRID_T> zz;
			
		{ //scope GIO
		gio::GenericIO GIO(comm, filename);
		GIO.openAndReadHeader(gio::GenericIO::MismatchRedistribute);

		total_Ng = GIO.readTotalNumElems();
		Ng = GIO.readNumElems();

		mag.resize(Ng);
		xx.resize(Ng);
		yy.resize(Ng);
		zz.resize(Ng);

		GIO.addVariable("mag", mag.data(), gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("x", xx.data(), gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("y", yy.data(), gio::GenericIO::VarHasExtraSpace);
		GIO.addVariable("z", zz.data(), gio::GenericIO::VarHasExtraSpace);

		GIO.readData();
		GIO.close();
		}

		params.ng=(int)std::cbrt(total_Ng); //assume cubic box
		params.mag=mag;
		params.xx=xx;
		params.yy=yy;
		params.zz=zz;
		params.rank=rank;

		std::vector<float> radiusBins;
		std::vector<float> rhoBins;
		std::vector<float> countBins;

		radiusBins.resize(nbins);
		rhoBins.resize(nbins);
		countBins.resize(nbins);
		for(int i=0; i<params.mag.size(); i++){
			params.mag.at(i) = pow(params.mag.at(i), 2.0);
		}
		calculateProfile(params, radiusBins, countBins, rhoBins);

		MPI_Barrier(MPI_COMM_WORLD);

		std::ostringstream ost;
		//format outfile name
		size_t pos1=0;
		size_t pos2=0;
		pos1=filename.find_first_of("/");
		pos2=filename.find_last_of(".");
		ost<<filename.substr(0,pos1);
		ost<<"/profiles";
		ost<<filename.substr(pos1,(pos2-pos1));
		ost << ".profile";
		
		{
		gio::GenericIO GIO(MPI_COMM_WORLD, ost.str());

		GIO.setNumElems(radiusBins.size());

		GIO.addVariable("r", radiusBins.data());	
		GIO.addVariable("count", countBins.data());
		GIO.addVariable("rho", rhoBins.data());

		GIO.write();
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
}
