#include "initread.h"
#include <iostream>

int main(int argc, char *argv[]){
	MPI::Init();

	int rank;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int argvi = optind;

	std::string pName;
	pName = argv[argvi];

	Parameters params(pName);

	std::string paramout = params.getParamsString();

	if(rank == 0){
		std::cout << paramout << std::endl;
	}

	//std::cout<<"End"<<std::endl;

	return 1;
}
