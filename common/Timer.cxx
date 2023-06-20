#include "Timer.h"

#include <cstdio>

using namespace std;

vector<string> Timer::names;
vector<double> Timer::start; 
vector<double> Timer::stop;
vector<double> Timer::accum; 
vector<double> Timer::laccum;

Timer::Timer() {}
Timer::~Timer() {}

Timer::TimerRef Timer::getTimer(const char *nm){
	string s(nm);
	TimerRef tr;

	for(tr=0; tr < static_cast<TimerRef>(names.size()); tr++){
		if(names[tr].compare(s) == 0){ break; }
	}
	if(names.size()==0 || tr == static_cast<TimerRef>(names.size()) ) {
		names.push_back(s);
		start.push_back(0.0);
		stop.push_back(0.0);
		accum.push_back(0.0);
		laccum.push_back(0.0);
	}
	return tr;
}

void Timer::startTimer(TimerRef tr){
	start[tr] = MPI_Wtime();
}

void Timer::stopTimer(TimerRef tr){
	stop[tr] = MPI_Wtime();
	accum[tr] += stop[tr]-start[tr];
	laccum[tr] += stop[tr]-start[tr];
}

double Timer::getCurrent(TimerRef tr){
	return MPI_Wtime() - start[tr];
}

void Timer::timingStats(double t, string name, int nameWidth){
	double tmax, tmin, tavg;
	int rank, nr;
	double numranks;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nr); 
	numranks=1.0*nr;

	MPI_Allreduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&t, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&t, &tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tavg/=numranks;

	if(rank==0){
		string tmpname(name.c_str());
		tmpname.resize(nameWidth, ' ');
		printf("%s  max  %.3e s  avg  %.3e s  min  %.3e s\n",
			tmpname.c_str(), tmax, tavg, tmin);
		fflush(stdout);

		fprintf(stderr, "%s done\n", tmpname.c_str());
		fflush(stderr);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	return;
}

void Timer::alltimingStats(double t, string name, int nameWidth){
	double tmax, tmin, tavg;
	int rank, nr;
	double numranks;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nr);
	numranks=1.0*nr;

	MPI_Allreduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&t, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&t, &tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	tavg/=numranks;

	double times[nr];
	if(rank==0){
		string tmpname(name.c_str());
		tmpname.resize(nameWidth, ' ');
		printf("%s  max  %.3e s  avg  %.3e s  min  %.3e s\n",
			tmpname.c_str(), tmax, tavg, tmin);
		fflush(stdout);

		fprintf(stderr, "%s done\n", tmpname.c_str());
		fflush(stderr);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	printf("%d, %10.5f\n",rank,t);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	return;
}

void Timer::timerStats(TimerRef tr, bool uselaccum, int nameWidth) {
	double t;
	t = uselaccum ? laccum[tr] : stop[tr] - start[tr];
	timingStats(t, names[tr], nameWidth);
	laccum[tr] = 0.0;
	return;
}

void Timer::stopTimerStats(TimerRef tr, bool uselaccum, int nameWidth){
	stopTimer(tr);
	timerStats(tr, uselaccum, nameWidth);
	return;
}

void Timer::alltimerStats(TimerRef tr, bool uselaccum, int nameWidth) {
	double t;
	t = uselaccum ? laccum[tr] : stop[tr] - start[tr];
	alltimingStats(t, names[tr], nameWidth);
	laccum[tr]=0.0;
	return;
}
