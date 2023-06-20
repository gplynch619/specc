#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <vector>

#include <mpi.h>

#define NAMEWIDTH 5

using std::string;
using std::vector;

class Timer{

	public:
		typedef int TimerRef;

		Timer();
		~Timer();

		static TimerRef getTimer(const char *nm);
		static void startTimer(TimerRef tr);
		static void stopTimer(TimerRef tr);
		static double getCurrent(TimerRef tr);
		static void timerStats(TimerRef, bool uselaccum=false, int nameWidth=NAMEWIDTH);
		static void alltimerStats(TimerRef, bool uselaccum=false, int nameWidth=NAMEWIDTH);
		static void stopTimerStats(TimerRef, bool uselaccum=false, int nameWidth=NAMEWIDTH);

		static void timingStats(double t, string name, int nameWidth=NAMEWIDTH);
		static void alltimingStats(double t, string name, int nameWidth=NAMEWIDTH);
	
	private:
		static vector<string> names;
		static vector<double> start;
		static vector<double> stop;
		static vector<double> accum;
		static vector<double> laccum;
};

#endif
