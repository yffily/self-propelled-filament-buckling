#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <chrono>
#include <thread>
#include "main/util.h"
#include "io/plotter.h"
#include "system/system.h"
#include "geomvec/geomvec2.h"
#include "pde/diff.h"
#include "main/checkpoint.h"

using std::vector;
using std::ofstream;
using std::cout;
using std::cin;
using std::clog;
using std::cerr;
using std::endl;
// cout, cerr and endl are already defined in one of the trilinos headers
namespace bfs = boost::filesystem;
using boost::filesystem::path;
using boost::filesystem::create_directory;

template <class TPoint> void run (Parameters & par);
template <class TPoint, class tag_unknown> void run (Parameters & par);
template <class TPoint, class tag_BC, class tag_unknown> void run (Parameters & par);

//===============================================================

int main (int argc, char *argv[])
	{
	if ( argc != 2 ) { cerr << "[Error] Please provide a parameter file." << endl; exit(1); }
	Parameters par(argv[1]);
	run <Geomvec2> (par);
	return 0;
	}

//===============================================================

template <class TPoint> void run (Parameters & par)
	{
//	if ( par.description_type == "position" )   run <TPoint,tag_position> (par);
//	else 
	if ( par.description_type == "angle" ) run <TPoint,tag_angle>    (par);
	else
		{
		cerr << "[Error] '"+par.description_type+"' is not a known description type." << endl;
		exit(2);
		}
	}

template <class TPoint, class tag_unknown> void run (Parameters & par)
	{
	if ( par.BC_type == "pivoting" )      run <TPoint,tag_pivoting_head,tag_unknown> (par);
	else if ( par.BC_type == "free" )     run <TPoint,tag_free_head,tag_unknown>     (par);
	else if ( par.BC_type == "clamped" )  run <TPoint,tag_clamped_head,tag_unknown>  (par);
	else if ( par.BC_type == "swimming" ) run <TPoint,tag_swimming_head,tag_unknown> (par);
	else
		{
		cerr << "[Error] '"+par.BC_type+"' is not a known boundary condition type." << endl;
		exit(2);
		}
	}

template <class TPoint, class tag_BC, class tag_unknown> void run (Parameters & par)
	{
	path outdir(par.outdir);
	create_directory(par.outdir);
	create_directory(par.posdir);
	create_directory(par.checkpoint_dir);
	
	diff_initialize(par);
	System <TPoint,tag_BC,tag_unknown> sys(par);
	par.print(par.outdir+"/parameters_final.dat");
	ofstream os((par.posdir+"/pos_step0.dat").c_str());
	os << sys.filament;
	os.close();
	
	Checkpoint chkpt(sys,par);
	
#if PLOTTER_TYPE==GNUPLOT
	Plotter gp(par);
	if ( par.n_plot>0 )  gp.plot(sys.filament.position,"time = "+to_string(sys.t));
#else
	Plotter_dumb gp;
#endif
	
	int progress_report = sys.step/par.n_report;
	cout << progress_report << "/" << par.report_freq << sep << get_date() << endl;
	
	
	auto m = (chkpt.used) ? fstream::app : fstream::trunc;
	ofstream os_max_displacement((par.outdir+"/max_displacement.dat").c_str(),m);
	ofstream os_head((par.outdir+"/time-headPos-headAng-curvEnergy.dat").c_str(),m);
	os_head.precision(8);
	if (!chkpt.used) os_head << "# time, head position, head angle, total curvature energy" << endl;
	while (sys.step<=par.n_step)
		{
		sys.evolve();
		sys.step++;
		if ( sys.step>=par.n_start )
			{
			if ( !gp.is_dumb() && par.n_plot>0 && sys.step % par.n_plot == 0 )
				{
				sys.filament.get_position();
				gp.plot(sys.filament.position,"time = "+to_string(sys.t));
				}
			if ( par.n_pos>0 && sys.step % par.n_pos == 0 )
				{
				unsigned int q = sys.step/par.n_pos;
				ofstream os((par.posdir+"/pos_step"+to_string(q)+".dat").c_str());
				os << sys.filament;
				os.close();
				}
			if ( par.n_print>0 && sys.step % par.n_print == 0 )
				{
				os_max_displacement << sys.t << sep << sqrt(sys.integrator.d2_max) << endl;
				sys.integrator.d2_max = 0;
				os_head << sys.t << sep << sys.filament.position[0] << sep << sys.filament.u[0] << sep << sys.filament.get_curvature_energy() << endl;
				}
			}
		if ( par.n_report>0 && sys.step % par.n_report == 0 )
			{
			cout << ++progress_report << "/" << par.report_freq << sep << get_date() << endl;
			if (par.use_checkpoint) chkpt.save_checkpoint(sys);
			}
		if ( par.t_sleep>0 )
			{
			std::this_thread::sleep_for(std::chrono::milliseconds(uint(1000*par.t_sleep)));
			}
		}
	}

