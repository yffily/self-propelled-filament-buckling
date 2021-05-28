#include "parameters.h"
using std::endl;
using std::max;


//__________________//
//                  //
// Parameters class //
//__________________//

Parameters::Parameters(char* filename) : conf(filename)
	{
	// general properties //
	read(verbose, "verbose");
	read(outdir, "outdir");
	read(use_checkpoint, "use_checkpoint");
	read(gnuplot_preamble_file, "gnuplot_preamble_file");

	// time and space discretization
	read(L, "L");
	read(dl, "dl");
	read(dl_print, "dl_print");
	read(dt, "dt");
	read(t_sim, "t_sim");
	read(t_start, "t_start");
	read(t_print, "t_print");
	read(t_pos, "t_pos");
	read(t_plot, "t_plot");
	read(t_sleep, "t_sleep");
	read(report_freq, "report_freq");
	
	// filament properties //
	read(description_type, "description_type");
	read(BC_type, "BC_type");
	read(sp, "sp");
	read(gamma, "gamma");
	read(A0, "A0");
	read(k0, "k0");
	
	// internal force properties //
	read(k, "k");
	read(w, "w");
	read(A, "A");
	
	// head properties //
	read(head_drag, "head_drag");
	read(head_rot_drag, "head_rot_drag");
	
	//----------------------------
	
//	// print what's been read before starting to mess around with the parameters
//	// !!! outdir may not exist yet!
//	print(outdir+"/parameters_initial.dat");
	
	//----------------------------
	
	if ( outdir == "auto" )
		{
		outdir = "data/bc_"+BC_type+"_sp"+to_string(sp)+"_gamma"+to_string(gamma)+"_A"+to_string(A)+"_k"+to_string(k);
//		if ( BC_type == "swimming" ) outdir += "_head-drag"+to_string(head_drag)+"_head-size"+to_string(head_size);
		if ( BC_type == "swimming" ) outdir += "_head-drag"+to_string(head_drag)+"_head-rot-drag"+to_string(head_rot_drag);
		}
	
	posdir = outdir+"/position";
	set(posdir,"posdir");
	checkpoint_dir = outdir+"/checkpoint";
	set(checkpoint_dir,"checkpoint_dir");
	
	N        = (unsigned int) (L/dl+1.5);
	dl       = L/(N-1);
	nl_print = max( 1, int(dl_print/dl + 0.5) );
	dl_print = dl * nl_print;
	n_step   = t_sim/dt + 0.5;
	t_sim    = dt * n_step;
	n_start  = (t_start<=0) ? 0 : t_start/dt + 0.5;
	n_print  = (t_print<=0) ? 0 : max( 1 , int(t_print/dt + 0.5) );
	n_pos    = t_pos/dt + 0.5;
	n_plot   = t_plot/dt + 0.5;
	n_report = max( n_step/report_freq , 1u );
	t_report = n_report * dt;
	set(N, "N");
	set(n_step, "n_step");
	set(n_start, "n_start");
	set(n_print, "n_print");
	set(n_pos, "n_pos");
	set(n_plot, "n_plot");
	set(n_report, "n_report");
	set(t_report, "t_report");
	
	sp4 = sp * sp * sp * sp;
	sp_inv = 1./sp;
	dl2 = dl  * dl;
	dl3 = dl2 * dl;
	dl4 = dl2 * dl2;
	dl_inv  = 1./dl;
	dl2_inv = 1./dl2;
	dl3_inv = 1./dl3;
	dl4_inv = 1./dl4;
	
	k_ = k*dblpi;
	w_ = w*dblpi;
	k0_ = k0*dblpi;
	
//	head_rot_drag = head_drag * pow(head_size,2);
//	set(head_rot_drag, "head_rot_drag");
}


void Parameters::print(std::ostream & os) {
	int w = 0;
	for (auto I=par_list.begin(); I!=par_list.end(); I++) w = max(w,int(I->name.size()));
	w += 1;
	os << std::setprecision(15);
	for (auto I=par_list.begin(); I!=par_list.end(); I++)
		{ os << std::setw(w) << I->name << " = " << to_string(*I) << endl; }
}

void Parameters::print(const std::string & filename)
	{
	std::ofstream os(filename.c_str());
	print(os);
	os.close();
	}

