#ifndef params_h
#define params_h
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <typeinfo>
#include "configFile.h"
#include "geomvec/geomvec2.h"
#include "geomvec/geomvec3.h"
#include "main/util.h"


class Parameter {
public:
	string name;
	string type;
	void * value;
	
	Parameter(string name, string type) : name(name), type(type) {};
//	~Parameter() { delete_value(); };
//	void delete_value();
};

// specialization of to_string that checks the type of the parameter first
// general definition is in 'util.h'
template <> inline std::string to_string(const Parameter & p) {
	std::stringstream ss;
	if      (p.type==typeid(bool).name())         ss << *(static_cast <bool*>         (p.value));
	else if (p.type==typeid(int).name())          ss << *(static_cast <int*>          (p.value));
	else if (p.type==typeid(unsigned int).name()) ss << *(static_cast <unsigned int*> (p.value));
	else if (p.type==typeid(double).name())       ss << *(static_cast <double*>       (p.value));
	else if (p.type==typeid(string).name())       ss << *(static_cast <string*>       (p.value));
	else if (p.type==typeid(Geomvec3).name()) (static_cast <Geomvec3*> (p.value))->write(ss,"(",",",")");
	else if (p.type==typeid(Geomvec2).name()) (static_cast <Geomvec2*> (p.value))->write(ss,"(",",",")");
	else ss << "unknown type";
	return ss.str();
};



class Parameters {
public:
	ConfigFile conf;
	std::vector <Parameter> par_list;
	
	// general properties //
	bool verbose;
	std::string outdir;           // output directory
	std::string posdir;
	std::string checkpoint_dir;
	bool use_checkpoint;
	string gnuplot_preamble_file;
	
	// time and space discretization
	double L;                     // filament length
	double dl;                    // spatial step
	double dl2, dl3, dl4;
	double dl_inv, dl2_inv, dl3_inv, dl4_inv;
	double dl_print;
	unsigned int nl_print;
	unsigned int N;
	double dt;                    // time step
	double t_sim;                 // duration of simulation
	double t_start;               // when to start recording
	double t_print;               // printing period
	double t_pos;                 // position printing period
	double t_plot;                // plotting period
	double t_sleep;               // wait for t_sleep after each iteration
	unsigned int report_freq;     // how many times progress and time are reported
	unsigned int n_step, n_start, n_print, n_pos, n_plot, n_report;
	double t_report;
	
	// filament properties //
	std::string description_type; // type of description (angle or position)
	std::string BC_type;          // type of boundary condition at the head
	double sp;                    // sperm number
	double sp4, sp_inv;
	double gamma;                 // ratio of frictions along and perp. to filament
	double A0;                    // amplitude of initial condition
	double k0, k0_;               // wave vector of initial condition
	
	// internal force properties //
	double k, k_;                 // wave vector of sinusoidal internal force
	double w, w_;                 // pusation of sinusoidal internal force
	                              // (k and w are frequencies, k_ and w_ are angular frequencies)
	double A;                     // amplitude of sinusoidal internal force
	
	// head properties //
	double head_drag, head_rot_drag;
	
	//----------------------------
	
	Parameters(char* filename);
	void print(std::ostream & os);
	void print(const std::string & filename);
	
	template <class T> void read(T & t, std::string name);
	template <class T> void read(T & t, std::string name, T default_value);
	template <class T> void set(T & t, std::string name);  // adds an entry in par_list and bind it to 't'
	template <class T> void reset(T & t, std::string name); // reset the value of a parameter already in par_list
};


template <class T> void Parameters::read(T & t, std::string name)
	{
	conf.readInto(t,name);
	this->set(t,name);
	};

template <class T> void Parameters::read(T & t, std::string name, T default_value)
	{
	conf.readInto(t,name,default_value);
	this->set(t,name);
	};

template <class T> void Parameters::set(T & t, std::string name)
	{
	Parameter p(name,typeid(T).name());
	// 'value' can be set to '&t' if 't' is going to stick around, e.g. if it's a member of 'Parameters'.
	// In this case 't' needs to be passed by reference, and to be assigned a value somewhere else.
	// As a result, a call like 'set(true,"verbose")' is not possible.
	// To create parameters at run-time, the assignment should be 'value = new T(t)'.
	// In this case, a 'delete' function is needed for 'value', that first determines the type.
	// I haven't managed to write that yet.
	// This method has two advantages: 't' can be passed by value, so that 'set(true,"verbose")' works,
	// and still works if 'verbose' was never declared.
	p.value = &t;
	par_list.push_back(p);
	};

// This is useless: it's equivalent to just 'par.name = t' !
// It would become useful with run-time parameter creation though (see 'Parameters::set').
template <class T> void Parameters::reset(T & t, std::string name)
	{
	Parameter p(name,typeid(T).name());
	std::vector <Parameter>::iterator I = par_list.begin();
	while ( I!=par_list.end() && I->name!=name ) I++;
	if (I==par_list.end()) std::cerr << "[Parameter::reset] Unable to find parameter \"" << name << "\'" << std::endl;
	else p.value = &t;
	};


#endif
