#ifndef filament_h
#define filament_h
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "main/tags.h"
//#include "pde/solver.h"
//#include "pde/solver_boost.h"

#include "main/util.h"
#if SOLVER_TYPE==TRILINOS
	#include "pde/solver.h"
#elif SOLVER_TYPE==BOOST
	#include "pde/solver_boost.h"
#endif


using std::istream;
using std::ostream;

class Internal_force
{
public:
	double operator () (const int & i, const double & t, const Parameters & p)
		{ return p.A; 	};
};


template <class TPoint, class tag_unknown> class Filament;

template <class TPoint> class Filament <TPoint,tag_angle>
{
public:
	typedef TPoint Point_type;
	typedef double Unknown_type;
	
	const unsigned int & N;
	Parameters * p;
	std::vector <Point_type> position;
	std::vector <Unknown_type> u;
	std::vector <Unknown_type> du, d2u, d3u, d4u;  // spatial derivatives of position (noted u for shortness)
	std::vector <double> tension;
	std::vector <double> f, df;              // internal force and its spatial derivative
	typedef std::pair <std::string, std::vector <Unknown_type> *> Item_type;
	typedef std::list <Item_type> Full_state_type;
	Full_state_type full_state;
	Internal_force internal_force;
	
	Filament (Parameters & par);
	void get_derivatives (const double & t);
	void move_head ();
	void get_position ();
	double get_curvature_energy ();
};


template <class TPoint, class tag_unknown>
ostream & operator << (ostream & os, Filament <TPoint,tag_unknown> & filament);

template <class TPoint, class tag_unknown>
istream & operator >> (istream & is, Filament <TPoint,tag_unknown> & filament);


#include "filament.hpp"

#endif
