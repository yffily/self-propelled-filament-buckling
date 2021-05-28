#ifndef sytem_h
#define sytem_h
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "system/filament.h"
#include "pde/pde.h"
#include "main/tags.h"

template <class TFilament, class tag_BC, class tag_unknown> class PDE_type;

template <class TFilament, class tag_BC> class PDE_type <TFilament,tag_BC,tag_angle>
{ public: 
	typedef PDE < PDE_base_angle <TFilament,tag_BC> > position;
	typedef PDE < PDE_base_tension_camalet <TFilament,tag_BC> > tension;
};

//============================================================================================

template <class TPoint, class tag_BC, class tag_unknown> class Integrator
{
public:
	typedef Filament <TPoint,tag_unknown> Filament_type;
	typedef typename Filament_type::Unknown_type Unknown_type;
	typedef typename PDE_type<Filament_type,tag_BC,tag_unknown>::position PDE_p;
	typedef typename PDE_type<Filament_type,tag_BC,tag_unknown>::tension PDE_t;
	
	const unsigned int & N;
	PDE_p pde_position;
	PDE_t pde_tension;
	Solver <PDE_p> solver_position;
	Solver <PDE_t> solver_tension;
	std::vector <Unknown_type> rhs;
	std::vector <Unknown_type> rhs_last;
	std::vector <Unknown_type> rhs_next;
	std::vector <Unknown_type> u_last;
	double d2_max;
	typedef std::pair <std::string, std::vector <Unknown_type> *> Item_type;
	typedef std::list <Item_type> Full_state_type;
	Full_state_type full_state;
	
	Integrator (Parameters & par, Filament_type & filament, const double & t);
	void initialize_rhs (Filament_type & filament, const double & t);
	void get_tension (Filament_type & filament, const double & t);
	void step_forward (Filament_type & filament, const double & t);
	double get_max_displacement2 (Filament_type & filament);
};

//============================================================================================

template <class TPoint, class tag_BC, class tag_unknown> class System
{
public:
	typedef Filament <TPoint,tag_unknown> Filament_type;
	typedef Integrator <TPoint,tag_BC,tag_unknown> Integrator_type;
	
	long int step;
	double t;
	Filament_type filament;
	Integrator_type integrator;
	
	System (Parameters & par);
	void evolve ();
	void print_full (ostream & os) const;
	void read_full (istream & is);
};

#include "system.hpp"

#endif
