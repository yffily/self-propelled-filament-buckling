#ifndef pde_camalet_tension_h
#define pde_camalet_tension_h
#include <vector>
#include <array>
#include <cmath>
#include "parameters/parameters.h"
#include "main/tags.h"
#include "main/traits.h"
#include "pde/pde_utilities.h"


template <class TFilament, class tag_BC>
class PDE_base_tension_camalet
{
public:
	typedef tag_non_constant_lhs tag_lhs;
	typedef tag_dont_expand_components tag_component_expansion;
	
	typedef double Unknown_type;
	typedef TFilament Filament_type;
	typedef LHS_arguments <TFilament,tag_component_expansion> LHS_arg_type;
	typedef Boundary_condition_base <Unknown_type,Filament_type,tag_component_expansion> BC_base_type;
	
	Parameters * p;
	Diff lhs_diff, lhs_diff_non_constant;
	std::array <BC_base_type *, 2> bc;
	std::array <unsigned int, 2> bc_id; // lines of the bulk pde that'll be replaced with a bc
	
	PDE_base_tension_camalet (Parameters & par) : p(&par), bc_id ({{0,par.N-1}})
		{
		lhs_diff  = par.gamma*D2; // + -1*D0;
		lhs_diff_non_constant = D0;
		
		// boundary conditions
		bc = {{ new typename BC_type<BC_base_type,tag_BC>::tension_camalet (par), 
		        new BCtc_force_free_tail <BC_base_type> (par) }};
		};
	
// ODE: gamma*d2T - (dpsi/dl)^2*T = gamma * d( dpsi/dl*f - dpsi/dl*d2psi/dl2 )/dl
//                                    + dpsi/dl * ( df/dl - d3psi/dl3 )
// =======================================================================================
	
	// left-hand side
	double lhs_bulk (LHS_arg_type args) const
		{
		
//		cout << p->gamma << sep << D2(args) << sep << args.filament.du[args.i1] << sep << D0(args) << endl;
		
		return p->gamma*D2(args) 
		       - norm2(args.filament.du[args.i1])*D0(args);
		};
	
	// right-hand side
	Unknown_type rhs (const unsigned int & i, const Filament_type & fil, const double & t) const
		{
		return - (p->gamma+1) * fil.du[i] * fil.d3u[i] 
		       - p->gamma * ( fil.d2u[i] * fil.d2u[i] + fil.df[i] );
		       
		};
	
};


#endif
