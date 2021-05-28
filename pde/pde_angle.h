#ifndef pde_angle_h
#define pde_angle_h
#include <vector>
#include <array>
#include "parameters/parameters.h"
#include "main/tags.h"
#include "main/traits.h"
#include "pde/pde_utilities.h"


template < class TFilament, class tag_BC>
class PDE_base_angle
{
public:
	typedef tag_non_constant_lhs tag_lhs;
	typedef tag_dont_expand_components tag_component_expansion;
	
	typedef double Unknown_type;
	typedef TFilament Filament_type;
	typedef LHS_arguments <TFilament,tag_component_expansion> LHS_arg_type;
	typedef Boundary_condition_base <Unknown_type,Filament_type,tag_component_expansion> BC_base_type;
	
	Parameters * p;
	Diff lhs_diff;
	double c4;
	Block <TFilament,tag_component_expansion> block;
	std::array <BC_base_type *, 4> bc;
	std::array <unsigned int, 4> bc_id; // lines of the bulk pde that'll be replaced with a bc
	
	PDE_base_angle(Parameters & par) : p(&par) //, bc_id ({{0,1,par.N-2,par.N-1}})
		{
		// The swimming head BC is equivalent to psi(0)=0 in the infinite rotational friction limit,
		// but it is enforced differently, by cancelling the right hand side of dpsi/dt=...
		// This only works, however, if the bulk equation of motion for s=0 is present.
		// As a result, we keep it and put the boundary condition on line 2 instead.
		bc_id = std::array <unsigned int, 4> ({{2,1,par.N-2,par.N-1}});
//		bc_id = std::array <unsigned int, 4> ({{0,1,par.N-2,par.N-1}});
		if ( std::is_same<tag_BC,tag_swimming_head>::value ) bc_id[0] = 2;
		c4 = par.dt * par.sp_inv;
		lhs_diff  = 1.5 * D0 + c4 * D4;
		
		bc = {{ new typename BC_type<BC_base_type,tag_BC>::angle1 (par), 
		        new typename BC_type<BC_base_type,tag_BC>::angle2 (par), 
		        new BCa_force_free_tail <BC_base_type> (par),
		        new BCa_torque_free_tail <BC_base_type> (par) }};
		};
	
// PDE: dpsi/dt + c4*d4psi/dl4 = d2f/dl2 - gamma*f*dpsi/dl 
//	                          + tension*d2psi/dl2 + (gamma+1)*d(tension)/dl*dpsi/dl
// =======================================================================
	
	// left-hand side
	double lhs_bulk (LHS_arg_type args) const
		{
		return 1.5 * D0(args) + c4 * D4(args);
		};
	
	// right-hand side
	Unknown_type rhs (const unsigned int & i, const Filament_type & fil, const double & t) const
		{
		return ( 
			    fil.du[i] * 
					 (
					 ( p->gamma + 1 ) * D1(fil.tension,i)
					 + p->gamma * fil.du[i] * fil.d2u[i] 
					 - p->gamma * fil.f[i] 
					 )
			    + fil.tension[i] * fil.d2u[i]
			    ) * p->sp_inv;
		};
	
};


#endif
