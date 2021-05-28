#ifndef boundary_condition_h
#define boundary_condition_h
#include <array>
#include <cmath>
#include "main/traits.h"
#include "main/tags.h"
#include "pde/diff.h"
#include "pde/pde_utilities.h"

//using namespace Diff;


// =============================
// class Boundary_condition_base
// =============================

template <class TUnknown, class TFilament, class tag_component_expansion>
class Boundary_condition_base
{
public:
	typedef LHS_arguments <TFilament,tag_component_expansion> LHS_arg_type;
	typedef TUnknown Unknown_type;
	typedef TFilament Filament_type;
	
	Parameters & p;
	int i_ev;     // point at which the BC is evaluated
	Row row;
	Block <TFilament,tag_component_expansion> block;
	
	Boundary_condition_base(Parameters & p, int i_ev, Diff lhs_diff) : p(p), i_ev(i_ev), row(lhs_diff.get_row(i_ev)) {};
	int get_nnz () const   { return row.indices.size();   };
	int get_k_min () const { return row.indices.front(); };
	int get_k_max () const { return row.indices.back(); };
	virtual double lhs(const LHS_arg_type & args) = 0;
	virtual TUnknown rhs(const TFilament & filament, const double & t) = 0;
};

//___________________________________________________________________________


#define create_bc(NAME,I_EV,LHS_DIFF,LHS,RHS) \
template <class BC_base_type> class NAME : public BC_base_type \
	{ \
	public: \
	using BC_base_type::i_ev; using BC_base_type::p; using BC_base_type::block; \
	typedef typename BC_base_type::LHS_arg_type LHS_arg_type; \
	typedef typename BC_base_type::Unknown_type Unknown_type; \
	typedef typename BC_base_type::Filament_type Filament_type; \
	NAME (Parameters & p) : BC_base_type( p , I_EV , LHS_DIFF ) {}; \
	double lhs (const LHS_arg_type & args) { LHS }; \
	Unknown_type rhs(const Filament_type & fil, const double & t) { RHS }; \
	};

// Note: some of the BC below could be made more implicit 
// by moving some non-linear terms from the rhs to the lhs,
// but it would require storing more quantities at current
// and previous time to keep the time integration second order.

//--------------------------------------------

// Boundary conditions for the angle
// =================================

// Force free head: d2psi/dl2(0) = 0
create_bc( BCa_force_free_head , 0 , D2 ,
           return D2(i_ev,args.j1); ,
           return 0.; )

// Torque free head: dpsi/dl(0) = 0
create_bc( BCa_torque_free_head , 0 , D1 ,
           return D1(i_ev,args.j1); ,
           return 0.; )

// Force free tail: d2psi/dl2(1) = 0
create_bc( BCa_force_free_tail , p.N-1 , D2 ,
           return D2(i_ev,args.j1); ,
           return 0.; )

// Torque free tail: dpsi/dl(1) = 0
create_bc( BCa_torque_free_tail , p.N-1 , D1 ,
           return D1(i_ev,args.j1); ,
           return 0.; )

// Pivoting head: d3psi/dl3 = f*dpsi/dl - d2psi/dl2*dpsi/dl
create_bc( BCa_no_head_translation , 0 , D3 ,
           return D3(i_ev,args.j1); ,
           return fil.tension[i_ev]*fil.du[i_ev]; )

// Clamped head: psi(0) = 0
create_bc( BCa_no_head_rotation , 0 , D0 ,
           return D0(i_ev,args.j1); ,
	        return 0; )

// Swimming head: 
create_bc( BCa_drag_force_head , 0 , D3+D2 ,
           return p.head_drag * D3(i_ev,args.j1) - D2(i_ev,args.j1); ,
           return p.head_drag * fil.du[i_ev]*fil.tension[i_ev]; )

create_bc( BCa_drag_torque_head , 0 , D4+D1 ,
           return p.head_rot_drag*D4(i_ev,args.j1) + D1(i_ev,args.j1); ,
           return p.head_rot_drag * (
                    p.gamma * fil.du[i_ev] * ( fil.du[i_ev] * fil.d2u[i_ev] - fil.f[i_ev] )
                    + fil.d2u[i_ev] * fil.tension[i_ev]
                    + (p.gamma+1) * fil.du[i_ev] * D1(fil.tension,i_ev)
                    ); )

// Boundary conditions for the tension (Camalet)
// =============================================

// Force free head: T(0) = 0
create_bc( BCtc_force_free_head , 0 , D0 ,
           return D0(i_ev,args.j1); ,
           return 0; )

// Force free tail: T(0) = 0
create_bc( BCtc_force_free_tail , p.N-1 , D0 ,
           return D0(i_ev,args.j1); ,
           return 0; )

// Pivoting head: dT/dl = f*dpsi/dl - d2psi/dl2*dpsi/dl
create_bc( BCtc_no_head_translation , 0 , D1 ,
           return D1(i_ev,args.j1); ,
	         return fil.f[i_ev] - fil.du[i_ev] * fil.d2u[i_ev]; )

// Swimming head: dT/dl = f*dpsi/dl - d2psi/dl2*dpsi/dl
create_bc( BCtc_drag_force_head , 0 , D0+D1 ,
           return p.gamma * p.head_drag * D1(i_ev,args.j1) - D0(i_ev,args.j1); ,
           return p.gamma * p.head_drag * ( fil.f[i_ev] - fil.du[i_ev] * fil.d2u[i_ev] ); )


#undef create_bc

//------------------------------------------------------------------------

template <class BC_base_type, class tag_BC> class BC_type;

template <class BC_base_type> class BC_type <BC_base_type,tag_free_head>
{
public:
	typedef BCa_force_free_head <BC_base_type> angle1;
	typedef BCa_torque_free_head <BC_base_type> angle2;
	typedef BCtc_force_free_head <BC_base_type> tension_camalet;
};

template <class BC_base_type> class BC_type <BC_base_type,tag_pivoting_head>
{
public:
	typedef BCa_no_head_translation <BC_base_type> angle1;
	typedef BCa_torque_free_head <BC_base_type> angle2;
	typedef BCtc_no_head_translation <BC_base_type> tension_camalet;
};

template <class BC_base_type> class BC_type <BC_base_type,tag_clamped_head>
{
public:
	typedef BCa_no_head_translation <BC_base_type> angle1;
	typedef BCa_no_head_rotation <BC_base_type> angle2;
	typedef BCtc_no_head_translation <BC_base_type> tension_camalet;
};

template <class BC_base_type> class BC_type <BC_base_type,tag_swimming_head>
{
public:
	typedef BCa_drag_force_head <BC_base_type> angle1;
	typedef BCa_drag_torque_head <BC_base_type> angle2;
	typedef BCtc_drag_force_head <BC_base_type> tension_camalet;
};

#endif
