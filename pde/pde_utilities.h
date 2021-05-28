#ifndef pde_utilities_h
#define pde_utilities_h
#include <array>
#include "main/traits.h"
#include "main/tags.h"
#include "pde/boundary_condition.h"


// ===================
// class LHS_arguments
// ===================
// Defines the arguments needed to evaluate the lhs of a PDE/ODE.
// (depends on dimensionality and whether components are expanded)

template <class TFilament, class tag_component_expansion> class LHS_arguments;

template <class TFilament> class LHS_arguments <TFilament,tag_dont_expand_components>
{
public:
	const unsigned int & i1;
	const unsigned int & j1;
	const TFilament & filament;
	const double & t;
};

template <class TFilament> class LHS_arguments <TFilament,tag_expand_components>
{
public:
	const unsigned int & i1;
	const unsigned int & j1;
	const unsigned int & i2;
	const unsigned int & j2;
	const TFilament & filament;
	const double & t;
};

//___________________________________________________________________________

// ===========
// class Block
// ===========
// Defines functions for simple block matrices (e.g. identity)
// (returns 1 if components are not expanded)

template <class TFilament, class tag_component_expansion> class Block;

template <class TFilament> class Block <TFilament,tag_dont_expand_components>
{
public:
	typedef LHS_arguments <TFilament,tag_dont_expand_components> LHS_arg_type;
	double identity(LHS_arg_type args) const { return 1; }
};

template <class TFilament> class Block <TFilament,tag_expand_components>
{
public:
	typedef LHS_arguments <TFilament,tag_expand_components> LHS_arg_type;
	double identity(LHS_arg_type args) const { return int(args.i2==args.j2); }
};


#endif
