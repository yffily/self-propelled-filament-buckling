#ifndef pde_h
#define pde_h
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include "parameters/parameters.h"
#include "main/traits.h"
#include "main/tags.h"
#include "pde/boundary_condition.h"
#include "pde/pde_angle.h"
#include "pde/pde_tension_camalet.h"
#include "pde/pde_angle.h"


template <class TPDE_base>
class PDE : public TPDE_base
{
public:
	typedef typename TPDE_base::Unknown_type Unknown_type;
	typedef typename TPDE_base::Filament_type Filament_type;
	typedef typename TPDE_base::tag_component_expansion tag_component_expansion;
	typedef LHS_arguments<Filament_type,tag_component_expansion> LHS_arg_type;
	
	using TPDE_base::lhs_diff;
	using TPDE_base::lhs_bulk;
	using TPDE_base::bc;
	using TPDE_base::bc_id;
	
	PDE(Parameters & par);
	int boundary_id (const unsigned int & i) const;
	double lhs (LHS_arg_type args) const;
	Unknown_type rhs (const unsigned int & i, const Filament_type & filament, const double & t) const;
	int get_nnz () const;                              // total number of non-zero elements
	int get_nnz (const unsigned int i) const;          // number of non-zero elements
	std::vector <int> get_nnz_vec () const { return get_nnz_vec(tag_component_expansion()); };
	std::vector <int> get_nnz_vec (tag_expand_components) const;
	std::vector <int> get_nnz_vec (tag_dont_expand_components) const;
	
	void get_row (const unsigned int & i, const Filament_type & filament, const double & t, unsigned int & n, std::vector <int> & indices, std::vector <double> & values) const;
	void get_row (const unsigned int & i, const Filament_type & filament, const double & t, unsigned int & n, std::vector <int> & indices, std::vector <double> & values, tag_dont_expand_components) const;
	void get_row (const unsigned int & i, const Filament_type & filament, const double & t, unsigned int & n, std::vector <int> & indices, std::vector <double> & values, tag_expand_components) const;
	
	template <class TFilament>
	void print_dense(std::ostream & os, const TFilament & filament, double t);
};


#include "pde.hpp"

#endif
