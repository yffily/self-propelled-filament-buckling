#ifndef solver_policies_h
#define solver_policies_h
#include <iostream>
#include <vector>
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "geomvec/geomvec2.h"
#include "geomvec/geomvec3.h"
#include "main/traits.h"
#include "main/tags.h"
#include "pde/pde.h"


template <class TPoint, class tag_component_expansion>
class Multivector_policy
{
public:
	// copy from a vector of points to a multivector (v2m) and vice versa (m2v)
	void v2m (TPoint u, Epetra_MultiVector & X, const size_t & i);
	void m2v (const Epetra_MultiVector & X, TPoint & u, const size_t & i);
	TPoint m2v (const Epetra_MultiVector & X, const size_t & i)
		{
		TPoint u;
		m2v(X,u,i);
		return u;
		};
};

//-----------------------------------------------------------

template <class TPDE>
class CrsMatrix_policy
{
public:
	typedef typename TPDE::Filament_type Filament_type;
	typedef typename TPDE::Unknown_type Unknown_type;
	typedef typename TPDE::tag_component_expansion tag_component_expansion;
//	typedef LHS_arguments <Filament_type,tag_component_expansion> LHS_arg_type;
	typedef typename TPDE::LHS_arg_type LHS_arg_type;
	
	const unsigned int & N;
	const unsigned int & D;
	unsigned int M;
	
	CrsMatrix_policy();
	CrsMatrix_policy(Parameters & p);
	void build_matrix(Epetra_CrsMatrix & A, const TPDE & pde, Filament_type & filament, const double & t);
	void build_matrix(Epetra_CrsMatrix & A, const TPDE & pde, Filament_type & filament, const double & t, tag_dont_expand_components);
	void build_matrix(Epetra_CrsMatrix & A, const TPDE & pde, Filament_type & filament, const double & t, tag_expand_components);
	void update_matrix(Epetra_CrsMatrix & A, const TPDE & pde, const Filament_type & filament, const double & t);
	void update_matrix(Epetra_CrsMatrix & A, const TPDE & pde, const Filament_type & filament, const double & t, tag_dont_expand_components);
	void update_matrix(Epetra_CrsMatrix & A, const TPDE & pde, const Filament_type & filament, const double & t, tag_expand_components);
	
};

//-----------------------------------------------------------



#include "solver_policies.hpp"

#endif
