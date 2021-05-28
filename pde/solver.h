#ifndef solver_h
#define solver_h
#include <iostream>
#include <vector>
#include "Epetra_MultiVector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "geomvec/geomvec2.h"
#include "geomvec/geomvec3.h"
#include "main/traits.h"
#include "main/tags.h"
#include "solver_policies.h"


template <class TPDE> class Solver : public Multivector_policy <typename TPDE::Unknown_type, typename TPDE::tag_component_expansion>, public CrsMatrix_policy <TPDE>
{
public:
	typedef typename TPDE::Filament_type Filament_type;
	typedef typename TPDE::Unknown_type Unknown_type;
	typedef typename TPDE::tag_lhs tag_lhs;
	typedef typename TPDE::tag_component_expansion tag_component_expansion;
	
	using Multivector_policy<typename TPDE::Unknown_type,typename TPDE::tag_component_expansion>::v2m;
	using Multivector_policy<typename TPDE::Unknown_type,typename TPDE::tag_component_expansion>::m2v;
	using CrsMatrix_policy<TPDE>::N;
	using CrsMatrix_policy<TPDE>::D;
	using CrsMatrix_policy<TPDE>::M;
	
	TPDE pde;
	Epetra_Map map;
	Epetra_MultiVector X;
	Epetra_MultiVector B;
	Epetra_CrsMatrix A;
	Epetra_LinearProblem problem;
	Amesos_BaseSolver * solver;
	Amesos factory;
	bool solver_allocated; // keeps track of whether to call delete in the destructor
	int width, precision;  // (printing) field width and number of decimal places
	
	Solver(TPDE & pde, Filament_type & filament, const double & t);
	~Solver() { if (solver_allocated) delete solver; };
	void update_lhs(const Filament_type & filament, const double & t);
	void update_lhs(const Filament_type & filament, const double & t, tag_constant_lhs) {};
	void update_lhs(const Filament_type & filament, const double & t, tag_non_constant_lhs);
	void solve(const Filament_type & filament, std::vector <Unknown_type> & u, const double & t);
	void solve(const std::vector <Unknown_type> & rhs, std::vector <Unknown_type> & u, const double & t);
	std::vector <Unknown_type> apply_lhs(const std::vector <Unknown_type> & u);
	void print_matrix(std::ostream & os);
	void print(std::ostream & os, const Filament_type & filament, double t);
};


#include "solver.hpp"

#endif
