#ifndef solver_boost_h
#define solver_boost_h
#define BOOST_UBLAS_NDEBUG
#include <iostream>
#include <vector>
#include <iomanip>
#include <initializer_list>
#include "geomvec/geomvec2.h"
#include "geomvec/geomvec3.h"
#include "main/traits.h"
#include "main/tags.h"
#include "pde/pde.h"
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 
#include <boost/numeric/ublas/matrix_sparse.hpp> 
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 

namespace ublas = boost::numeric::ublas; 
using std::cout;
using std::clog;
using std::endl;

template <class TPDE> class Solver
{
public:
	typedef typename TPDE::Filament_type Filament_type;
	typedef typename TPDE::Unknown_type Unknown_type;
	typedef typename TPDE::tag_lhs tag_lhs;
	typedef typename TPDE::tag_component_expansion tag_component_expansion;
	typedef ublas::compressed_matrix <double> Matrix;
	typedef ublas::vector <double> Vector;
	typedef ublas::permutation_matrix <size_t> Permutation; 
	
	unsigned int N, D, M;
	TPDE pde;
	Matrix A, A_;
	Permutation P;
	Vector X, B;
	int width, precision;  // (printing) field width and number of decimal places
	
	Solver(TPDE & pde, Filament_type & filament, const double & t);
	unsigned int get_M(TPDE & pde);
	void update_lhs(const Filament_type & filament, const double & t);
	void update_lhs(const Filament_type & filament, const double & t, tag_constant_lhs) {};
	void update_lhs(const Filament_type & filament, const double & t, tag_non_constant_lhs);
	void solve(const Filament_type & filament, std::vector <Unknown_type> & u, const double & t);
	void solve(const Filament_type & filament, std::vector <Unknown_type> & u, std::vector <Unknown_type> & v, const double & t);
	void solve(const std::vector <Unknown_type> & rhs, std::vector <Unknown_type> & u, const double & t);
	std::vector <Unknown_type> apply_lhs(const std::vector <Unknown_type> & u);
	template <class E> void print_matrix(std::ostream & os, const ublas::matrix_expression <E> & mat);
	void print_matrix(std::ostream & os) { print_matrix(os,A_); };
	template <class E> void print_mathematica_matrix(std::ostream & os, const ublas::matrix_expression <E> & mat);
	void print_mathematica_matrix(std::ostream & os) { print_mathematica_matrix(os,A_); };
	void print(std::ostream & os, const Filament_type & filament, double t);
};


#include "solver_boost.hpp"

#endif
