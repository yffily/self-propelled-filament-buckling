#ifndef diff_h
#define diff_h
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <tuple>
#include <cmath>  // contains abs function
#include "main/util.h"
#include "parameters/parameters.h"


//namespace Diff
//{


class Row
{
public:
	std::vector <int> indices;
	std::vector <double> values;
	
	void shift_indices(const int & i);
};

std::ostream & operator << (std::ostream & os, Row row);


class Diff
{
public:
	static unsigned int N;        // size of the grid
	
	std::vector <unsigned int> boundary_ids; 
	std::vector <Row> boundary_rows; 
	Row bulk_row;
	
	Row & get_row (const unsigned int & i);
	const Row & get_row (const unsigned int & i) const;
	double operator () (const unsigned int & i, const unsigned int & j) const;
	template <class TLHS_args> double operator () (const TLHS_args & args) const;
	template <class TPoint> TPoint operator () (const std::vector <TPoint> & u, const unsigned int i) const;
	unsigned int get_nnz (const unsigned int i) const;          // number of non-zero elements
	int get_k_min (const unsigned int i) const;        // smallest value of j-i that returns a non-zero element
	int get_k_max (const unsigned int i) const;        // smallest value of j-i that returns a non-zero element
	void print_dense(std::ostream & os) const;
};

template <class TLHS_args> double Diff::operator () (const TLHS_args & args) const
	{
	return (*this)(args.i1,args.j1);
	};

template <class TPoint> TPoint Diff::operator () (const std::vector <TPoint> & u, const unsigned int i) const
	{
	TPoint v = 0.;
	const Row & row = get_row(i);
	for (unsigned int j=0; j<row.indices.size(); j++) v += row.values[j] * u[i+row.indices[j]];
	return v;
	};

std::ostream & operator << (std::ostream & os, Diff diff);


//------------------------------------------------------------------------------

//namespace E_Diff
//{
extern const Diff D0;
extern Diff D1;
extern Diff D2;
extern Diff D3;
extern Diff D4;
//extern void diff_initialize(unsigned int & N);
void diff_initialize(const Parameters & par);
//};


//------------------------------------------------------------------------------


Row   operator - (Row row);
Row   operator - (Row row1);

Row  & operator *= (Row   & row,  const double & f);
Diff & operator *= (Diff & diff, const double & f);
Diff   operator *  (Diff   diff, const double & f);
Diff   operator *  (const double & f, Diff   diff);

Row  & operator /= (Row   & row,  const double & f);
Diff & operator /= (Diff & diff, const double & f);
Diff   operator /  (Diff   diff, const double & f);

template <class Op> Row add   (const Row   & row1 , const Row   & row2 );
template <class Op> Diff add (const Diff & diff1, const Diff & diff2);
Row     operator +  (const Row   & row1 , const Row   & row2 );
Row     operator -  (const Row   & row1 , const Row   & row2 );
Diff   operator +  (const Diff & diff1, const Diff & diff2);
Diff   operator -  (const Diff & diff1, const Diff & diff2);

//Diff operator /  (const Diff & diff, const double & f);

//Diff test (const Diff & diff1, const Diff & Diff);

//};

#endif
