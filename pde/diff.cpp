#include "diff.h"


//namespace Diff
//{

using std::vector;
using std::ostream;
using std::cout;
using std::endl;


void Row::shift_indices(const int & i)
	{
	for (size_t k=0; k<indices.size(); k++) indices[k] += i;
	}


ostream & operator << (std::ostream & os, Row row)
	{
	for (size_t k=0; k<row.indices.size(); k++)
		{
		os << "(" << row.indices[k] << "," << row.values[k] << ")     ";
		}
	return os;
	}


//_______________________________________________________________________________


unsigned int Diff::N; // = E_Diff::N;


Row & Diff::get_row (const unsigned int & i)
	{
	for (size_t j=0; j<boundary_ids.size(); j++)
		{
		if ( boundary_ids[j] == i ) return boundary_rows[j];
		}
	return bulk_row;
	}

const Row & Diff::get_row (const unsigned int & i) const
	{
	for (size_t j=0; j<boundary_ids.size(); j++)
		{
		if ( boundary_ids[j] == i ) return boundary_rows[j];
		}
	return bulk_row;
	}

double Diff::operator () (const unsigned int & i, const unsigned int & j) const
	{
	const Row & row = get_row(i);
	for (size_t k=0; k<row.indices.size(); k++)
		{
		if ( row.indices[k] == (int) j - (int) i ) return row.values[k];
		}
	return 0;
	}

unsigned int Diff::get_nnz (const unsigned int i) const
	{
	return get_row(i).indices.size();
	}

int Diff::get_k_min (const unsigned int i) const
	{
	return get_row(i).indices.back()-i;
	}

int Diff::get_k_max (const unsigned int i) const
	{
	return get_row(i).indices.front()-i;
	}


void Diff::print_dense(std::ostream & os) const
	{
	os << std::setiosflags(std::ios::fixed) << std::setprecision(3);
	for (unsigned int i=0; i<N; i++)
		{
		for (unsigned int j=0; j<N; j++) os << std::setw(10) << (*this)(i,j) << " ";
		os << endl;
		}
	}

ostream & operator << (std::ostream & os, Diff diff)
	{
	diff.print_dense(os);
	return os;
	}


//---------------------------------------------------------------------------------


// Definitions of the elementary differential operators
// (i.e. first, second... derivative)
//namespace E_Diff
//{
const Diff D0 ( { {}, {}, { vector <int> ({0}), vector <double> ({1}) } } );
Diff D1;
Diff D2;
Diff D3;
Diff D4;
void diff_initialize(const Parameters & par)
	{
	Diff::N = par.N;
	const int & N = par.N;
	D1 = { {0,N-1u},
		    { Row( { vector <int> ({ 0, 1, 2}), vector <double> ({-1.5, 2,-0.5}) } ),
		      Row( { vector <int> ({-2,-1, 0}), vector <double> ({ 0.5,-2, 1.5}) } ) },
		    { vector <int> ({-1,1}), vector <double> ({-0.5,0.5}) } };
	D2 = { {0,N-1u},
		    { Row( { vector <int> ({ 0, 1, 2, 3}), vector <double> ({2,-5,4,-1}) } ),
		      Row( { vector <int> ({-3,-2,-1, 0}), vector <double> ({-1,4,-5,2}) } ) },
		    { vector <int> ({-1,0,1}), vector <double> ({1,-2,1}) } };
	D3 = { {0,1,N-2u,N-1u},
		    { Row( { vector <int> ({ 0, 1, 2, 3, 4}), vector <double> ({-2.5,9,-12,7,-1.5}) } ),
		      Row( { vector <int> ({-1, 0, 1, 2, 3}), vector <double> ({-1.5,5,-6,3,-0.5}) } ),
		      Row( { vector <int> ({-3,-2,-1, 0, 1}), vector <double> ({0.5,-3,6,-5,1.5}) } ),
		      Row( { vector <int> ({-4,-3,-2,-1, 0}), vector <double> ({1.5,-7,12,-9,2.5}) } ) },
		    { vector <int> ({-2,-1,1,2}), vector <double> ({-0.5,1,-1,0.5}) } };
	D4 = { {0,1,N-2u,N-1u},
		    { Row( { vector <int> ({ 0, 1, 2, 3, 4, 5}), vector <double> ({3,-14,26,-24,11,-2}) } ),
		      Row( { vector <int> ({-1, 0, 1, 2, 3, 4}), vector <double> ({2,-9,16,-14,6,-1})   } ),
		      Row( { vector <int> ({-4,-3,-2,-1, 0, 1}), vector <double> ({-1,6,-14,16,-9,2})   } ),
		      Row( { vector <int> ({-5,-4,-3,-2,-1, 0}), vector <double> ({-2,11,-24,26,-14,3}) } ) },
		    { vector <int> ({-2,-1,0,1,2}), vector <double> ({1,-4,6,-4,1}) } };
	D1 /= par.dl;
	D2 /= par.dl2;
	D3 /= par.dl3;
	D4 /= par.dl4;
	}
//}


//------------------------------------------------------------------------------

// Additive inverse
Row operator - (Row row)
	{
	for (size_t j=0; j<row.values.size(); j++) row.values[j] = -row.values[j];
	return row;
	}

Diff operator - (Diff diff)
	{
	for (size_t i=0; i<diff.boundary_rows.size(); i++) diff.boundary_rows[i] = -diff.boundary_rows[i];
	diff.bulk_row = -diff.bulk_row;
	return diff;
	}

// Multiplication by a scalar
Row & operator *= (Row & row, const double & f)
	{
	for (size_t j=0; j<row.values.size(); j++) row.values[j] *= f;
	return row;
	}

Diff & operator *= (Diff & diff, const double & f)
	{
	for (size_t i=0; i<diff.boundary_rows.size(); i++) diff.boundary_rows[i] *= f;
	diff.bulk_row *= f;
	return diff;
	}

Diff operator * (Diff diff, const double & f)
	{
	return diff *= f;
	}

Diff operator * (const double & f, Diff diff)
	{
	return diff *= f;
	}

// Division by a scalar
Row & operator /= (Row & row, const double & f)
	{
	for (size_t j=0; j<row.values.size(); j++) row.values[j] /= f;
	return row;
	}

Diff & operator /= (Diff & diff, const double & f)
	{
	for (size_t i=0; i<diff.boundary_rows.size(); i++) diff.boundary_rows[i] /= f;
	diff.bulk_row /= f;
	return diff;
	}

Diff operator / (Diff diff, const double & f)
	{
	return diff /= f;
	}

// Addition and subtraction

template <class Op> Row add  (const Row   & row1 , const Row   & row2 )
	{
	Row row;
	unsigned int k1 = 0, k2 = 0;
	int j1 = row1.indices[0], j2 = row2.indices[0];
	while ( k1 < row1.indices.size() && k2 < row2.indices.size() )
		{
		if ( j1 < j2 )
			{
			row.indices.push_back(j1);
			row.values.push_back( Op() ( row1.values[k1] , 0. ) );
			++k1;
			j1 = row1.indices[k1];
			}
		else if ( j1 > j2 )
			{
			row.indices.push_back(j2);
			row.values.push_back( Op() ( 0. , row2.values[k2] ) );
			++k2;
			j2 = row2.indices[k2];
			}
		else // j1 == j2
			{
			row.indices.push_back(j1);
			row.values.push_back( Op() ( row1.values[k1] , row2.values[k2] ) );
			++k1;
			j1 = row1.indices[k1];
			++k2;
			j2 = row2.indices[k2];
			}
		}
	while ( k1 < row1.indices.size() )
		{
		row.indices.push_back(j1);
		row.values.push_back( Op() ( row1.values[k1] , 0. ) );
		++k1;
		j1 = row1.indices[k1];
		}
	while ( k2 < row2.indices.size() )
		{
		row.indices.push_back(j2);
		row.values.push_back( Op() ( 0. , row2.values[k2] ) );
		++k2;
		j2 = row2.indices[k2];
		}
	return row;
	}


template <class Op> Diff add (const Diff & diff1, const Diff & diff2)
	{
	Diff diff;
	unsigned int k1 = 0, k2 = 0;
	unsigned int i1 = 0, i2 = 0;
	if ( diff1.boundary_ids.size() > 0 ) i1 = diff1.boundary_ids[0];
	if ( diff2.boundary_ids.size() > 0 ) i2 = diff2.boundary_ids[0];
	while ( k1 < diff1.boundary_ids.size() && k2 < diff2.boundary_ids.size() )
		{
		if ( i1 < i2 )
			{
			diff.boundary_ids.push_back(i1);
			diff.boundary_rows.push_back( add <Op> ( diff1.boundary_rows[k1] , diff2.bulk_row ) );
			++k1;
			i1 = diff1.boundary_ids[k1];
			}
		else if ( i1 > i2 )
			{
			diff.boundary_ids.push_back(i2);
			diff.boundary_rows.push_back( add <Op> ( diff2.boundary_rows[k2] , diff1.bulk_row ) );
			++k2;
			i2 = diff2.boundary_ids[k2];
			}
		else // i1 == i2
			{
			diff.boundary_ids.push_back(i1);
			diff.boundary_rows.push_back( add <Op> ( diff1.boundary_rows[k1] , diff2.boundary_rows[k2] ) );
			++k1;
			i1 = diff1.boundary_ids[k1];
			++k2;
			i2 = diff2.boundary_ids[k2];
			}
		}	
	while ( k1 < diff1.boundary_ids.size() )
		{
		diff.boundary_ids.push_back(i1);
		diff.boundary_rows.push_back( add <Op> ( diff1.boundary_rows[k1] , diff2.bulk_row ) );
		++k1;
		i1 = diff1.boundary_ids[k1];
		}
	while ( k2 <diff2.boundary_ids.size() )
		{
		Row r = diff2.boundary_rows[k2];
		diff.boundary_ids.push_back(i2);
		diff.boundary_rows.push_back( add <Op> ( diff2.boundary_rows[k2] , diff1.bulk_row ) );
		++k2;
		i2 = diff2.boundary_ids[k2];
		}
	diff.bulk_row = add <Op> ( diff1.bulk_row , diff2.bulk_row );
	return diff;
	}

Row operator + (const Row & row1, const Row & row2)
	{
	return add < std::plus <double> > (row1,row2);
	}

Row operator - (const Row & row1, const Row & row2)
	{
	return add < std::minus <double> > (row1,row2);
	}

Diff operator + (const Diff & diff1, const Diff & diff2)
	{
	return add < std::plus <double> > (diff1,diff2);
	}

Diff operator - (const Diff & diff1, const Diff & diff2)
	{
	return add < std::minus <double> > (diff1,diff2);
	}


//}

