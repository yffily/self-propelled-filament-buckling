

template <class TPDE_base>
PDE<TPDE_base>::PDE(Parameters & par) : TPDE_base(par)
	{
	// Delete the (potential) boundary terms in lhs_diff that would overlap with a boundary condition.
	// Here we assume that both bc_id and lhs_diff.boundary_ids are sorted.
	for ( size_t k=0; k<bc_id.size(); k++ )
		{
		bc[k]->row.shift_indices( bc[k]->i_ev - bc_id[k] );
		for ( size_t l=0; l<lhs_diff.boundary_ids.size(); l++ )
			{
			if ( lhs_diff.boundary_ids[l] == bc_id[k] )
				{
				lhs_diff.boundary_rows[l] = bc[k]->row;
				break;
				}
			else if ( lhs_diff.boundary_ids[l] > bc_id[k] )
				{
				lhs_diff.boundary_rows.insert(lhs_diff.boundary_rows.begin()+l,bc[k]->row);
				break;
				}
			}
		}
	}


template <class TPDE_base>
int PDE<TPDE_base>::boundary_id (const unsigned int & i) const
	{
	for ( size_t k=0; k<bc_id.size(); k++ ) { if ( i == bc_id[k] ) return k; }
	return -1;
	}

template <class TPDE_base>
double PDE<TPDE_base>::lhs (LHS_arg_type args) const
	{
	int k = boundary_id(args.i1);
	if ( k < 0 ) return lhs_bulk(args);
	return bc[k]->lhs(args);
	}

template <class TPDE_base>
auto PDE<TPDE_base>::rhs (const unsigned int & i, const Filament_type & filament, const double & t) const -> Unknown_type
	{
	int k = boundary_id(i);
	if ( k < 0 ) return TPDE_base::rhs(i,filament,t);
	return bc[k]->rhs(filament,t);
	}

template <class TPDE_base>
int PDE<TPDE_base>::get_nnz() const
	{
	size_t n(0);
	for (size_t i=0; i<Diff::N; i++) n += lhs_diff.get_row(i).indices.size();
	return n;
	}

template <class TPDE_base>
int PDE<TPDE_base>::get_nnz(const unsigned int i) const
	{
	return lhs_diff.get_row(i).indices.size();
	}

template <class TPDE_base>
std::vector <int> PDE<TPDE_base>::get_nnz_vec(tag_expand_components) const
	{
	const unsigned int & N = Diff::N;
	const unsigned int & D = Traits<Unknown_type>::dim;
	std::vector <int> nnz(N*D);
	for (size_t i=0; i<N; i++)
		{
		nnz[2*i] = get_nnz(i) * D; 
		nnz[2*i+1] = get_nnz(i) * D; 
		}
	return nnz;
	}

template <class TPDE_base>
std::vector <int> PDE<TPDE_base>::get_nnz_vec(tag_dont_expand_components) const
	{
	std::vector <int> nnz(Diff::N);
	for (unsigned int i=0; i<Diff::N; i++) nnz[i] = get_nnz(i); 
	return nnz;
	}

//------------------------------------------------------------

// get a row of the matrix -- used to build a CrsMatrix

template <class TPDE_base>
void PDE<TPDE_base>::get_row (const unsigned int & i, const Filament_type & filament, const double & t, unsigned int & n, std::vector <int> & indices, std::vector <double> & values) const
	{
	get_row(i,filament,t,n,indices,values,tag_component_expansion());
	}

template <class TPDE_base>
void PDE<TPDE_base>::get_row (const unsigned int & i, const Filament_type & filament, const double & t, unsigned int & n, std::vector <int> & indices, std::vector <double> & values, tag_dont_expand_components) const
	{
	const Row & row = lhs_diff.get_row(i);
	n = row.indices.size();
	indices.resize(n);
	values.resize(n);
	for (unsigned int k=0; k<n; k++)
		{
		indices[k] = i + row.indices[k];
		values[k]  = lhs(LHS_arg_type({i,indices[k],filament,t}));
		}
	}

template <class TPDE_base>
void PDE<TPDE_base>::get_row (const unsigned int & i, const Filament_type & filament, const double & t, unsigned int & n, std::vector <int> & indices, std::vector <double> & values, tag_expand_components) const
	{
	const unsigned int & D = Traits<Unknown_type>::dim;
	int i1 = i / D;
	int i2 = i % D;
	const Row & row = lhs_diff.get_row(i1);
	n = D * row.indices.size();
	indices.resize(n);
	values.resize(n);
	for (unsigned int k=0; k<n; k++)
		{
		int j1 = i1 + row.indices[ k / D ];
		int j2 = k % D;
		indices[k] = j1 * D + j2;
		values[k]  = lhs(LHS_arg_type({i1,j1,i2,j2,filament,t}));
		}
	}


//------------------------------------------------------------

template <class TPDE_base>
template <class TFilament>
void PDE<TPDE_base>::print_dense(std::ostream & os, const TFilament & filament, double t)
	{
	os << std::setiosflags(std::ios::fixed) << std::setprecision(3);
	for (size_t i=0; i<filament.N; i++)
		{
		for (size_t j=0; j<filament.N; j++) os << std::setw(10) << (*this)(i,j,filament,t) << " ";
		os << sep+sep << rhs(i,filament,t) << std::endl;
		}
	}

