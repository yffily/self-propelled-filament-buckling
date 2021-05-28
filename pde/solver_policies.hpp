
// ========================
// class Multivector_policy
// ========================

template <>
void Multivector_policy<double,tag_dont_expand_components>::v2m (double u, Epetra_MultiVector & B, const size_t & i)
	{ B[0][i] = u; };

template <>
void Multivector_policy<Geomvec2,tag_dont_expand_components>::v2m (Geomvec2 u, Epetra_MultiVector & B, const size_t & i)
	{ B[0][i] = u.x; B[1][i] = u.y; };

template <>
void Multivector_policy<Geomvec3,tag_dont_expand_components>::v2m (Geomvec3 u, Epetra_MultiVector & B, const size_t & i)
	{ B[0][i] = u.x; B[1][i] = u.y; B[2][i] = u.z; };


template <>
void Multivector_policy<double,tag_dont_expand_components>::m2v (const Epetra_MultiVector & X, double & u, const size_t & i)
	{ u = X[0][i]; };

template <>
void Multivector_policy<Geomvec2,tag_dont_expand_components>::m2v (const Epetra_MultiVector & X, Geomvec2 & u, const size_t & i)
	{ u.x = X[0][i]; u.y = X[1][i]; };

template <>
void Multivector_policy<Geomvec3,tag_dont_expand_components>::m2v (const Epetra_MultiVector & X, Geomvec3 & u, const size_t & i)
	{ u.x = X[0][i]; u.y = X[1][i]; u.z = X[2][i]; };

//---------------------------------------------------

template <>
void Multivector_policy<double,tag_expand_components>::v2m (double u, Epetra_MultiVector & B, const size_t & i)
	{ B[0][i] = u; }

template <>
void Multivector_policy<Geomvec2,tag_expand_components>::v2m (Geomvec2 u, Epetra_MultiVector & B, const size_t & i)
	{ B[0][2*i] = u.x; B[0][2*i+1] = u.y; }

template <>
void Multivector_policy<Geomvec3,tag_expand_components>::v2m (Geomvec3 u, Epetra_MultiVector & B, const size_t & i)
	{ B[0][3*i] = u.x; B[0][3*i+1] = u.y; B[0][3*i+2] = u.z; }


template <>
void Multivector_policy<double,tag_expand_components>::m2v (const Epetra_MultiVector & X, double & u, const size_t & i)
	{ u = X[0][i]; };

template <>
void Multivector_policy<Geomvec2,tag_expand_components>::m2v (const Epetra_MultiVector & X, Geomvec2 & u, const size_t & i)
	{ u.x = X[0][2*i]; u.y = X[0][2*i+1]; }

template <>
void Multivector_policy<Geomvec3,tag_expand_components>::m2v (const Epetra_MultiVector & X, Geomvec3 & u, const size_t & i)
	{ u.x = X[0][3*i]; u.y = X[0][3*i+1]; u.z = X[0][3*i+2]; }



// ======================
// class CrsMatrix_policy
// ======================

template <class TPDE>
CrsMatrix_policy<TPDE>::CrsMatrix_policy() : N(Diff::N), D(Traits<Unknown_type>::dim)
	{
	if ( std::is_same<tag_component_expansion,tag_expand_components>::value ) M = N * D;
	else M = N;
	}

template <class TPDE>
void CrsMatrix_policy<TPDE>::build_matrix(Epetra_CrsMatrix & A, const TPDE & pde, Filament_type & filament, const double & t)
	{
	for (unsigned int i=0; i<M; i++)
		{
		unsigned int n;
		std::vector <int> indices;
		std::vector <double> values;
		pde.get_row(i,filament,t,n,indices,values);
		A.InsertMyValues((int) i, (int) n, &values[0], &indices[0]);
		}
	A.FillComplete();
	}
	
template <class TPDE>
void CrsMatrix_policy<TPDE>::update_matrix(Epetra_CrsMatrix & A, const TPDE & pde, const Filament_type & filament, const double & t)
	{
	for (unsigned int i=0; i<M; i++)
		{
		unsigned int n;
		std::vector <int> indices;
		std::vector <double> values;
		pde.get_row(i,filament,t,n,indices,values);
		A.ReplaceMyValues(i, n, &values[0], &indices[0]);
		}
	A.FillComplete();
	}
	
	
