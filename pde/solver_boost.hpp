
// ============
// class Solver
// ============

template <class TPDE>
Solver<TPDE>::Solver(TPDE & pde, Filament_type & filament, const double & t) : 
                    N(Diff::N), D(Traits<Unknown_type>::dim), M(get_M(pde)),
                    pde(pde), A(M,M,pde.get_nnz()), P(M), X(M), B(M),
                    width(8), precision(3)
	{
	// build the matrix and perform LU factorization
	for (unsigned int i=0; i<M; i++)
		{
		unsigned int n;
		std::vector <int> indices;
		std::vector <double> values;
		pde.get_row(i,filament,t,n,indices,values);
		for (size_t j=0; j<indices.size(); j++) A.push_back(i,indices[j],values[j]);
		}
	
	// make a copy of A for printing before lu_factorize changes it
//	if (std::is_same<typename TPDE::tag_lhs,tag_constant_lhs>::value) 
	A_ = Matrix(A);
	
////	pde.print_dense(clog,filament,t);
//	this->print_matrix(clog);
//	clog << "-------------------------------------------------------------\n";
	
	lu_factorize(A,P);
	
//	this->print_matrix(clog,A_);
//	clog << "-------------------------------------------------------------\n";
//	this->print_matrix(clog,A);
//	clog << "-------------------------------------------------------------\n";
	
//	Matrix B = prod(A,P);
//	auto B = prod(A,P);
//	clog << demangled_type(B) << endl;
////	clog << B.size1() << endl;
//	clog << B(0,0) << endl;
//	clog << B(1,1) << endl;
	
//	Permutation IP = inverse_permutation(P);
//	
//	swap_rows(P,A_);
//	
//	this->print_matrix(clog,A_);
//	clog << "-------------------------------------------------------------\n";
//	
//	
//	exit(3);
	}

template <class TPDE>
unsigned int Solver<TPDE>::get_M(TPDE & pde)
	{
	if ( std::is_same<tag_component_expansion,tag_expand_components>::value ) return Diff::N * Traits<Unknown_type>::dim;
	else return Diff::N;
	}

template <class TPDE>
void Solver<TPDE>::update_lhs(const Filament_type & filament, const double & t)
	{
	update_lhs(filament,t,tag_lhs());
	}

template <class TPDE>
void Solver<TPDE>::update_lhs(const Filament_type & filament, const double & t, tag_non_constant_lhs)
	{
	for (unsigned int i=0; i<M; i++)
		{
		unsigned int n;
		std::vector <int> indices;
		std::vector <double> values;
		pde.get_row(i,filament,t,n,indices,values);
		for (size_t j=0; j<indices.size(); j++) A_(i,indices[j]) = values[j];
		}
	A = A_;
	lu_factorize(A,P);
//	cout << demangled_type(*this) << endl;
//	this->print_matrix(cout,A_); //cout << "\n";
//	cout << "-------------------------------------------------------------\n";
	}


// solve the equation using the current state in filament, write the output in u
template <class TPDE>
void Solver<TPDE>::solve(const Filament_type & filament, std::vector <Unknown_type> & u, const double & t) // Solve A.X=B
	{
	// only for Unkown_type=double for now !!
	for (size_t i=0; i<u.size(); i++) B(i) = pde.rhs(i,filament,t);
	X.assign(B);
	lu_substitute(A,P,X);
	for (size_t i=0; i<u.size(); i++) u[i] = X(i);
	};

// Solve the equation using the current state in filament, write the output in u and v.
// Same as above but with two vectors passed for composite pde.
template <class TPDE>
void Solver<TPDE>::solve(const Filament_type & filament, std::vector <Unknown_type> & u, std::vector <Unknown_type> & v, const double & t) // Solve A.X=B
	{
//	for (size_t i=0; i<u.size(); i++) v2m(pde.rhs(i,filament,t),B,i);
//	for (size_t i=u.size(); i<u.size()+v.size(); i++) v2m(pde.rhs(i,filament,t),B,i);
//	solver->NumericFactorization();
//	solver->Solve();
//	for (size_t i=0; i<u.size(); i++) m2v(X,u[i],i);
//	for (size_t i=0; i<v.size(); i++) m2v(X,v[i],i+u.size());
	};

// solve the equation using a pre-computed rhs, write the output in u
template <class TPDE>
void Solver<TPDE>::solve(const std::vector <Unknown_type> & rhs, std::vector <Unknown_type> & u, const double & t) // Solve A.X=B
	{
	// If the lhs doesn't change, we need to make a fresh copy of A every time
	// otherwise A gets changed by lu_factorize.
	for (size_t i=0; i<u.size(); i++) B(i) = rhs[i];
	X.assign(B);
	lu_substitute(A,P,X);
	for (size_t i=0; i<u.size(); i++) u[i] = X(i);
//	clog << rhs << endl;
//	clog << u << endl << endl;
	};


template <class TPDE>
auto Solver<TPDE>::apply_lhs(const std::vector <Unknown_type> & u) -> std::vector <Unknown_type>
	{
	for (size_t i=0; i<u.size(); i++) X(i) = u[i];
	B = prod(A,X);
	std::vector <Unknown_type> v(u.size());
	for (size_t i=0; i<u.size(); i++) v[i] = B(i);
	return v;
	}


//______________________________________________________________________________________________

// printing

template <class TPDE>
template <class E>
void Solver<TPDE>::print_matrix(std::ostream & os, const ublas::matrix_expression <E> & mat)
//template <class Mat>
//void Solver<TPDE>::print_matrix(std::ostream & os, const Mat & mat)
	{
	os << std::fixed << std::setprecision(precision);
//	os << std::scientific << std::setprecision(precision);
	for (size_t i=0; i<M; i++)
		{
		for (size_t j=0; j<M; j++) os << std::setw(width) << mat()(i,j) << " ";
		os << endl;
		}
	}

template <class TPDE>
template <class E>
void Solver<TPDE>::print_mathematica_matrix(std::ostream & os, const ublas::matrix_expression <E> & mat)
	{
//	os.unsetf(ios_base::floatfield);
	os << std::fixed;
	os << std::setprecision(9);
	os << "{ ";
	for (size_t i=0; i<M; i++)
		{
		os << "{";
		for (size_t j=0; j<M; j++) { os << mat()(i,j); if (j<M-1) os << ","; }
		os << "}";
		if (i<M-1) os << ", ";
		}
	os << " }";
	}


template <class TPDE>
void Solver<TPDE>::print(std::ostream & os, const Filament_type & filament, double t)
	{
//	os << std::fixed << std::setprecision(precision);
//	os << std::scientific << std::setprecision(precision);
	for (size_t i=0; i<M; i++)
		{
		for (size_t j=0; j<M; j++) os << std::setw(width) << A_(i,j) << " ";
		os << sep+"| " << std::setw(width) << B(i);
		os << endl;
		}
	}
