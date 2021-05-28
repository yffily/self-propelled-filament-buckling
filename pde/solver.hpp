
// ============
// class Solver
// ============


template <class TPDE>
Solver<TPDE>::Solver(TPDE & pde, Filament_type & filament, const double & t) : CrsMatrix_policy<TPDE>(), pde(pde), map((int)M,0,Epetra_SerialComm()), X(map,N*D/M), B(map,N*D/M), A(Copy, map, &pde.get_nnz_vec()[0] ), solver_allocated(true), width(8), precision(3)
	{
	this->build_matrix(A,pde,filament,t);
	
	// Set up the linear problem.
	problem.SetOperator(&A);
	problem.SetLHS(&X);
	problem.SetRHS(&B);
	
	// Initialize the solver and perform the step that only depends on the sparsity pattern
	solver = factory.Create("Klu", problem);
	solver->SymbolicFactorization();  // only uses sparsity pattern of A
	};


template <class TPDE>
void Solver<TPDE>::update_lhs(const Filament_type & filament, const double & t)
	{
	update_lhs(filament,t,tag_lhs());
	}

template <class TPDE>
void Solver<TPDE>::update_lhs(const Filament_type & filament, const double & t, tag_non_constant_lhs)
	{
	this->update_matrix(A,pde,filament,t);
	}


// solve the equation using the current state in filament, write the output in u
template <class TPDE>
void Solver<TPDE>::solve(const Filament_type & filament, std::vector <Unknown_type> & u, const double & t) // Solve A.X=B
	{
	for (size_t i=0; i<u.size(); i++) v2m(pde.rhs(i,filament,t),B,i);
	solver->NumericFactorization();
	solver->Solve();
	for (size_t i=0; i<u.size(); i++) m2v(X,u[i],i);
	};


// solve the equation using a pre-computed rhs, write the output in u
template <class TPDE>
void Solver<TPDE>::solve(const std::vector <Unknown_type> & rhs, std::vector <Unknown_type> & u, const double & t) // Solve A.X=B
	{
	for (size_t i=0; i<u.size(); i++) v2m(rhs[i],B,i);
	solver->NumericFactorization();
	solver->Solve();
	for (size_t i=0; i<u.size(); i++) m2v(X,u[i],i);
	};


template <class TPDE>
auto Solver<TPDE>::apply_lhs(const std::vector <Unknown_type> & u) -> std::vector <Unknown_type>
	{
	for (size_t i=0; i<u.size(); i++) v2m(u[i],X,i);
	A.Multiply(false,X,B);
	std::vector <Unknown_type> v(u.size());
	for (size_t i=0; i<u.size(); i++) m2v(B,v[i],i);
	return v;
	}


//______________________________________________________________________________________________

// printing

template <class TPDE>
void Solver<TPDE>::print_matrix(std::ostream & os)
	{
//	os << std::fixed << std::setprecision(precision);
//	os << std::scientific << std::setprecision(precision);
	for (unsigned int i=0; i<M; i++)
		{
		int n;
		int * indices;
		double * values;
		A.ExtractMyRowView(i,n,values,indices);
		for (int j=0; j<indices[0]; j++) os << std::setw(width) << 0. << " ";
		for (int k=0; k<n-1; k++)
			{
			os << std::setw(width) << values[k] << " ";
			for (int j=indices[k]+1; j<indices[k+1]; j++) os << std::setw(width) << 0. << " ";
			}
		os << std::setw(width) << values[n-1] << " ";
		for (unsigned int j=indices[n-1]+1; j<M; j++) os << std::setw(width) << 0. << " ";
		os << endl;
		}
	}


template <class TPDE>
void Solver<TPDE>::print(std::ostream & os, const Filament_type & filament, double t)
	{
//	os << std::fixed << std::setprecision(precision);
//	os << std::scientific << std::setprecision(precision);
	for (unsigned int i=0; i<M; i++)
		{
		int n;
		int * indices;
		double * values;
		A.ExtractMyRowView(i,n,values,indices);
		for (int j=0; j<indices[0]; j++) os << std::setw(width) << 0. << " ";
		for (int k=0; k<n-1; k++)
			{
			os << std::setw(width) << values[k] << " ";
			for (int j=indices[k]+1; j<indices[k+1]; j++) os << std::setw(width) << 0. << " ";
			}
		os << std::setw(width) << values[n-1] << " ";
		for (unsigned int j=indices[n-1]+1; j<M; j++) os << std::setw(width) << 0. << " ";
		os << sep+"| ";
		for (int l=0; l<B.NumVectors(); l++) os << std::setw(width) << B[l][i] << " ";
		os << endl;
		}
	}
