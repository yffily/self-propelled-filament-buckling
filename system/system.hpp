
// ================
// class Integrator
// ================

template <class TPoint, class tag_BC, class tag_unknown>
Integrator<TPoint,tag_BC,tag_unknown>::Integrator (Parameters & par, Filament_type & filament, const double & t) : N(par.N), 
       pde_position(par), pde_tension(par), 
       solver_position(pde_position,filament,t), solver_tension(pde_tension,filament,t), 
       rhs(N,0.), rhs_last(N,0.), rhs_next(N,0.), 
       u_last(filament.u), d2_max(0),
       full_state({Item_type("u_last",&u_last),Item_type("rhs_last",&rhs_last)})
	{
	}

template <class TPoint, class tag_BC, class tag_unknown>
void Integrator<TPoint,tag_BC,tag_unknown>::initialize_rhs (Filament_type & filament, const double & t)
	{
	filament.get_derivatives(t);
	get_tension(filament,t);
	u_last = filament.u;
	for (size_t i=0; i<N; i++) rhs_last[i] = pde_position.rhs(i,filament,t);
	}

template <class TPoint, class tag_BC, class tag_unknown>
void Integrator<TPoint,tag_BC,tag_unknown>::get_tension (Filament_type & filament, const double & t)
	{
	solver_tension.update_lhs(filament,t);
	solver_tension.solve(filament,filament.tension,t);
	}

template <class TPoint, class tag_BC, class tag_unknown>
void Integrator<TPoint,tag_BC,tag_unknown>::step_forward (Filament_type & filament, const double & t)
	{
	for (size_t i=0; i<N; i++)
		{
		rhs[i] = pde_position.rhs(i,filament,t);
		if ( pde_position.boundary_id(i) >= 0 ) rhs_next[i] = rhs[i];
		else
			{
			rhs_next[i] = 2 * filament.u[i] - 0.5 * u_last[i]
		                   + filament.p->dt * ( 2 * rhs[i] - rhs_last[i] );
			}
		}
	solver_position.update_lhs(filament,t);
	solver_position.solve(rhs_next,filament.u,t);
	filament.move_head();
	d2_max = std::max( d2_max , get_max_displacement2(filament) );
	swap(rhs,rhs_last);
	u_last = filament.u;
	
//	solver_tension.print(cout,filament,t); cout << endl;
//	cout << filament.tension << endl << endl;
//	solver_position.print(cout,filament,t); cout << endl;
//	cout << filament.u << endl << endl;
//	exit(3);
	}

template <class TPoint, class tag_BC, class tag_unknown>
double Integrator<TPoint,tag_BC,tag_unknown>::get_max_displacement2 (Filament_type & filament)
	{
	double d2 = 0.;
	for (size_t i=0; i<N; i++) d2 = std::max( d2 , norm2(filament.u[i]-u_last[i]) );
	return d2;
	}


// ============
// class System
// ============

template <class TPoint, class tag_BC, class tag_unknown>
System<TPoint,tag_BC,tag_unknown>::System (Parameters & par) : step(0), t(0), filament(par), integrator(par,filament,t)
	{
	integrator.initialize_rhs(filament,t);
//	integrator.solver_tension.print(cout,filament,t);
//	integrator.solver_position.print(cout,filament,t);
	}

template <class TPoint, class tag_BC, class tag_unknown>
void System<TPoint,tag_BC,tag_unknown>::evolve ()
	{
	integrator.step_forward(filament,t);
	t += filament.p->dt;
	filament.get_derivatives(t);
	integrator.get_tension(filament,t);
	}

template <class TPoint, class tag_BC, class tag_unknown>
void System<TPoint,tag_BC,tag_unknown>::print_full (ostream & os) const
	{
	os << std::setprecision(16);
	os << "# step" << endl << step << endl << endl;
	os << "# time" << endl << t << endl << endl;
	os << "# head position" << endl << filament.position[0] << endl << endl;
	const typename Filament_type::Full_state_type & filament_state = filament.full_state;
	for (auto i=filament_state.begin(); i!=filament_state.end(); i++)
		{
		os << "# " << i->first << endl;
		for (auto j=i->second->begin(); j!=i->second->end(); j++) os << *j << endl;
		os << endl;
		}
	const typename Integrator_type::Full_state_type & integrator_state = integrator.full_state;
	for (auto i=integrator_state.begin(); i!=integrator_state.end(); i++)
		{
		os << "# " << i->first << endl;
		for (auto j=i->second->begin(); j!=i->second->end(); j++) os << *j << endl;
		os << endl;
		}
	}

template <class TPoint, class tag_BC, class tag_unknown>
void System<TPoint,tag_BC,tag_unknown>::read_full (istream & is)
	{
	std::string line;
	getline(is,line);           // grab the comment line
	is >> step;
	getline(is,line);           // grab the end-of-line character after t
	getline(is,line);           // grab the empty line
	getline(is,line);           // grab the comment line
	is >> t;
	getline(is,line);           // grab the end-of-line character after t
	getline(is,line);           // grab the empty line
	getline(is,line);           // ...
	is >> filament.position[0];
	getline(is,line);
	getline(is,line);
	typename Filament_type::Full_state_type & filament_state = filament.full_state;
	for (auto i=filament_state.begin(); i!=filament_state.end(); i++)
		{
		getline(is,line);
		for (auto j=i->second->begin(); j!=i->second->end(); j++) is >> *j;
		getline(is,line);
		getline(is,line);
		}
	typename Integrator_type::Full_state_type & integrator_state = integrator.full_state;
	for (auto i=integrator_state.begin(); i!=integrator_state.end(); i++)
		{
		getline(is,line);
		for (auto j=i->second->begin(); j!=i->second->end(); j++) is >> *j;
		getline(is,line);
		getline(is,line);
		}
	filament.get_position();
	filament.get_derivatives(t);
	integrator.solver_tension.update_lhs(filament,t);
	integrator.get_tension(filament,t);
	}
