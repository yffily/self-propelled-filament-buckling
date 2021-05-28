
// ==============
// class Filament
// ==============

template <class TPoint>
Filament<TPoint,tag_angle>::Filament (Parameters & par) : N(par.N), p(&par), 
                      position(N,0.), u(N,0.), 
                      du(N,0.), d2u(N,0.), d3u(N,0.), d4u(N,0.), 
                      tension(N,0.), 
                      f(N,0.), df(N,0.),
                      full_state({Item_type("u",&u)})
	{
	// Initial condition
	// =================
	for (size_t i=0; i<N; i++) u[i] = p->A0*cos(p->k0_*p->dl*i);
	double t = 0;
	get_derivatives(t);
	}


template <class TPoint>
void Filament<TPoint,tag_angle>::get_derivatives (const double & t)
	{
	for (size_t i=0; i<N; i++)
		{
		du[i]  = D1(u,i);
		d2u[i] = D2(u,i);
		d3u[i] = D3(u,i);
		d4u[i] = D4(u,i);
		f[i]   = internal_force(i,t,*p);
//		f[i]   = p->A * sin(p->k*i*p->dl-p->w*t);
		}
	for (size_t i=0; i<N; i++)
		{
		df[i]   = D1(f,i);
		}
	}

template <class TPoint>
void Filament<TPoint,tag_angle>::move_head ()
	{
	TPoint t = unit_vector_from_angle(u[0]);
	TPoint n = t.perp();
	position[0] += p->dt * p->sp_inv * (
	                  ( - d3u[0] + du[0]*tension[0] ) * n
	                 + p->gamma * ( du[0] * d2u[0] + D1(tension,0) - f[0]  ) * t
	               );
	}

template <class TPoint>
void Filament<TPoint,tag_angle>::get_position ()
	{
	for (size_t i=1; i<N; i++)
		{
		position[i] = position[i-1] + p->dl*unit_vector_from_angle(u[i-1]);
		}
	}

template <class TPoint>
double Filament<TPoint,tag_angle>::get_curvature_energy ()
	{
	double e = 0;
	for (size_t i=0; i<N; i++) { e += du[i]*du[i]; }
	return e*p->dl;
	}

//---------------------------------------------------------------------------

template <class TPoint, class tag_unknown>
ostream & operator << (ostream & os, Filament <TPoint,tag_unknown> & filament)
	{
	filament.get_position();
	for (size_t i=0; i<filament.N; i+=filament.p->nl_print)
		{
		os << filament.position[i];
		os << sep << filament.tension[i];
		if ( std::is_same<tag_unknown,tag_angle>::value ) os << sep << filament.u[i];
		os << endl;
		}
	return os;
	}

template <class TPoint, class tag_unknown>
istream & operator >> (istream & is, Filament <TPoint,tag_unknown> & filament)
	{
	unsigned int i = 0;
	std::string line;
	while ( getline(is,line) )
		{
		std::istringstream iss(line);
		iss >> filament.position[i];
		iss >> filament.tension[i];
		if ( std::is_same<tag_unknown,tag_angle>::value ) iss >> filament.u[i];
		i++;
		}
	if ( i != filament.N )
		{
		std::cerr << "[Error] Input file has " << i << " lines; filament expected " << filament.N << "." << std::endl;
		exit(2);
		}
	return is;
	}


