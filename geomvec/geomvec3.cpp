#include "geomvec3.h"
using std::ostream;
using std::istream;
using std::endl;
using std::string;
using std::cerr;

const unsigned int Geomvec3::dim;

Geomvec3::Geomvec3() 
	{}

Geomvec3::Geomvec3(const double & x_, const double & y_, const double & z_) : x(x_), y(y_), z(z_) 
	{}

Geomvec3::Geomvec3(const Geomvec3 & v) : x(v.x), y(v.y), z(v.z) 
	{}

Geomvec3::Geomvec3(const double & d) : x(d), y(d), z(d) 
	{}


double Geomvec3::X(int i)
	{
	if (i==0)      return x;
	else if (i==1) return y;
	else if (i==2) return z;
	else { cerr << "vector has no component #" << i << endl; exit(1); }
	}

void Geomvec3::set_X(int i, const double & d)
	{
	if (i==0)      x=d;
	else if (i==1) y=d;
	else if (i==2) z=d;
	else { cerr << "vector has no component #" << i << endl; exit(1); }
	}

double & Geomvec3::operator [] (int i)
	{
	if (i==0)      return x;
	else if (i==1) return y;
	else if (i==2) return z;
	else { cerr << "vector has no component #" << i << endl; exit(1); }
	}

const double & Geomvec3::operator [] (int i) const
	{
	if (i==0)      return x;
	else if (i==1) return y;
	else if (i==2) return z;
	else { cerr << "vector has no component #" << i << endl; exit(1); }
	}


Geomvec3 & Geomvec3::operator = (const Geomvec3 & v)
	{
	x = v.x; y = v.y; z = v.z; 
	return *this; 
	}

Geomvec3 & Geomvec3::operator += (const Geomvec3 & v)
	{
	x += v.x; y += v.y; z += v.z; 
	return *this; 
	}

Geomvec3 & Geomvec3::operator -= (const Geomvec3 & v)
	{
	x -= v.x; y -= v.y; z -= v.z; 
	return *this; 
	}

Geomvec3 & Geomvec3::operator *= (const double & d)
	{ 
	x *= d; y *= d; z *= d; 
	return *this; 
	}

Geomvec3 & Geomvec3::operator /= (const double & d)
	{
	double a = 1./d;
	x *= a; y *= a; z *= a; 
	return *this; 
	}


Geomvec3 Geomvec3::operator + (const Geomvec3 & v) const
		{ 
		return Geomvec3(x+v.x,y+v.y,z+v.z);
		};

Geomvec3 Geomvec3::operator - (const Geomvec3 & v) const
	{ 
	return Geomvec3(x-v.x,y-v.y,z-v.z);
	};


Geomvec3 Geomvec3::operator * (const double & d) const
	{ 
	return Geomvec3(x*d,y*d,z*d); 
	}

Geomvec3 Geomvec3::operator / (const double & d) const
	{ 
	double a = 1./d;
	return Geomvec3(x*a,y*a,z*a); 
	}


double Geomvec3::operator * (const Geomvec3 & v) const
	{ 
	return x*v.x+y*v.y+z*v.z; 
	}

Geomvec3 Geomvec3::operator ^ (const Geomvec3 & v) const
	{
	return Geomvec3( y*v.z-z*v.y , z*v.x-x*v.z , x*v.y-y*v.x );
	}

Geomvec3 & Geomvec3::invert_ew (void)   // element inversion
	{
	x=1./x; y=1./y; z=1./z;
	return *this;
	}


bool Geomvec3::operator == (const Geomvec3 & v) const
	{
	return ( x==v.x && y==v.y && z==v.z );
	}

bool Geomvec3::operator != (const Geomvec3 & v) const
	{
	return ( x!=v.x || y!=v.y || z!=v.z );
	}


double Geomvec3::norm(void) const
	{ 
	return sqrt(x*x+y*y+z*z);
	};

double Geomvec3::norm2(void) const
	{ 
	return x*x+y*y+z*z;
	};

Geomvec3 & Geomvec3::normalize(void)
	{ 
	double n = 1./sqrt(x*x+y*y+z*z);
	x*=n; y*=n; z*=n;
	return *this;
	};

Geomvec3 Geomvec3::unit(void) const
	{ 
	double n = 1./sqrt(x*x+y*y+z*z);
	return Geomvec3(x*n,y*n,z*n);
	};

Geomvec3 Geomvec3::unit_or_zero(void) const
	{ 
	if ( x==0 && y==0 && z==0 ) return Geomvec3(0.);
	double n = 1./sqrt(x*x+y*y+z*z);
	return Geomvec3(x*n,y*n,z*n);
	};

ostream & Geomvec3::write(ostream& stm, const string & left, const string & inner, const string & right)		// stm << v
	{ 
	stm << left << x << inner << y << inner << z << right;
	return stm; 
	};
	

//-----------------------------------------------------------------

double norm (const Geomvec3 & v)
	{ 
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	}

double norm2 (const Geomvec3 & v)
	{ 
	return v.x*v.x + v.y*v.y + v.z*v.z;
	}


Geomvec3 operator - (const Geomvec3 & v)
	{
	return Geomvec3(-v.x,-v.y,-v.z);
	}

Geomvec3 operator * (const double & f, const Geomvec3 & v)
	{
	return Geomvec3(f*v.x,f*v.y,f*v.z);
	}


Geomvec3 multiply_ew (const Geomvec3 & v, const Geomvec3 & w) // element wise multiplication
	{
	return Geomvec3(v.x*w.x,v.y*w.y,v.z*w.z);
	}

Geomvec3 divide_ew (const Geomvec3 & v, const Geomvec3 & w)   // element wise division
	{
	return Geomvec3(v.x/w.x,v.y/w.y,v.z/w.z);
	}

Geomvec3 sqrt_ew (const Geomvec3 & v)
	{
	return Geomvec3(sqrt(v.x),sqrt(v.y),sqrt(v.z));
	}


ostream & operator << (ostream& stm, const Geomvec3 & v)		// stm << v
	{ 
	stm << v.x << sep << v.y << sep << v.z << sep;
	return stm; 
	};
	
// reading: first character and delimiters have to be exactly one (non number or dot) character.
istream & operator >> (istream & stm, Geomvec3 & v)
	{ 
	read_until_number(stm);
	stm >> v.x;
	read_until_number(stm);
	stm >> v.y;
	read_until_number(stm);
	stm >> v.z;
	return stm; 
	}


Geomvec3 random_vec(RNG_taus & rng, Geomvec3 amplitude) 
	{
	amplitude.x*=rng.get_double();
	amplitude.y*=rng.get_double();
	amplitude.z*=rng.get_double();
	return amplitude;
	}

Geomvec3 random_normal_vec(RNG_taus & rng, Geomvec3 amplitude) 
	{
	amplitude.x*=random_normal(rng);
	amplitude.y*=random_normal(rng);
	amplitude.z*=random_normal(rng);
	return amplitude;
	}



