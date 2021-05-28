#ifndef geomvec3_h
#define geomvec3_h
#include <iostream>
#include <fstream>
#include <cmath>
#include "main/util.h"
#include "RNG/RNG_taus.h"


class Geomvec3 {

public:
	static const unsigned int dim = 3;
	double x,y,z;
	
	Geomvec3 ();
	Geomvec3 (const double & x_, const double & y_, const double & z_);
	Geomvec3 (const Geomvec3 & v);
	Geomvec3 (const double & d);
	
	double X(int i);
	void set_X(int i, const double & d);
	double & operator [] (int i);
	const double & operator [] (int i) const;
	
	Geomvec3 & operator = (const Geomvec3 & v);
	
	Geomvec3 & operator += (const Geomvec3 & v);
	Geomvec3 & operator -= (const Geomvec3 & v);
	Geomvec3 & operator *= (const double & d);
	Geomvec3 & operator /= (const double & d);
	
	Geomvec3 operator + (const Geomvec3 & v) const;
	Geomvec3 operator - (const Geomvec3 & v) const;
	
	Geomvec3 operator * (const double & d) const;
	Geomvec3 operator / (const double & d) const;
	
	double operator * (const Geomvec3 & v) const;      // dot product
	Geomvec3 operator ^ (const Geomvec3 & v) const;  // cross product
	
	Geomvec3 & invert_ew (void);   // element wise inversion
	
	bool operator == (const Geomvec3 & v) const;
	bool operator != (const Geomvec3 & v) const;
	
	double norm (void) const;
	double norm2 (void) const;
	Geomvec3 & normalize (void);
	Geomvec3 unit (void) const;
	Geomvec3 unit_or_zero (void) const;

	std::ostream & write(std::ostream & stm, const std::string & left, const std::string & inner, const std::string & right);
};

// same as v.norm() and v.norm2()
double norm (const Geomvec3 & v);
double norm2 (const Geomvec3 & v);

Geomvec3 operator - (const Geomvec3 & v);
Geomvec3 operator * (const double & f, const Geomvec3 & v);

Geomvec3 multiply_ew (const Geomvec3 & v, const Geomvec3 & w);
Geomvec3 divide_ew (const Geomvec3 & v, const Geomvec3 & w);
Geomvec3 sqrt_ew (const Geomvec3 & v);

std::ostream & operator << (std::ostream & stm, const Geomvec3 & v);
std::istream & operator >> (std::istream & stm, Geomvec3 & v);

//Geomvec3 random_vec(RNG_taus & rng, double amplitude[]);
Geomvec3 random_vec(RNG_taus & rng, Geomvec3 amplitude=Geomvec3(1.,1.,1.));
Geomvec3 random_normal_vec(RNG_taus & rng, Geomvec3 amplitude=Geomvec3(1.,1.,1.));


#endif


