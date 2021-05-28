#ifndef geomvec2_h
#define geomvec2_h
#include <iostream>
#include <fstream>
#include <cmath>
#include "main/util.h"
#include "RNG/RNG_taus.h"

class Geomvec2 {

public:
	static const unsigned int dim = 2;
	double x,y;
	
	Geomvec2 ();
	Geomvec2 (const double & x_, const double & y_);
	Geomvec2 (const Geomvec2 & v);
	Geomvec2 (const double & d);
	
	double X(int i);
	void set_X(int i, const double & d);
	double & operator [] (int i);
	const double & operator [] (int i) const;
	
	Geomvec2 & operator = (const Geomvec2 & v);
	Geomvec2 & operator = (const double & d);
	
	Geomvec2 & operator += (const Geomvec2 & v);
	Geomvec2 & operator -= (const Geomvec2 & v);
	Geomvec2 & operator *= (const double & d);
	Geomvec2 & operator /= (const double & d);
	
	Geomvec2 operator + (const Geomvec2 & v) const;
	Geomvec2 operator - (const Geomvec2 & v) const;
	
	Geomvec2 operator * (const double & d) const;
	Geomvec2 operator / (const double & d) const;
	
	double operator * (const Geomvec2 & v) const;   // dot product
	double operator ^ (const Geomvec2 & v) const;   // 2D "cross product"
	
	Geomvec2 & invert_ew (void);                    // element wise inversion
	
	bool operator == (const Geomvec2 & v) const;
	bool operator != (const Geomvec2 & v) const;
	
	double norm (void) const;
	double norm2 (void) const;
	Geomvec2 & normalize (void);
	Geomvec2 unit (void) const;
	Geomvec2 unit_or_zero (void) const;
	
	Geomvec2 perp(void) const;              // orthogonal vector
	double angle(void) const;             // angle with x axis

	std::ostream & write(std::ostream & stm, const std::string & left, const std::string & inner, const std::string & right);
};

// same as v.norm() and v.norm2()
double norm (const Geomvec2 & v);
double norm2 (const Geomvec2 & v);

Geomvec2 operator - (const Geomvec2 & v);

Geomvec2 operator * (const double & f, const Geomvec2 & v);

Geomvec2 multiply_ew (const Geomvec2 & v, const Geomvec2 & w);
Geomvec2 divide_ew (const Geomvec2 & v, const Geomvec2 & w);

std::ostream & operator << (std::ostream & stm, const Geomvec2 & v);
std::istream & operator >> (std::istream & stm, Geomvec2 & v);

Geomvec2 unit_vector_from_angle(const double & angle);

//Geomvec2 random_vec(RNG_taus & rng, double amplitude[]);
Geomvec2 random_vec(RNG_taus & rng, Geomvec2 amplitude=Geomvec2(1.,1.));
Geomvec2 random_normal_vec(RNG_taus & rng, Geomvec2 amplitude=Geomvec2(1.,1.));


#endif


