#include "traits.h"

const unsigned int Traits<double>::dim;


double norm(const double & d)
	{
	return std::abs(d);
	}

double norm2(const double & d)
	{
	return d*d;
	}

