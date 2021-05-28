#ifndef traits_h
#define traits_h
#include <cmath>

// The class 'Solver' needs to access the dimension of its 'Point_type'.
// Geomvec2 and Geomvec3 contain a 'static const dim' with the appropriate value,
// but double doesn't, hence the use of a trait strategy.

template <class T> class Traits { public: static const unsigned int dim = T::dim; };
template <> class Traits <double> { public: static const unsigned int dim = 1u; };

template <class T> const unsigned int Traits<T>::dim;

// provided for compatibility with vector types
double norm(const double & d);
double norm2(const double & d);

#endif


