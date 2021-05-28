#ifndef util_h
#define util_h
#include <iostream>
//#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/types.h>
#include <iterator>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
using namespace std;

#define BOOST 1
#define TRILINOS 2
//#define SOLVER_TYPE TRILINOS // BOOST
#define NONE 1
#define GNUPLOT 2
//#define PLOTTER_TYPE NONE // GNUPLOT

//-------------------------------------------------------------

//const int nd=2;
const double pi(4.0*atan2(1.0, 1.0));
const double dblpi(8.0*atan2(1.0, 1.0));
const std::string sep="    ";
const mode_t directory_permissions(0775);

//-------------------------------------------------------------

// convert any type to a string via a stringstream
template < class T > inline std::string to_string (const T & t)
	{
	std::stringstream ss;
	ss << t;
	return ss.str();
	};
	
// extend the definition to std::vector with the separator as an optional argument
template < class T > inline std::string to_string (const std::vector <T> & V, std::string s = sep)
	{
	std::stringstream ss;
	for (size_t i=0; i<V.size(); i++) ss << V[i] << s;
	return ss.str();
	};

// print a std::vector on a single line using the << operator
template <class T> std::ostream & operator << (std::ostream & os, const std::vector <T> & V)
	{
	for (size_t i=0; i<V.size(); i++) os << V[i] << sep;
	return os;
	};

// print a std::vector on a single line with default separator = sep
template <class T> std::ostream & print (const std::vector <T> & V, std::ostream & os, std::string s = sep)
	{
	for (size_t i=0; i<V.size(); i++) os << V[i] << s;
	return os;
	};


// find the bin a value is in in an ordered std::vector containing bin edges
template <class T> unsigned int find_bin (const std::vector <T> & V, const T & v)
	{
	if ( v<*V.begin() || v>*V.rbegin() ) return V.size()-1;
	unsigned int i1,i2,i;
	i1=0;
	i2=V.size()-1;
	while (i2-i1>1)
		{
		i=i1+(i2-i1)/2;
		if (v<V[i]) i2=i;
		else i1=i;
		}
	return i1;
	};

// Defines a string with holes to fill in by calling the () operator.
class Namer
{
public:
	string s1, s2, s3;
	
	Namer(string s1, string s2 = "", string s3 = "") : s1(s1), s2(s2), s3(s3) {};
	template <class T1, class T2> string operator () (T1 & t1, T2 & t2 = "")
		{
		return s1 + to_string(t1) + s2 + to_string(t2) + s3;
		}
	// default 2nd argument above doesn't work (why?) -> overload
	template <class T1> string operator () (T1 & t1)
		{
		return s1 + to_string(t1) + s2;
		}
};

inline double & bring_back_angle(double & angle)
	{
	while (angle>pi) { angle-=dblpi; }
	while (angle<=-pi) { angle+=dblpi; }
	return angle;
	}

inline void wait() { std::cin.ignore(); }

inline void read_until_number(std::istream & is)
	{
	bool done = false;
	char a;
	while ( !done )
		{
		a = is.get();
		if ( isdigit(a) || a=='-' || a=='.' ) done = true;
		}
	is.putback(a);
	}

inline std::string TimeToDate(double s)
	{
	int m=int(s/60); s-=60*m;
	int h=int(m/60); m-=60*h;
	int d=int(h/24); h-=24*d;
	stringstream ss;
	if (d>0) { ss << d << "days "; }
	if (h>0) { ss << h << "h "; }
	if (m>0) { ss << m << "min "; }
	ss << s << "s";
	return ss.str();
	};

inline std::string get_date()
	{
	time_t t;
	time(&t); 
	std::string s = asctime(localtime(&t));
	return s.substr(0,s.length()-1);
	};


//-------------------------------------------------------------

#include "cxxabi.h"  // provides demangler for typeid.name()

inline const std::string demangle(const char* name)
	{
	int status = -4;
	char* res = abi::__cxa_demangle(name, NULL, NULL, &status);
	const char* const demangled_name = (status==0)?res:name;
	std::string ret_val(demangled_name);
	free(res);
	return ret_val;
	}

template <class T> std::string demangled_type(T &)
	{
	return demangle(typeid(T).name());
	}

template <class T> std::string demangled_type()
	{
	return demangle(typeid(T).name());
	}

//------------------------------------------------------------------

inline void copy_overwrite(boost::filesystem::path from, boost::filesystem::path to)
	{
	if ( boost::filesystem::exists(to) ) boost::filesystem::remove(to);
	copy(from,to);
	}

//------------------------------------------------------------------

#endif
