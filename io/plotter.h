#include <fstream>
#include <string>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"
#include "geomvec/geomvec2.h"
#include "geomvec/geomvec3.h"
#include "parameters/parameters.h"

class Plotter_dumb
{
public:
/*	Plotter_dumb() {};*/
	bool is_dumb() { return true; };
	template <class T>
	void plot (std::vector <T> & u, std::string title ) {};
};

//==========================================================

// [Warning] Gnuplot::send accepts vectors of boost::tuple
// or std::pair, but not std::tuple.

class Plotter : public Gnuplot
{
public:
	std::vector <double> l;
	
	Plotter(Parameters & par);
	virtual ~Plotter() {};
	bool is_dumb() { return false; };
	template <class T>
	void plot (std::vector <T> & u, std::string title ) {};
};


Plotter::Plotter(Parameters & par) : Gnuplot(), l(par.N)
	{
	std::ifstream is(par.gnuplot_preamble_file);
	std::string line;
	while ( getline(is,line) )
		{
//		std::clog << line << std::endl;
		(*this) << line << std::endl;
		}
	for (size_t i=0; i<l.size(); i++) l[i] = par.dl * i;
	};

template <>
void Plotter::plot <double> (std::vector <double> & u, std::string title)
	{
	(*this) << "set title '"+title+"'\n";
	(*this) << "plot '-' @style1, 0 lt 0\n";
	
	typedef boost::tuple <double,double> Point;
	std::vector <Point> v(l.size());
	for (size_t i=0; i<l.size(); i++) v[i] = Point(l[i],u[i]);
	send(v);
	};

template <>
void Plotter::plot <Geomvec2> (std::vector <Geomvec2> & u, std::string title)
	{
	(*this) << "set title '"+title+"'\n";
	(*this) << "plot '-' @style1, 0 lt 0\n";
	
//	clog << u << endl;
	
	typedef boost::tuple <double,double> Point;
	std::vector <Point> v(l.size());
	for (size_t i=0; i<l.size(); i++) v[i] = Point(u[i].x,u[i].y);
	send(v);
//	exit(3);
	};

template <>
void Plotter::plot <Geomvec3> (std::vector <Geomvec3> & u, std::string title)
	{
	(*this) << "set title '"+title+"'\n";
	(*this) << "splot '-' @style1\n";
	
	typedef boost::tuple <double,double,double> Point;
	std::vector <Point> v(l.size());
	for (size_t i=0; i<l.size(); i++) v[i] = Point(u[i].x,u[i].y,u[i].z);
	send(v);
	};


