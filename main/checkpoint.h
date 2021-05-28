#ifndef checkpoint_h
#define checkpoint_h
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include "main/util.h"
#include "system/system.h"
#include "parameters/parameters.h"
using namespace std;
namespace bfs = boost::filesystem;

class Checkpoint
{
public:
	bool used;
	string checkpoint_dir;
	
	template <class TSystem> Checkpoint(TSystem & sys, Parameters & par);
	template <class TSystem> void save_checkpoint(const TSystem & sys);
	template <class TSystem> void save_checkpoint(const TSystem & sys, ostream & stm);
	template <class TSystem> void load_checkpoint(TSystem & sys);
	template <class TSystem> void load_checkpoint(TSystem & sys, istream & stm);
	string find_last_checkpoint();
	bool is_valid_checkpoint(istream & stm);
	bool is_valid_checkpoint(const bfs::path & p);
};

#include "checkpoint.hpp"

#endif
