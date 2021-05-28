
template <class TSystem>
Checkpoint::Checkpoint(TSystem & sys, Parameters & par) : used(false), checkpoint_dir(par.checkpoint_dir)
	{
	mkdir(checkpoint_dir.c_str(),directory_permissions);
	string checkpoint_file = find_last_checkpoint();
	if ( par.use_checkpoint && checkpoint_file!="" ) 
		{
		cout << "[Warning] Loading checkpoint file \""+checkpoint_file+"\". " ;
		ifstream stm_chk(checkpoint_file.c_str());
		load_checkpoint(sys,stm_chk);
		cout << "Resuming at iteration " << sys.step << "/" << par.n_step << "." << endl;
		
		// Read 'max_displacement.dat' up to the last checkpoint.
		vector <string> filenames = {"max_displacement.dat","time-headPos-headAng-curvEnergy.dat"};
		for (auto f: filenames)
			{
			string filename = par.outdir+"/"+f;
			string file_contents("");
			fstream stm(filename.c_str(),fstream::in);
			string line;
			while ( getline(stm,line) )
				{
				istringstream iss(line); double t; iss >> t;
//				if ( line[0] == '#' || ( t<sys.t || iss.fail() ) ) break;
//				if ( line[0] != '#' && iss.fail() ) break;
				file_contents += line + "\n";
				if ( line[0] != '#' && ( t>sys.t || iss.fail() ) ) break;
				}
			stm.close();
			// Clear [filename] and rewrite contents up to last checkpoint.
			stm.open(filename.c_str(),fstream::out | fstream::trunc);
			stm << file_contents << flush;
			stm.close();
			}
		used = true;
		}
	else
		{
		if ( par.use_checkpoint ) cout << "[Warning] No available checkpoint file; starting over." << endl;
		}		
	};

template <class TSystem>
void Checkpoint::save_checkpoint(const TSystem & sys)
	{
	string checkpoint_file = checkpoint_dir+"/checkpoint";
	// Rewrite over checkpoints. Only keep current and previous.
	if ( bfs::exists(checkpoint_file) ) bfs::rename(checkpoint_file,checkpoint_file+"_-1");
	ofstream stm(checkpoint_file.c_str());
	save_checkpoint(sys,stm);
	stm.close();
	}
	
template <class TSystem>
void Checkpoint::save_checkpoint(const TSystem & sys, ostream & stm)
	{
	sys.print_full(stm);
	stm << endl << "END";
	}

template <class TSystem>
void Checkpoint::load_checkpoint(TSystem & sys)
	{
	string checkpoint_file = checkpoint_dir+"/checkpoint"+to_string(sys.step);
	ifstream stm(checkpoint_file.c_str());
	load_checkpoint(sys,stm);
	stm.close();
	}

template <class TSystem>
void Checkpoint::load_checkpoint(TSystem & sys, istream & stm)
	{
	sys.read_full(stm);
	}

//------------------------------------------------------------

string Checkpoint::find_last_checkpoint()
	{
	bfs::path p(checkpoint_dir);
	if ( ! bfs::is_directory(p) ) return "";
	
	bfs::directory_iterator dir_first(p), dir_last;
	std::vector <bfs::path> files;
	
	// retrieve the list of valid checkpoints
	auto pred = [this](const bfs::directory_entry& p)
		{ return this->is_valid_checkpoint(bfs::path(p)); };
	std::copy(boost::make_filter_iterator(pred, dir_first, dir_last),
		boost::make_filter_iterator(pred, dir_last, dir_last),
		std::back_inserter(files)
		);
	// find the most recent one
	auto newest = std::max_element ( files.begin(), files.end(),
	                                [](const bfs::path& p1, const bfs::path& p2)
		                             { return bfs::last_write_time(p1) < bfs::last_write_time(p2); } );
	
	if ( newest == files.end() ) return "";
	else return newest->string();
	}

// check that a checkpoint was written entirely (i.e. last line is "END")
bool Checkpoint::is_valid_checkpoint(istream & stm)
	{
	string line;
	while ( getline(stm,line) )
		{
		if (line=="END") return true;
		}
	return false;
	}

bool Checkpoint::is_valid_checkpoint(const bfs::path & p)
	{
	bfs::ifstream stm(p);
	return bfs::is_regular_file(p) && is_valid_checkpoint(stm);
	}

