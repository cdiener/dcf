/*
 * main.cxx
 * This file is part of modes
 *
 * Copyright (C) 2013 - Christian Diener
 *
 * dcf is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * dcf is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cppdesigner. If not, see <http://www.gnu.org/licenses/>.
 */
 
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "iztli.h"
#include "trainer.h"
#include "design.h"

/**
modes program

@section DESCRIPTION

The program takes a text or fasta file as argument followed by the maximum linker size that will be used to join the
sequences in the file and the number of iterations. The probability to be a feasible is calculated by a random forest
classfication based on positive and negative examples of the desired peptides. Those are assumed to be found in files named 
pos_CRIT.txt  and neg_CRIT.txt, where CRIT denotes the classification criteria. Models are automatically cashed, but only used 
if the criteria are passed as arguments.
 
Usage:
To run an optimization joining the sequences found in *sequence_file* with a maximum linker size of *n_linker* using
*n_iter* iterations by using the criteria C1 (and C2) contained in folder data, run the following:
@code{.sh}
> modes sequence_file n_linker n_iter ./data C1 (C2)
@endcode 
*/

const int NTREE = 200; 		// Number of decicion trees in the forest
const int CV = 1; 			// use cross-validation? 0 - use bootstrap, 1 - use cross-validation
const int BEST_FILE = 1; 	// Save best linkers to a file?
const int E_LOG = 1; 		// Save an energy value log?
const int OUT_IT = 100; 	// Output information each OUT_IT iterations

using namespace std;

int main (int argc, char* argv[])
{	
	// Check command line arguments
	if(argc<6)
	{		
		cerr<<"Missing arguments! Need at least 5!"<<endl;
		cout<<"Usage: ./modes sequence_file \'l1 l2 ...\' n_iter data_dir C1 (C2)"<<endl;

		return 1;
	}
	
	string seq_path = argv[1];
	string l_sizes = argv[2];
	int n_iter = atoi(argv[3]);
	string data_dir = argv[4];
	
	double start, end, time_s;
	start = omp_get_wtime();	//set timer
	
	//Classification section
	double mean_err, sd_err, train_err;
	string df_serialized, tmp, est_type, pos_path, neg_path, mod_path;
	alglib::decisionforest df;
	vector<alglib::decisionforest> dfs;
	alglib::real_2d_array data;
	ifstream in_mod_file;
	ofstream out_mod_file;
	trainer rf(NTREE);
	alglib::dfreport rep;
	alglib::ae_int_t info;
	
	for(int i=5; i<argc; i++)
	{
		mod_path = data_dir;
		mod_path += "/model_";
		mod_path += argv[i];
		mod_path += ".txt";
		df_serialized="";
		
		in_mod_file.open(mod_path);
		if(in_mod_file.is_open())
		{
			cout<<"Found saved model for "<<argv[i]<<"."<<endl;
			in_mod_file>>est_type;
			in_mod_file>>mean_err;
			in_mod_file>>sd_err;
			in_mod_file>>train_err;
			
			while(!in_mod_file.eof())
			{
				std::getline(in_mod_file, tmp);
				df_serialized += tmp;
			}
			
			in_mod_file.close();
			df = alglib::decisionforest();
			dfunserialize(df_serialized, df);
			dfs.push_back(df);
		}
		else
		{
			pos_path = data_dir;
			pos_path += "/pos_";
			pos_path += argv[i];
			pos_path += ".txt";
			neg_path = data_dir;
			neg_path += "/neg_";
			neg_path += argv[i];
			neg_path += ".txt";
			data = read_vars( pos_path, neg_path ); 
			if(data.rows() == 0)
			{
				cerr<<"Fatal error: Nothing to classify!"<<endl;
				return 0;
			}
			
			cout<<"Classifying "<<argv[i]<<" on "<<n_var<<" variables over "<<data.rows()<<" peptides";
			cout.flush();
			
			rf = trainer(NTREE);
			if(CV) rf.rep_cv(data, 4, 4);
			else rf.bootstrap(data, 16);
			
			mean_err = rf.mean_error();
			sd_err = rf.sd_error();
			
			// Build complete forest and check for success
			df = alglib::decisionforest();
			dfbuildrandomdecisionforest( data, data.rows(), n_var, 2, NTREE, R, info, df, rep);
			dfs.push_back(df);
			
			if(info==1) cout<<"done."<<endl; else cout<<"failed."<<endl;
			train_err = rep.relclserror;
			
			// Save model
			mod_path = data_dir;
			mod_path += "/model_";
			mod_path += argv[i];
			mod_path += ".txt";
			cout<<"Saved model and error estimates to "<<mod_path<<"."<<endl;
			out_mod_file.open(mod_path);
			est_type = CV?"cv":"bootstrap";
			out_mod_file<<est_type<<'\t';
			out_mod_file<<mean_err<<'\t';
			out_mod_file<<sd_err<<'\t';
			out_mod_file<<train_err<<endl;
			dfserialize(df, df_serialized);
			out_mod_file<<df_serialized;
			out_mod_file.close();
		}
		
		// Some diagnosis
		cout<<"Training set classification error: "<<train_err<<endl;
		cout<<est_type<<" error estimate: "<<mean_err<<" +- "<<sd_err<<endl<<endl;
	}
	
	end = omp_get_wtime();
	time_s = end-start;
	cout<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for classification."<<endl<<endl;

	// Optimization part	
	cout<<"Optimizing on "<<omp_get_max_threads()<<" threads."<<endl;
	start = omp_get_wtime();

	vector<string> seqs = read_seq(seq_path);

	// Parse the linker sizes
	stringstream s_stream(l_sizes);
	std::vector<int> linker_sizes(seqs.size()+1, 0);
	for(unsigned int i=0; i<seqs.size()+1; i++) s_stream >> linker_sizes[i];
	
	sann opt(dfs, seqs, linker_sizes, n_iter, 100, 32);

	cout<<"Initial solution:"<<endl<<opt<<endl;
	
	
	ofstream elog;
	if(E_LOG)
	{
		elog.open("log.txt");
		if(!elog.is_open())
		{	
			cerr<<"Log file could not be opened!"<<endl; 
			return 0;
		}
		elog<<"iter\tsamples\taccepted\tcurrent_E\tbest_E"<<endl; 
		elog<<"0\t0\t0\t"<<opt.current_energy()<<'\t'<<opt.get_best_energy()<<endl;
	}
	
	#ifdef CURSES_HAVE_CURSES_H
	int curses_success = opt.init_curses();
	if(!curses_success)
	{ 
		opt.end_curses();
		cout<<"Terminal window too small for nice output - using fallback..."<<endl; 
		cout.flush();
	} 
	
	for(int i=0; i<n_iter; i++)
	{
		opt.anneal();
		
		if( (i+1)%OUT_IT==0 )
		{ 
			if(curses_success) opt.update_curses();
			else
			{
				cout<<"        \r";
				cout<<opt;
				cout.flush();
			}
				
		}
		
		if(E_LOG)
			elog<<i+1<<'\t'<<(i+1)*100<<'\t'<<opt.get_accepted()<<'\t'<<opt.current_energy()<<'\t'<<opt.get_best_energy()<<endl;
	}
	if(curses_success) opt.end_curses();
	#else
	for(unsigned int i=0; i<n_iter; i++)
	{
		opt.anneal();
		
		if( (i+1)%OUT_IT==0 )
		{ 
			
			cout<<"        \r";
			cout<<opt;
			
			//cout<<opt.get_state()<<endl<<endl;
			cout.flush();
		}
		
		if(E_LOG)
			elog<<i+1<<'\t'<<(i+1)*100<<'\t'<<opt.get_accepted()<<'\t'<<opt.current_energy()<<'\t'<<opt.get_best_energy()<<endl;
	}
	#endif
	
	
	if(E_LOG) elog.close();

	if(BEST_FILE)
	{
		ofstream out("best_link.txt");
		if(out.is_open()) out <<opt.get_best();
		out.close();
		cout<<endl<<endl;
	}
	
	cout<<"\n\nBest sequences found: "<<endl<<opt.get_best(0.0, 0)<<endl;
	
	end = omp_get_wtime();	
	time_s = end-start;
	cout<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for optimization."<<endl;

	return 0;
}

  
 
