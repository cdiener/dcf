/*
 * main.cxx
 * This file is part of cppdesigner
 *
 * Copyright (C) 2013 - Christian Diener
 *
 * cppdesigner is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * cppdesigner is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cppdesigner. If not, see <http://www.gnu.org/licenses/>.
 */
 
#include <iostream>
#include <fstream>
#include <omp.h>
#include "iztli.h"
#include "design.h"

/**
cppdesigner program

@section DESCRIPTION

The program takes a text or fasta file as argument followed by the maximum linker size that will be used to join the
sequences in the file to cell penetrating fusion peptide. The probability to be a CPP is calculated by a random forest
classfication based on positive and negative examples of CPPs. Those are assumed to be found in the files pos.txt and
neg.txt. Or, if not, must be passed as two additional arguments. 
*/

const int NTREE = 64;
const char* MOD_PATH = "model.txt";
const int CV = 0; // use cross-validation? 0 - use bootstrap, 1 - use cross-validation
const int BEST_FILE = 0; // Save best linkers to a file?
const int E_LOG = 0; // Save an energy value log?

using namespace std;

int main (int argc, char* argv[])
{	
	// Check command line arguments
	if(argc<4)
	{		
		cerr<<"Missing arguments! Need at least 3!"<<endl;
		cout<<"Usage: ./cppdesigner sequence_file n_linker n_iter (pos_CPPs neg_CPPs) (high_CPPs lowCPPs)"<<endl;

		return 1;
	}

	int n_linker, n_iter;
	string seq_path, pos_path, neg_path, high_path, low_path;
	pos_path = "examples/pos.txt";
	neg_path = "examples/neg.txt";
	high_path = "examples/high.txt";
	low_path = "examples/low.txt";

	switch(argc)
	{
		case 4:	seq_path = argv[1];
				n_linker = atoi(argv[2]);
				n_iter = atoi(argv[3]);
				break;
		case 6: seq_path = argv[1];
				n_linker = atoi(argv[2]);
				n_iter = atoi(argv[3]);
				pos_path = argv[4];
				neg_path = argv[5];
				break;
		case 8: seq_path = argv[1];
				n_linker = atoi(argv[2]);
				n_iter = atoi(argv[3]);
				pos_path = argv[4];
				neg_path = argv[5];
				high_path = argv[6];
				low_path = argv[7];
				break;
		default: cerr<<"Wrong number of arguments! Need 3, 5 or 7!"<<endl;
				 cout<<"Usage: ./cppdesigner sequence_file n_linker (pos_CPPs neg_CPPs) (high_CPPs lowCPPs)"<<endl;
				 return 1;
			
	}

	double start, end;
	double time_s;
	//clock_gettime(CLOCK_MONOTONIC, &start); 
	
	//Classification section
	double mean_err, sd_err, train_err;
	std::string df_serialized, tmp, est_type;
	alglib::decisionforest df;
	
	ifstream mod_file(MOD_PATH);
	if(mod_file.is_open())
	{
		cout<<"Found saved model."<<endl;
		mod_file>>est_type;
		mod_file>>mean_err;
		mod_file>>sd_err;
		mod_file>>train_err;
		
		while(!mod_file.eof())
		{
			std::getline(mod_file, tmp);
			df_serialized += tmp;
		}
		
		mod_file.close();
		dfunserialize(df_serialized, df);
	}
	else
	{
		start = omp_get_wtime();	//set timer
		
		alglib::real_2d_array data = read_vars( pos_path, neg_path ); 
		if(data.rows() == 0)
		{
			cerr<<"Fatal error: Nothing to classify!"<<endl;
			return 0;
		}
		
		cout<<"Classifying on "<<n_var<<" variables...";
		cout.flush();
		
		trainer rf(NTREE);
		if(CV) rf.rep_cv(data, 8, 8);
		else rf.bootstrap(data, 64);
		
		mean_err = rf.mean_error();
		sd_err = rf.sd_error();

		alglib::dfreport rep;
		alglib::ae_int_t info;
		
		// Build complete forest and check for success
		dfbuildrandomdecisionforest( data, data.rows(), n_var, 2, NTREE, R, info, df, rep);
		if(info==1) cout<<"done."<<endl; else cout<<"failed."<<endl;
		train_err = rep.relclserror;
		
		// Save model
		cout<<"Saved model and error estimates to "<<MOD_PATH<<"."<<endl;
		ofstream mod_file(MOD_PATH);
		est_type = CV?"cv":"bootstrap";
		mod_file<<est_type<<'\t';
		mod_file<<mean_err<<'\t';
		mod_file<<sd_err<<'\t';
		mod_file<<train_err<<endl;
		dfserialize(df, df_serialized);
		mod_file<<df_serialized;
		mod_file.close();
		
		//clock_gettime(CLOCK_MONOTONIC, &end);	
		end = omp_get_wtime();
		
		//time_s = end.tv_sec-start.tv_sec + 1.0e-9*end.tv_nsec - 1.0e-9*start.tv_nsec;
		time_s = end-start;
		cout<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for classification."<<endl<<endl;
	}
	
	// Some diagnosis
	cout<<"Training set classification error: "<<train_err<<endl;
	cout<<est_type<<" error estimate: "<<mean_err<<" +- "<<sd_err<<endl<<endl;

	// Optimization part	
	cout<<"Optimizing on "<<omp_get_max_threads()<<" threads."<<endl;
	//clock_gettime(CLOCK_MONOTONIC, &start);
	start = omp_get_wtime();

	vector<string> seqs = read_seq(seq_path);
	sann opt(df, seqs, n_iter, n_linker, 100, NTREE);
	cout<<"Initial solution:"<<endl<<opt<<endl;
	
	
	ofstream elog;
	if(E_LOG)
	{
		elog.open("log.txt");
		if(!elog.is_open())
		{	
			cout<<"Log file could not be opened!"<<endl; 
			return 0;
		}
		elog<<"iter\tsamples\taccepted\tcurrent_E\tbest_E"<<endl; 
		elog<<"0\t0\t0\t"<<opt.current_energy()<<'\t'<<opt.get_best_energy()<<endl;
	}
	
	for(unsigned int i=0; i<n_iter; i++)
	{
		opt.anneal();
		
		if( (i+1)%100==0 )
		{ 
			cout<<"        \r";
			cout<<opt;
			cout.flush();
		}
		
		if(E_LOG)
			elog<<i+1<<'\t'<<(i+1)*100<<'\t'<<opt.get_accepted()<<'\t'<<opt.current_energy()<<'\t'<<opt.get_best_energy()<<endl;
	}
	
	if(E_LOG) elog.close();

	if(BEST_FILE)
	{
		ofstream out("best_link.txt");
		if(out.is_open()) out <<opt.get_best();
		out.close();
		cout<<endl<<endl;
	}
	else cout<<"\n\nBest sequences found: "<<endl<<opt.get_best()<<endl;
	
	end = omp_get_wtime();	
	//time_s = end.tv_sec-start.tv_sec + 1.0e-9*end.tv_nsec - 1.0e-9*start.tv_nsec;
	time_s = end-start;
	cout<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for optimization."<<endl;

	return 0;
}

  
 
