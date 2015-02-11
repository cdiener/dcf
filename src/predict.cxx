/*
 * predict.cxx
 * This file is part of modes
 *
 * Copyright (C) 2014 - Christian Diener
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
#include <sstream>
#include <omp.h>
#include "iztli.h"
#include "trainer.h"
#include "design.h"

/**
modes program

@section DESCRIPTION

Predicts the avergae probability of the criteria C1..CN for a list of sequence
files. It will always report the individual probabilities P(A), P(B), etc, as 
well as the joint probability P(A & B & ...). 
 
Usage:
To run a prediction of sequences found in *sequence_file*, use the following:
@code{.sh}
> ./predict data_dir 'C1 (C2) ...' data_file1 data_file2 ...
@endcode 
*/

const int NTREE = 64; 		// Number of decision trees in the forest
const int CV = 1; 			// use cross-validation? 0 - use bootstrap, 1 - use cross-validation

using namespace std;

/**
 * Predicts the probability of a sequence to constitute any of the tested classes.
 * 
 * */
vector<double> predict(vector<alglib::decisionforest> dfs, const string& seq)
{
	vector<int> i_seq = string_to_array(seq);
	vector<int> counts;
	vector<double> preds(dfs.size()+1);
	preds[dfs.size()] = 1.0;
	alglib::real_1d_array vars;
	alglib::real_1d_array probs;
	vars.setlength(n_var);
	
	vars[0] = charge(i_seq);
	vars[1] = Hm(i_seq);
	vars[2] = pI(i_seq);
	vars[3] = rangeHM(i_seq);
	vars[4] = charge_var(i_seq);
	counts = countAA(i_seq);
	vars[5] = logP(counts);
	vars[6] = alpha(i_seq);
	for(unsigned int j=n_var-20; j<n_var; j++) vars[j] = 1.0*counts[j-n_var+20]/i_seq.size();
	
	for(unsigned int i=0; i<dfs.size(); i++)
	{
		dfprocess(dfs[i], vars, probs);
		preds[i] = probs[0];
		preds[dfs.size()] *= probs[0];
	}
	
	return preds;
}

int main (int argc, char* argv[])
{	
	// Check command line arguments
	if(argc<4)
	{		
		cerr<<"Missing arguments! Need at least 3!"<<endl;
		cout<<"Usage: ./predict data_dir 'C1 (C2) ...' data_file1 data_file2 ..."<<endl;

		return 1;
	}
	
	ofstream log_file("predict_log.txt");

	string data_dir = argv[1];
	string crit_list = argv[2];
	
	// Parse the criteria
	stringstream s_stream(crit_list);
	vector<string> criteria;
	while(s_stream){
		string c;
		s_stream >> c;
		if(c=="") break;
		criteria.push_back(c);
	}
	
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
	
	for(int i=0; i<criteria.size(); i++)
	{
		mod_path = data_dir;
		mod_path += "/model_";
		mod_path += criteria[i];
		mod_path += ".txt";
		df_serialized="";
		
		in_mod_file.open(mod_path);
		if(in_mod_file.is_open())
		{
			log_file<<"Found saved model for "<<criteria[i]<<"."<<endl;
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
			pos_path += criteria[i];
			pos_path += ".txt";
			neg_path = data_dir;
			neg_path += "/neg_";
			neg_path += criteria[i];
			neg_path += ".txt";
			data = read_vars( pos_path, neg_path ); 
			if(data.rows() == 0)
			{
				cerr<<"Fatal error: Nothing to classify!"<<endl;
				return 0;
			}
			
			log_file<<"Classifying "<<criteria[i]<<" on "<<n_var<<" variables over "<<data.rows()<<" peptides";
			
			rf = trainer(NTREE);
			if(CV) rf.rep_cv(data, 4, 4);
			else rf.bootstrap(data, 16);
			
			mean_err = rf.mean_error();
			sd_err = rf.sd_error();
			
			// Build complete forest and check for success
			df = alglib::decisionforest();
			dfbuildrandomdecisionforest( data, data.rows(), n_var, 2, NTREE, R, info, df, rep);
			dfs.push_back(df);
			
			if(info==1) log_file<<"done."<<endl; else log_file<<"failed."<<endl;
			train_err = rep.relclserror;
			
			// Save model
			mod_path = data_dir;
			mod_path += "/model_";
			mod_path += criteria[i];
			mod_path += ".txt";
			log_file<<"Saved model and error estimates to "<<mod_path<<"."<<endl;
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
		log_file<<"Training set classification error: "<<train_err<<endl;
		log_file<<est_type<<" error estimate: "<<mean_err<<" +- "<<sd_err<<endl<<endl;
	}
	
	end = omp_get_wtime();
	time_s = end-start;
	log_file<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for classification."<<endl<<endl;

	// Prediction part	
	start = omp_get_wtime();
	
	#pragma omp parallel for
	for(int a_idx=3; a_idx<argc; a_idx++)
	{ 
		vector<string> seqs = read_seq(argv[a_idx]);
		vector<double> all(dfs.size()+1, 0.0);
		vector<double> probs;
		
		//Get family name
		string fam_name(argv[a_idx]);
		int name_start = fam_name.find('/');
		int name_end = fam_name.find('.');
		fam_name = fam_name.substr(name_start+1, name_end-name_start-1); 
		
		// Predict
		for(int i=0; i<seqs.size(); i++)
		{
			probs = predict(dfs, seqs[i]);
			
			for(int j=0; j<probs.size(); j++) all[j] += probs[j]; 
		}
		
		#pragma omp critical
		{
			cout<<fam_name;
			for(int j=0; j<all.size(); j++) cout<<'\t'<<all[j]/seqs.size();
			cout<<endl;
		}
	}
	end = omp_get_wtime();	
	time_s = end-start;
	log_file<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for prediction.";

	return 0;
}

  
 
