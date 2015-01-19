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

The program takes a text or fasta file as argument followed by the maximum linker size that will be used to join the
sequences in the file and the number of iterations. The probability to be a feasible is calculated by a random forest
classfication based on positive and negative examples of the desired peptides. Those are assumed to be found in files named 
pos_CRIT.txt  and neg_CRIT.txt, where CRIT denotes the classification criteria. Models are automatically cashed, but only used 
if the criteria are passed as arguments.
 
Usage:
To run a prediction of sequences found in *sequence_file*, use the following:
@code{.sh}
> predict sequence_file n_iter ./data C1 (C2)
@endcode 
*/

const int NTREE = 50; 		// Number of decicion trees in the forest
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
		cout<<"Usage: ./predict sequence_file data_dir C1 (C2)"<<endl;

		return 1;
	}
	
	string seq_path = argv[1];
	string data_dir = argv[2];
	
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
	
	for(int i=3; i<argc; i++)
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

	// Prediction part	
	start = omp_get_wtime();
	
	vector<string> seqs = read_seq(seq_path);
	vector<double> all(dfs.size()+1, 0.0);
	vector<double> current_seq;
	
	// Predict
	#pragma omp parallel for
	for(int i=0; i<seqs.size(); i++)
	{
		current_seq = predict(dfs, seqs[i]);
		for(int j=0; j<current_seq.size(); j++) all[j] += current_seq[j]; 
	}
	
	end = omp_get_wtime();	
	time_s = end-start;
	cout<<"Needed "<<time_s<<" s (+- "<<omp_get_wtick()<<" s) for prediction."<<endl;
	cout<<"Predicted probabilities"<<endl;
	for(int i=0; i<all.size(); i++) cout<<all[i]/seqs.size()<<'\t';
	cout<<endl;

	return 0;
}

  
 
