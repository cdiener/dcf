/*
 * trainer.cxx
 * 
 * Copyright 2014 Christian Diener <ch.diener@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include "trainer.h"
#include "iztli.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>

// Classification functions

alglib::real_2d_array read_vars(const std::string pos, const std::string neg)
{
	std::vector<std::string> pos_seq, neg_seq;
	std::vector<int> seq, counts;
	alglib::real_2d_array full;

	//Read sequences
	pos_seq = read_seq(pos);
	neg_seq = read_seq(neg);
	
	if(pos_seq.size() == 0 || neg_seq.size() == 0) 
	{
		std::cerr<<"No sequences loaded!"<<std::endl;
		throw 0;
	}
	
	full.setlength(pos_seq.size()+neg_seq.size(), n_var+1);
	
	for(unsigned int i=0; i<pos_seq.size(); i++)
	{
		// Calculate iztli parameters
		seq = string_to_array(pos_seq[i]);
		full[i][0] = charge(seq);
		full[i][1] = Hm(seq);
		full[i][2] = pI(seq);
		full[i][3] = rangeHM(seq);
		full[i][4] = charge_var(seq);
		counts = countAA(seq);
		full[i][5] = logP(counts);
		full[i][6] = alpha(seq);
		for(unsigned int j=7; j<27; j++) full[i][j] = 1.0*counts[j-7]/seq.size();


		// Mark as positive class
		full[i][n_var] = 0.0;
	}

	for(unsigned int i=pos_seq.size(); i<pos_seq.size()+neg_seq.size(); i++)
	{
		// Calculate iztli parameters
		seq = string_to_array( neg_seq[i-pos_seq.size()] );
		full[i][0] = charge(seq);
		full[i][1] = Hm(seq);
		full[i][2] = pI(seq);
		full[i][3] = rangeHM(seq);
		full[i][4] = charge_var(seq);
		counts = countAA(seq);
		full[i][5] = logP(counts);
		full[i][6] = alpha(seq);
		for(unsigned int j=7; j<27; j++) full[i][j] = 1.0*counts[j-7]/seq.size();


		// Mark as negative class
		full[i][n_var] = 1.0;
	}

	return full;
} 

class_set::class_set()
{
	srand( 1e6*omp_get_wtime() + alglib::randominteger(RAND_MAX/2)+omp_get_thread_num() );
}

void class_set::generate(const alglib::real_2d_array data, double train_frac)
{
	std::vector<int> idx( data.rows() );
	for(unsigned int i=0; i<data.rows(); i++) idx[i] = i;
  	
	std::random_shuffle( idx.begin(), idx.end() );

	int n_train = std::round( train_frac*data.rows() );
	this->train.setlength( n_train, data.cols() );
	this->test.setlength( data.rows()-n_train, data.cols() );

	for(int i=0; i<n_train; i++)
		for(unsigned int j=0; j<data.cols(); j++) this->train[i][j] = data[ idx[i] ][j];

	for(unsigned int i=0; i<data.rows()-n_train; i++)
		for(unsigned int j=0; j<data.cols(); j++) this->test[i][j] = data[ idx[i+n_train] ][j];

	this->n_train = n_train;
	this->n_test = data.rows()-n_train;
}

void class_set::generate_bootstrap(const alglib::real_2d_array data, int n)
{
	std::vector<int> idx( data.rows(), 0 );
	this->train.setlength( n, data.cols() );
	int not_used=0;
	
	// We want to unsure that the test set has at least 3 entries
	while(not_used < 3)
	{
		for(unsigned int i=0; i<data.rows(); i++) idx[alglib::randominteger(idx.size())]++;
		not_used = 0;
		for(unsigned int i=0; i<data.rows(); i++) if(idx[i]==0) not_used++;
		if(not_used<3) continue;
	
		this->test.setlength( not_used, data.cols() );
		
		int train_i = 0, test_i = 0;
		for(unsigned int i=0; i<data.rows(); i++)
		{
			if( idx[i]>0 )
			{
				for(int j=0; j<idx[i]; j++) 
					for(unsigned int k=0; k<data.cols(); k++) this->train[train_i+j][k] = data[i][k];
				train_i += idx[i];
			}
			else
			{
				for(unsigned int k=0; k<data.cols(); k++) test[test_i][k] = data[i][k];
				test_i++;
			}
		}
	}
		
}

trainer::trainer(int n_tree)
{
	this->n_tree = n_tree;
	// Innitialize random number seeds (different for each thread)
	#pragma omp parallel
	{
		srand( 1e6*omp_get_wtime() + alglib::randominteger(RAND_MAX/2)+omp_get_thread_num() );
	}
}

void trainer::bootstrap(const alglib::real_2d_array data, int rep)
{
	if(reset)
	{
		errors.resize(0);
		models.resize(0);
	}
	
	#pragma omp parallel
	{
		alglib::decisionforest df;
		alglib::dfreport dfrep;
		alglib::ae_int_t info;
		class_set boot;
		double err;
		
		#pragma omp for
		for(int i=0; i<rep; i++)
		{
			boot.generate_bootstrap(data, (int)data.rows());
			dfbuildrandomdecisionforest( boot.train, boot.train.rows(), n_var, 2, n_tree, R, info, df, dfrep);
			err = 0.632*dfrelclserror(df, boot.test, boot.test.rows()) + 0.368*dfrep.relclserror;
			#pragma omp critical 
			{
				errors.push_back(err);
				models.push_back(df);
			}
			
			std::cout<<".";
			std::cout.flush();
		}
	}
}

void trainer::cv(const alglib::real_2d_array data, int folds)
{
	srand( 1e6*omp_get_wtime() + alglib::randominteger(RAND_MAX/2)+omp_get_thread_num() );
	
	if(reset)
	{
		errors.resize(0);
		models.resize(0);
	}
		
	std::vector<int> idx( data.rows() );
	for(unsigned int i=0; i<data.rows(); i++) idx[i] = i;
  	std::random_shuffle( idx.begin(), idx.end() );
	
	#pragma omp parallel
	{
		alglib::decisionforest df;
		alglib::dfreport dfrep;
		alglib::ae_int_t info;
		alglib::real_2d_array test, train;
		std::vector<int> tmp_idx;
		double err;
		int start, end;
		
		#pragma omp for
		for(int i=0; i<folds; i++)
		{
			tmp_idx = idx;
			// Decide for the test set
			start = i*floor(data.rows()/folds);
			if(i<folds-1) end = (i+1)*floor(data.rows()/folds);
			else end = data.rows();
			test.setlength(end-start, data.cols());

			for(unsigned int k=0; k<test.rows(); k++)
			{
				for(unsigned int j=0; j<data.cols(); j++) test[k][j] = data[ idx[start+k] ][j];
			}
			tmp_idx.erase(tmp_idx.begin()+start, tmp_idx.begin()+end);
				
			// Generate the training set
			train.setlength(tmp_idx.size(), data.cols());
			
			for(unsigned int k=0; k<tmp_idx.size(); k++)
				for(unsigned int j=0; j<data.cols(); j++) train[k][j] = data[ tmp_idx[k] ][j]; 
			
			dfbuildrandomdecisionforest( train, train.rows(), n_var, 2, n_tree, R, info, df, dfrep);
			err = dfrelclserror(df, test, test.rows());
			
			#pragma omp critical 
			{
				errors.push_back(err);
				models.push_back(df);
			}
			
			std::cout<<".";
			std::cout.flush();
		}
	}
}

void trainer::rep_cv(const alglib::real_2d_array data, int rep, int folds)
{
	if(reset)
	{
		errors.resize(0);
		models.resize(0);
		reset = 0;
	}
	
	for(int i=0; i<rep; i++) this->cv(data, folds);
	reset = 1;
}

double trainer::mean_error()
{
	double me=0;
	
	for(unsigned int i=0; i<errors.size(); i++) me += errors[i];
	
	return me/errors.size();
}

double trainer::sd_error()
{
	double sde=0, me = mean_error();
	if(errors.size()<2) return -1;
	
	for(unsigned int i=0; i<errors.size(); i++) sde += (errors[i] - me)*(errors[i] - me);
	
	return sqrt(sde)/(errors.size()-1);
}

double trainer::min_error()
{
	double mine=HUGE_VAL;
	
	for(unsigned int i=0; i<errors.size(); i++) if(errors[i] < mine) mine = errors[i];
	
	return mine;
}

alglib::decisionforest trainer::best_model() const
{
	double mine = HUGE_VAL;
	alglib::decisionforest best = models[0];
	for(unsigned int i=1; i<errors.size(); i++) if(errors[i]<mine) best = models[i];
	
	return best;
}
