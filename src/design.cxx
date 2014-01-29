/*
 * design.cxx
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

#include "design.h"
#include "iztli.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

// Classification functions

alglib::real_2d_array read_vars(std::string pos, std::string neg)
{
	std::vector<std::string> pos_seq, neg_seq;
	double collector;
	int *seq, *counts;
	int j;
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
		full[i][3] = maxHM(seq);
		full[i][4] = charge_var(seq);
		counts = countAA(seq);
		full[i][5] = logP(seq, counts);
		for(unsigned int j=6; j<26; j++) full[i][j] = counts[j-6];

		delete seq, counts;

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
		full[i][3] = maxHM(seq);
		full[i][4] = charge_var(seq);
		counts = countAA(seq);
		full[i][5] = logP(seq, counts);
		for(unsigned int j=6; j<26; j++) full[i][j] = counts[j-6];

		delete seq, counts;

		// Mark as negative class
		full[i][n_var] = 1.0;
	}

	return full;
} 

class_set::class_set()
{
	srand( 1e6*omp_get_wtime() + alglib::randominteger(RAND_MAX/2)+omp_get_thread_num() );
}

void class_set::generate(alglib::real_2d_array data, double train_frac)
{
	std::vector<int> idx( data.rows() );
	for(unsigned int i=0; i<data.rows(); i++) idx[i] = i;
  	
	std::random_shuffle( idx.begin(), idx.end() );

	int n_train = std::round( train_frac*data.rows() );
	this->train.setlength( n_train, data.cols() );
	this->test.setlength( data.rows()-n_train, data.cols() );

	for(unsigned int i=0; i<n_train; i++)
		for(unsigned int j=0; j<data.cols(); j++) this->train[i][j] = data[ idx[i] ][j];

	for(unsigned int i=0; i<data.rows()-n_train; i++)
		for(unsigned int j=0; j<data.cols(); j++) this->test[i][j] = data[ idx[i+n_train] ][j];

	this->n_train = n_train;
	this->n_test = data.rows()-n_train;
}

void class_set::bootstrap(alglib::real_2d_array data, int n)
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
				for(unsigned int j=0; j<idx[i]; j++) 
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

void trainer::bootstrap(alglib::real_2d_array data, int rep)
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
		for(unsigned int i=0; i<rep; i++)
		{
			boot.bootstrap(data, (int)data.rows());
			dfbuildrandomdecisionforest( boot.train, boot.train.rows(), n_var, 2, n_tree, R, info, df, dfrep);
			err = 0.632*dfrelclserror(df, boot.test, boot.test.rows()) + 0.368*dfrep.relclserror;
			#pragma omp critical 
			{
			errors.push_back(err);
			models.push_back(df);
			}
		}
	}
}

void trainer::cv(alglib::real_2d_array data, int folds)
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
		for(unsigned int i=0; i<folds; i++)
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
		}
	}
}

void trainer::rep_cv(alglib::real_2d_array data, int rep, int folds)
{
	if(reset)
	{
		errors.resize(0);
		models.resize(0);
		reset = 0;
	}
	
	for(unsigned int i=0; i<rep; i++) this->cv(data, folds);
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

alglib::decisionforest trainer::best_model()
{
	double mine = HUGE_VAL;
	alglib::decisionforest best = models[0];
	for(unsigned int i=1; i<errors.size(); i++) if(errors[i]<mine) best = models[i];
	
	return best;
}

// SANN implementations

template <class number_type>
int sann::sample_idx(number_type* cweights, int n)
{	
		double rn = alglib::randomreal();
		
		//Special case for zero weights - just return a random index
		if( cweights[n-1] == (number_type)0 ) return floor(rn*n);
		
		rn *= cweights[n-1];
		
		// Binary search to find index
		int low = 0;
		int up = n-1;
		int look = floor(up*0.5);
		
		while(up-low>1)
		{
				if( cweights[look] < rn ) low = look;
				else up = look;
				look = floor( 0.5*(low+up) );
		}
		
		return (cweights[low]<rn) ? up : low;
}

int sann::ELP_idx(double energy)
{	
	if(energy < E_low) return 0;
	else if(energy >= E_up) return n_bins-1;
	else return floor( (energy-E_low)/step );
}

int sann::init()
{
	//Check some allocations and reset optimizer
	
	if( hist == nullptr ) return 0;
	int n_bins = (int) (E_up-E_low)/step + 1;
	for(unsigned int i=0; i<n_bins; i++) hist[i] = 0;
	if( link_idx == nullptr ) return 0;
	for(unsigned int i=0; i<n_link; i++) link_idx[i] = 0;
	action_idx[0] = action_idx[1] = action_idx[2] = 0;
	accepted = iter = 0;
	if(candidates.empty()) return 0;
	if(subs.empty()) return 0;
	
	return 1;
}

sann::sann(alglib::decisionforest df, std::vector<std::string> seqs, int iter_max, int max_link, int n_c, int n_bins, int n_best)
{	
	this->df = df;
	
	// read sequence file
	for(unsigned int i=0; i<seqs.size(); i++) this->seqs.push_back( string_to_array(seqs[i]) );
	n_link = this->seqs.size()+1; 
	
	this->max_link = max_link;
	n_candidates = n_c;
	this->n_best = n_best;
	
	// ELP init
	E_up = 1.0;		// for now we choose the entire range of E
	E_low = 0.0;
	this->n_bins = n_bins;
	step = (E_up - E_low)/n_bins;
	hist = new int[n_bins];
	
	// solutions init
	current.links = new int[n_link*max_link];
	for(unsigned int i=0; i<n_link*max_link; i++) current.links[i] = 20;
	
	//Get initial solution
	current.energy = energy(current.links);
	best_energy = current.energy;
	update_best(current.links, current.energy);
	
	//Boltzmann init
	link_idx = new int[n_link];
	
	this->iter_max = iter_max;
	
	#pragma omp parallel
	{
	srand( 1e6*omp_get_wtime() + alglib::randominteger(RAND_MAX/2)+omp_get_thread_num() );
	}
	
	candidates.resize(n_candidates);
	for(unsigned int i=0; i<n_candidates; i++)
	{
		candidates[i] = new int[n_link*max_link+2];
	}
	energies.resize(n_candidates+1);
	
	// subs init
	read_blosum("blosum80.qij");
	
	if( !init() )
	{ 
		std::cerr<<"Could not initialize!"<<std::endl;
		return;
	}
}

sann::~sann()
{
	delete[] this->hist;
	delete[] this->link_idx;
	if(subs.size()>0) for(unsigned int i=0; i<20; i++) delete[] subs[i];
	for(unsigned int i=0; i<n_candidates; i++) delete[] (candidates[i]);
	
	solution tmp;
	while(!best.empty())
	{
		tmp = best.top();
		delete[] tmp.links;
		
		best.pop();
	}
}

void sann::read_blosum(std::string file)
{
	std::ifstream in(file);
	std::string tmp, AA_order;
	int m_idx, a, b;
	double dval, sum;
	
	// Init subs
	subs.resize(20);
	for(unsigned int i=0; i<20; i++) subs[i] = new double[20];

	// skip comments
	while(in.peek() == '#')
	{
		std::getline(in, tmp);
		continue;
	}

	// read aa ordering used
	for(unsigned int i=0; i<20; i++)
	{
		in>>tmp;
		AA_order += tmp;
	}

	// fill up substitution matrix
	for(unsigned int i=0; i<20; i++)
	{
		for(unsigned int j=0; j<=i; j++)
		{
			in >> dval;
			
			if(i != j)
			{
				a = AAmap.find(AA_order[i])->second;
				b = AAmap.find(AA_order[j])->second;
				subs[a][b] = subs[b][a] = dval;
			}
			else subs[i][j] = 0.0; 
		}
	}
	
	// Normalize and cumulate
	
	for(unsigned int i=0; i<20; i++)
	{
		sum = 0.0;
		for(unsigned int j=0; j<20; j++) 
		{
			sum += subs[i][j];
			subs[i][j] = sum;
		}
		for(unsigned int j=0; j<20; j++) subs[i][j] /= sum;
	}
}

int sann::link_size(linker_set links, int idx)
{
	int stride = idx*max_link;
	
	// Catch some common cases:
	if(links[stride] > 19) return 0;
	if(links[stride+max_link-1] < 20) return max_link;
	
	int i = 0;
	
	while( links[stride+i] < 20 ) i++;
	
	return i;
}

int* sann::assemble_seq(linker_set link)
{   
	int length = 0;
	int write_idx = 1;
	int* l_link = new int[n_link];
	
	for(unsigned int i=0; i<n_link; i++) l_link[i] = link_size(link, i);
	
	for(unsigned int i=0; i<n_link; i++)
	{
		length += l_link[i];
		if(i<n_link-1) length += seqs[i][0];
	}
	
	int* seq = new int[length+1];
	seq[0] = length;
	
	for(unsigned int i=0; i<n_link; i++)
	{
		for(unsigned int j=0; j<l_link[i]; j++) seq[write_idx+j] = link[i*max_link + j];
		write_idx += l_link[i];
		
		if(i<n_link-1)
		{
			for(unsigned int j=1; j<=seqs[i][0]; j++) seq[write_idx+j-1] = seqs[i][j];
			write_idx += seqs[i][0];
		}
	}
	
	delete[] l_link;
	
	return seq;
}

std::string sann::get_seq(int* int_seq)
{
	std::string out;
	
	for(unsigned int i=1; i<=int_seq[0]; i++) out += AA[ int_seq[i] ];
	return out;
}

double sann::energy(linker_set link)
{
	int* seq = assemble_seq(link);
	int* counts;
	double e = 1.0;
	alglib::real_1d_array vars, probs;
	vars.setlength(n_var);
	
	vars[0] = charge(seq);
	vars[1] = Hm(seq);
	vars[2] = pI(seq);
	vars[3] = maxHM(seq);
	vars[4] = charge_var(seq);
	counts = countAA(seq);
	vars[5] = logP(seq, counts);
	for(unsigned int j=n_var-20; j<n_var; j++) vars[j] = counts[j-n_var+20];
	
	for(unsigned int i=0; i<dfs.size(); i++)
	{
		dfprocess(df, vars, probs);
		e *= probs[1];
	}
	
	delete[] seq;
	delete[] counts;
	
	return e;
}

int sann::update_best(linker_set links, double energy)
{
	solution worst;
	if(!best.empty()) worst = best.top();
	
	if( (best.size() < n_best) || (energy < worst.energy) )
	{
		solution new_best;
		
		new_best.links = new int[n_link*max_link];
		for(unsigned int i=0; i<n_link*max_link; i++) new_best.links[i] = links[i];
		new_best.energy = energy;
		
		best.push(new_best);
		if(best.size() > n_best)
		{
			 worst = best.top();
			 best.pop();
			 delete[] worst.links;
		}
		
		if(best_energy > energy) best_energy = energy;
		
		return 1;
	}
	
	return 0;
}

std::string sann::get_best(double error, int full_seq)
{
	std::ostringstream out;
	std::vector<solution> sols;
	std::priority_queue<solution, std::vector<solution>, solution_compare> tmp(best);
	int* int_seq = nullptr;
	
	while(!tmp.empty())
	{ 
		sols.push_back( tmp.top() );
		tmp.pop();
	}
	
	
	if(!full_seq) for(unsigned int i=0; i<n_link; i++) out<<"Link "<<i<<"\t ";
	else out<<"peptide \t\t ";
	out<<"P(CPP)"<<std::endl;
	
	for(int i=sols.size()-1; i>=0; i--)
	{
		if(!full_seq)
		{
			for(unsigned int j=0; j<n_link; j++)
			{
				for(unsigned int k=0; k<max_link; k++)
				{
					if( sols[i].links[j*max_link+k]<20 ) out<<AA[ sols[i].links[j*max_link+k] ];
					else out<<'-';
				}
				out<<"\t ";
			}
		}
		else
		{
			int_seq = assemble_seq(sols[i].links);
			out<<get_seq(int_seq)<<"\t ";
			delete int_seq;
		}
		out<<(1.0-error)*(1.0-sols[i].energy)<<std::endl;
	}
	
	return out.str();
}

std::string sann::sub_to_string()
{
	std::ostringstream out;
	out.precision(3);
	
	if(subs.size()<20)
	{
		out<<"Not loaded!";
		return out.str();
	}
	
	for(unsigned int i=0; i<20; i++) out<<'\t'<<AA[i];
	out<<std::endl;
	
	for(unsigned int i=0; i<20; i++)
	{
		out<<AA[i];
		for(unsigned int j=0; j<20; j++)
		{
			out<<'\t'<<subs[i][j];
		}
		out<<std::endl;
	}
	
	return out.str();
}

std::string sann::get_state()
{
	std::ostringstream out;
	out<<"--> Iteration #"<<iter<<" ("<<iter*n_candidates<<" configurations tested, "<<accepted<<" accepted)"<<std::endl;
	
	double scale = 1.0/n_candidates + (1.0 - 2.0/n_candidates)*iter/iter_max; 
	scale = ( log(scale*n_candidates) - log(1.0-scale) )/E_c;
	out<<"Energy scaling: "<<scale<<", ELP influence: "<<ELP*100.0<<"\%"<<std::endl;
	
	int max_hist=0;
	for(unsigned int i=0; i<(E_up-E_low)/step; i++) if(hist[i]>max_hist) max_hist = hist[i];
	
	out<<"ELP histogram ("<<E_low<<"->"<<E_up<<"):"<<std::endl;
	for(unsigned int i=6; i>0; i--)
	{
		for(unsigned int j=0; j<(E_up-E_low)/step; j++)
		{
			if(hist[j]>i*max_hist/6.0) out<<"*";
			else out<<" ";
		}
		out<<std::endl;
	} 
	for(unsigned int j=0; j<(E_up-E_low)/step; j++) out<<'-';
	out<<std::endl;
	
	out<<"Modifications: ";
	for(unsigned int i=0; i<n_link; i++) out<<"linker "<<i<<": "<<100.0*link_idx[i]/accepted<<"\% ";
	out<<" (";
	
	out<<100.0*action_idx[0]/accepted<<"\% additions, ";
	out<<100.0*action_idx[1]/accepted<<"\% changes, ";
	out<<100.0*action_idx[2]/accepted<<"\% deletions";
	out<<")"<<std::endl;
	
	for(unsigned int i=0; i<n_link; i++)
	{
		for(unsigned int j=0; j<max_link; j++)
		{
			if( current.links[i*max_link+j]<20 ) out<<AA[ current.links[i*max_link+j] ];
			else out<<'-';
		}
		out<<"  ";
	}
	
	out<<"E: "<<current.energy<<" (best: "<<best_energy<<")";
	
	
	return out.str();
}

std::ostream& operator<<(std::ostream& out, sann& opt)
{
	int* int_seq = opt.assemble_seq( opt.current.links );
	
	out<<"Iter: "<<opt.iter<<" ("<<opt.accepted<<" accepted) | Current: "<<opt.get_seq(int_seq)<<
			" (E = "<<opt.current.energy<<") | best P(CPP): "<<1.0-opt.best_energy;
	delete[] int_seq;
	
	return out;
}

int sann::mutate(int aa)
{
	
	double p_blosum = (double)iter/iter_max;	//probability that we will use the blosum sample
	int new_aa;
	
	if( alglib::randomreal()<p_blosum ) new_aa = sample_idx<double>(subs[aa], 20);
	else
	{
		new_aa = alglib::randominteger(19);
		if(new_aa >= aa) new_aa++;
	} 
		
	return new_aa;
}

void sann::cut(linker_set links, int l_idx, int pos)
{
	unsigned int i;
	
	for(i=l_idx*max_link+pos; i<(l_idx+1)*max_link-1; i++)
	{
		links[i] = links[i+1];
		if( links[i] == 20 ) break;
	} 
	links[(l_idx+1)*max_link-1] = 20;
	
	return;
}

void sann::sample_linker(linker_set new_link, linker_set old_link)
{
	for(unsigned int i=0; i<n_link*max_link; i++) new_link[i] = old_link[i];
	
	int cum_action[3];
	
	// First sample which linker to modify
	int l_i = alglib::randominteger(n_link);
	int length_i = link_size(new_link, l_i);
	new_link[n_link*max_link] = l_i;
	
	// Choose the action to perform
	if(new_link[l_i*max_link] > 19)	// empty linker - can only add
	{
		new_link[l_i*max_link] = alglib::randominteger(20);
		new_link[n_link*max_link+1] = 0;
		
		return;
	}
	else if(length_i == max_link) //full linker - can only change or delete
	{
		cum_action[0] = 0;
		cum_action[1] = action_idx[1] + entropy_add;
		cum_action[2] = action_idx[1] + action_idx[2] + entropy_add;
		
		int a_i = sample_idx<int>(cum_action, 3);
		int pos = alglib::randominteger(length_i);
		
		switch(a_i)
		{
			case 1: new_link[l_i*max_link+pos] = mutate( new_link[l_i*max_link+pos] );
					new_link[n_link*max_link+1] = 1;
					break;
			case 2: cut(new_link, l_i, pos); 
					new_link[n_link*max_link+1] = 2;
					break;
		}
		
		return;
	}
	else //neither full nor empty linker - can add, change or delete
	{
		cum_action[0] = action_idx[0] + entropy_add;
		cum_action[1] = action_idx[0] + action_idx[1] + entropy_add;
		cum_action[2] = action_idx[0] + action_idx[1] + action_idx[2] + entropy_add;
		
		int a_i = sample_idx<int>(cum_action, 3);
		int pos = alglib::randominteger(length_i);
		
		switch(a_i)
		{
			case 0:	new_link[l_i*max_link + length_i] = alglib::randominteger(20);
					new_link[n_link*max_link+1] = 0;
					break;
			case 1: new_link[l_i*max_link+pos] = mutate( new_link[l_i*max_link+pos] );
					new_link[n_link*max_link+1] = 1;
					break;
			case 2: cut(new_link, l_i, pos);
					new_link[n_link*max_link+1] = 2;
					break;
		}
		
		return;
	}	
}

void sann::anneal()
{
	int h_idx; 
	int min_idx=-1;
	double min_E = current.energy;
	
	energies[n_candidates] = min_E;
	double *cum_e = new double[n_candidates+1];
	double scale = 1.0/n_candidates + (1.0 - 2.0/n_candidates)*iter/iter_max; 
	scale = ( log(scale*n_candidates) - log(1.0-scale) )/E_c; 
	double w;
	
	#pragma omp parallel for private(w, h_idx)
	for(unsigned int i=0; i<n_candidates; i++)
	{
		sample_linker(candidates[i], current.links);
		energies[i] = energy(candidates[i]);
		
		
		#pragma omp critical
		{
		if(energies[i] < min_E)
		{
			min_E = energies[i];
			min_idx = i;
		}
		}
		
		// Calculate weights
		h_idx = ELP_idx(energies[i]);
		
		w = energies[i] + (1.0 - (double)iter/iter_max)*ELP*hist[h_idx]/iter;
		cum_e[i] = exp(-w*scale); 
	}
	
	// Add current solution to the list of candidates
	h_idx = ELP_idx( energies[n_candidates] );
	w = energies[n_candidates] + (1.0 - (double)iter/iter_max)*ELP*hist[h_idx]/iter;
	cum_e[n_candidates] = exp(-w*scale);
	
	if(min_idx >= 0) update_best(candidates[min_idx], energies[min_idx]);
	
	for(unsigned int i=0; i<n_candidates+1; i++)
		if(i>0) cum_e[i] += cum_e[i-1];
	
	int new_sol = sample_idx<double>(cum_e, n_candidates+1);
	delete[] cum_e;
	
	if(new_sol < n_candidates)
	{
		accepted++;
		link_idx[ candidates[new_sol][n_link*max_link] ]++;
		action_idx[ candidates[new_sol][n_link*max_link+1] ]++;
		for(unsigned int i=0; i<n_link*max_link; i++) current.links[i] = candidates[new_sol][i]; 
		current.energy = energies[new_sol];
	}
	
	h_idx = ELP_idx( current.energy );
	hist[h_idx]++;
	
	iter++;
	
}
