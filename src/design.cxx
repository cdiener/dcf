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
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

// SANN implementations

template <class number_type>
int sann::sample_idx(const std::vector<number_type>& cweights, int n)
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

sann::sann(const std::vector<alglib::decisionforest>& dfs, const std::vector<std::string>& seqs, const std::vector<int>& max_link, int iter_max, int n_c, int n_bins, int n_best)
{	
	this->dfs = dfs;
	
	// read sequence file
	for(unsigned int i=0; i<seqs.size(); i++) this->seqs.push_back( string_to_array(seqs[i]) );
	n_link = this->seqs.size()+1; 
	
	this->vmax_link = std::vector<int>(max_link);
	for(unsigned int i=0; i<n_link; i++) if(vmax_link[i] > 0) mod_link.push_back(i);
	if( mod_link.empty() )
	{
		std::cerr<<"All linkers have zero length. Nothing to optimize!"<<std::endl;
		throw "error";
	}
	
	n_candidates = n_c;
	this->n_best = n_best;
	
	// ELP init
	E_up = 1.0;		// for now we choose the entire range of E
	E_low = 0.0;
	this->n_bins = n_bins;
	step = (E_up - E_low)/n_bins;
	hist = std::vector<int>(n_bins, 0);
	
	//Get initial solution
	current.links = std::vector<std::vector<int> >(n_link, std::vector<int>());
	current.energy = energy(current.links);
	best_energy = current.energy;
	update_best(current.links, current.energy);
	
	//Boltzmann init
	link_idx = std::vector<int>(n_link, 0);
	
	//Action init
	action_idx = std::vector<int>(3,0);
	
	this->iter_max = iter_max;
	
	// Innitialize random number seeds (different for each thread)
	#pragma omp parallel
	{
		srand( 1e6*omp_get_wtime() + alglib::randominteger(RAND_MAX/2)+omp_get_thread_num() );
	}
		
	// Initialize candidate set
	candidates.resize(n_c);
	for(unsigned int i=0; i<n_c; i++) 
		candidates[i].links = linker_set( n_link, std::vector<int>() );
	
	energies.resize(n_candidates+1);
	
	// subs init
	read_blosum("blosum80.qij");
	
	accepted = iter = 0;
	for(unsigned int i=0; i<n_link; i++) std::cout<<"vmax "<<i<<": "<<vmax_link[i]<<std::endl;
}

void sann::read_blosum(std::string file)
{
	std::ifstream in(file);
	std::string tmp, AA_order;
	int m_idx, a, b;
	double dval, sum;
	
	if(!in.good())
	{
		std::cerr<<"Could not open blosum file!!!"<<std::endl;
		throw "error";
	}
	
	// Init subs
	subs.resize(20);
	for(unsigned int i=0; i<20; i++) subs[i] = std::vector<double>(20);

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

std::vector<int> sann::assemble_seq(const linker_set& link)
{   
	int length = 0;
	int write_idx = 0;
	std::vector<int> l_link = std::vector<int>(n_link);
	
	for(unsigned int i=0; i<n_link; i++) l_link[i] = link[i].size();
	
	for(unsigned int i=0; i<n_link; i++)
	{
		length += l_link[i];
		if(i<n_link-1) length += seqs[i].size();
	}
	
	std::vector<int> seq = std::vector<int>(length);
	
	for(unsigned int i=0; i<n_link; i++)
	{
		for(unsigned int j=0; j<l_link[i]; j++) seq[write_idx+j] = link[i][j];
		write_idx += l_link[i];
		
		if(i<n_link-1)
		{
			for(unsigned int j=0; j<seqs[i].size(); j++) seq[write_idx+j] = seqs[i][j];
			write_idx += seqs[i].size();
		}
	}
	
	return seq;
}

std::string sann::get_seq(const std::vector<int>& int_seq)
{
	std::string out;
	
	for(unsigned int i=0; i<int_seq.size(); i++) out += AA[ int_seq[i] ];
	return out;
}

double sann::energy(const linker_set& link)
{
	std::vector<int> seq = assemble_seq(link);
	std::vector<int> counts;
	double p = 1.0;
	alglib::real_1d_array vars;
	alglib::real_1d_array probs;
	vars.setlength(n_var);
	
	vars[0] = charge(seq);
	vars[1] = Hm(seq);
	vars[2] = pI(seq);
	vars[3] = maxHM(seq);
	vars[4] = charge_var(seq);
	counts = countAA(seq);
	vars[5] = logP(seq, counts);
	for(unsigned int j=n_var-20; j<n_var; j++) vars[j] = 1.0*counts[j-n_var+20]/seq.size();
	
	for(unsigned int i=0; i<dfs.size(); i++)
	{
		dfprocess(dfs[i], vars, probs);
		p *= probs[0];
	}
	
	return 1.0-p;
}

int sann::update_best(const linker_set& links, double energy)
{
	solution worst;
	if(!best.empty()) worst = best[0];
	
	if( (best.size() < n_best) || (energy < worst.energy) )
	{
		if(best_energy > energy) best_energy = energy;
		
		solution new_best;
		new_best.links = links;
		new_best.energy = energy;
		
		// Check for duplicates
		for(unsigned int i=0; i<best.size(); i++)
			if(best[i].links == new_best.links) return 0;
		
		if(best.size() < n_best)
		{
			best.push_back(std::move(new_best));
			std::push_heap(best.begin(), best.end(), solution_compare());
		}
		else 
		{
			std::pop_heap(best.begin(), best.end(), solution_compare());
			best[best.size()-1] = std::move(new_best);
			std::push_heap(best.begin(), best.end(), solution_compare());
		}
		
		return 1;
	}
	
	return 0;
}

std::string sann::get_best(double error, int full_seq)
{
	std::ostringstream out;
	std::vector<solution> sols;
	std::vector<solution> tmp(best);
	std::vector<int> int_seq;
	
	while(!tmp.empty())
	{ 
		std::pop_heap(tmp.begin(), tmp.end(), solution_compare());
		sols.push_back( tmp.back() );
		tmp.pop_back();
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
				for(unsigned int k=0; k<sols[i].links[j].size(); k++) out<<AA[ sols[i].links[j][k] ];
				out<<"\t ";
			}
		}
		else
		{
			int_seq = assemble_seq(sols[i].links);
			out<<get_seq(int_seq)<<"\t ";
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

#ifdef CURSES_HAVE_CURSES_H
void sann::get_curses()
{
	int rows, cols, lower, upper, left, right, nbin;
	nbin = (E_up-E_low)/step;
	getmaxyx(stdscr, rows, cols);
	lower = std::min(rows-10, 16);
	upper = 5;
	left = 5;
	right = std::min(cols-left, nbin+left);
	clear();
	refresh();
	
	int max_hist=0;
	for(unsigned int i=0; i<nbin; i++) if(hist[i]>max_hist) max_hist = hist[i];
	
	mvaddch(0, left, ACS_DIAMOND);
	mvprintw(0, left+1, "Iteration #%d (%d configurations tested, %d accepted)", iter, iter*n_candidates, accepted);
	double scale = 1.0/n_candidates + (1.0 - 2.0/n_candidates)*iter/iter_max; 
	scale = ( log(scale*n_candidates) - log(1.0-scale) )/E_c;
	mvprintw(1, left+1, "Energy scaling: %g, ELP influence: %d%%", scale, (int)round(ELP*100.0));
	mvprintw(3, left, "ELP histogram");
	
	mvaddch(upper, left, ACS_UARROW);
	mvaddstr(upper, left-2, "E");
	for(unsigned int i=upper+1; i<lower; i++) mvaddch(i, left, ACS_VLINE);
	mvaddch(lower,left, ACS_LLCORNER);
	mvaddstr(lower+1, left, "0");
	mvaddstr(lower+1, right, "1");
	for(unsigned int i=left+1; i<right; i++) mvaddch(lower, i, ACS_HLINE);
	mvaddch(lower, right, ACS_RARROW);
	
	double w = right-left;
	double h = lower-upper;
	
	for(unsigned int i=h; i>0; i--)
	{
		for(unsigned int j=0; j<nbin; j++)
		{
			if( hist[j]>=i*max_hist/(double)h ) mvaddch(lower-i, left+1+j/(double)nbin*w, (char)254);
		}
	}	
	
	std::ostringstream out;
	out<<"Modifications: ";
	for(unsigned int i=0; i<n_link; i++) out<<"linker "<<i<<": "<<100.0*link_idx[i]/accepted<<"%% ";
	out<<std::endl<<"\t\t(";
	out<<100.0*action_idx[0]/accepted<<"%% additions, ";
	out<<100.0*action_idx[1]/accepted<<"%% changes, ";
	out<<100.0*action_idx[2]/accepted<<"%% deletions";
	out<<")";
	std::string tmp_str = out.str(); 
	mvprintw(lower+3, left, tmp_str.c_str() );
	mvprintw(lower+6, left, "E: %g (best: %g)", current.energy, best_energy);
	
	refresh();
}
#endif

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
		for(unsigned int j=0; j<vmax_link[i]; j++) out<<AA[ current.links[i][j] ];
		out<<"  ";
	}
	
	out<<"E: "<<current.energy<<" (best: "<<best_energy<<")";
	
	
	return out.str();
}

std::ostream& operator<<(std::ostream& out, sann& opt)
{
	std::vector<int> int_seq = opt.assemble_seq( opt.current.links );
	
	out<<"Iter: "<<opt.iter<<" ("<<opt.accepted<<" accepted) | Current: "<<opt.get_seq(int_seq)<<
			" (E = "<<opt.current.energy<<") | best P(CPP): "<<1.0-opt.best_energy;
	
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

void sann::cut(linker_set& links, int l_idx, int pos)
{
	links[l_idx].erase( links[l_idx].begin()+pos );
	
	return;
}

void sann::sample_linker(cand& c, const linker_set& old_link)
{
	c.links = old_link;
	
	std::vector<int> cum_action(3);
	
	// First sample which linker to modify
	int l_i = 0;
	if(mod_link.size() == 1) l_i = mod_link.back(); 
	else l_i = mod_link[ alglib::randominteger( mod_link.size() ) ];
	
	int length_i = c.links[l_i].size();
	if(length_i > vmax_link[l_i]) std::cout<<"Something went wrong!"<<std::endl;
	c.li = l_i;
	
	// Choose the action to perform
	if( length_i == 0 )	// empty non-zero linker - can only add
	{
		c.links[l_i].emplace_back( alglib::randominteger(20) );
		c.ai = 0;
		
		return;
	}
	else if(length_i == vmax_link[l_i]) //full linker - can only change or delete
	{
		cum_action[0] = 0;
		cum_action[1] = action_idx[1] + entropy_add;
		cum_action[2] = action_idx[1] + action_idx[2] + entropy_add;
		
		int a_i = sample_idx<int>(cum_action, 3);
		int pos = alglib::randominteger(length_i);
		
		switch(a_i)
		{
			case 1: c.links[l_i][pos] = mutate( c.links[l_i][pos] );
					c.ai = 1;
					break;
			case 2: cut(c.links, l_i, pos); 
					c.ai = 2;
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
			case 0:	c.links[l_i].emplace_back( alglib::randominteger(20) );
					c.ai = 0;
					break;
			case 1: c.links[l_i][pos] = mutate( c.links[l_i][pos] );
					c.ai = 1;
					break;
			case 2: cut(c.links, l_i, pos);
					c.ai = 2;
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
	std::vector<double> cum_e(n_candidates+1);
	double scale = 1.0/n_candidates + (1.0 - 2.0/n_candidates)*iter/iter_max; 
	scale = ( log(scale*n_candidates) - log(1.0-scale) )/E_c; 
	double w;
	
	#pragma omp parallel for private(w, h_idx)
	for(unsigned int i=0; i<n_candidates; i++)
	{
		sample_linker(candidates[i], current.links);
		energies[i] = energy(candidates[i].links);
		
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
	
	if(min_idx >= 0) update_best(candidates[min_idx].links, energies[min_idx]);
	
	for(unsigned int i=0; i<n_candidates+1; i++)
		if(i>0) cum_e[i] += cum_e[i-1];
	
	int new_sol = sample_idx<double>(cum_e, n_candidates+1);
	
	if(new_sol < n_candidates)
	{
		accepted++;
		link_idx[ candidates[new_sol].li ]++;
		action_idx[ candidates[new_sol].ai ]++;
		current.links = candidates[new_sol].links; 
		current.energy = energies[new_sol];
	}
	
	h_idx = ELP_idx( current.energy );
	hist[h_idx]++;
	
	iter++;
	
}
