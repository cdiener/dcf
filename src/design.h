/*
 * design.h
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
 
#ifndef __DESIGN_H__
#define __DESIGN_H__

#include "alglib/ap.h"
#include "alglib/dataanalysis.h"
#include "sys_config.h"
#include <string>
#include <vector>
#include <queue>
//#include <memory>

//Whether we can use curses for nicer output
#ifdef CURSES_HAVE_CURSES_H
#include <curses.h>
#endif

/**
General prototypes and definitions.
*/

const int entropy_add = 8;
const double E_c = 0.05; // The largest competing energy difference that we still want to accept with P = 1 - min_acc
const double ELP = 1.5; // Factor how strong the EP should influence the energy values (0-none, 1-as much as the energy values)  

/** 
Helper functions and types
*/

/**
 * Represents a linker set
 */
typedef std::vector<std::vector<int> > linker_set;

/**
 * Represents a candidate for the optimization
 */
typedef struct{ linker_set links; int li; int ai; } cand;

/**
 * A solution for the optimization
 */ 
typedef struct{ linker_set links; double energy; } solution;

/**
 * compare class for solution constructs
 */
class solution_compare
{
	public:
		bool operator() (const solution& a, const solution& b) { return (a.energy<b.energy); }
};

/* Optimization */

/**
 * Class for the simulated annealing optimizer looking for the optimal linker
 * configuration
 */
class sann
{
	private:
		// General variables
		std::vector<std::vector<int> > seqs;
		int n_link;
		std::vector<int> vmax_link;		// maximum lengths of the linkers
		std::vector<int> mod_link;		// indices of non-zero linkers
		int n_candidates;
		int n_best;
		int accepted, iter, iter_max;
		std::vector<alglib::decisionforest> dfs;
		std::vector<cand> candidates;
		std::vector<double> energies;
		
		// Current solution
		solution current;
		
		// Best solutions queue
		std::vector<solution> best;
		double best_energy;
		
		// Boltzmann sample histograms
		std::vector<int> link_idx;		//index distributions - this is currently *not* used for sampling 
		std::vector<int> action_idx;	//action distributions - actions are [ 0=add, 1=change, 2=delete ]
		
		// BLOSUM-based substitution
		std::vector<std::vector<double> > subs;
		
		// ELP histogram
		double E_low; //lower energy value
		double E_up;  //upper energy value
		double step;  //step size for bins
		int n_bins;	  //number of bins
		std::vector<int> hist;	  //actual histogram
		
		// Helper functions
		/**
		 * Function to sample an index from a weighted list of indices.
		 * 
		 * @param cweights Array of ints which denote the cumulative weights of the entries.
		 * @param n Number of indices.
		 * @return An index ranging from 
		 */
		template <class number_type> int sample_idx(const std::vector<number_type>& cweights, int n);
		
		/**
		 * Function to get the index of an energy value in ELP histogram
		 * 
		 * @param energy The energy value
		 * @return The index
		 */
		int ELP_idx(double energy);
		
		// Curses variables
		#ifdef CURSES_HAVE_CURSES_H
		WINDOW* curses_wins[3];
		#endif
		
	public:
		/**
		 * Constructor
		 * 
		 * @param seq_file List of sequences to be linked into a CPP
		 * @param max_link Maximum length of the individual linkers
		 * @param n_c Number of candidates to be generated in each annealing step
		 * @param n_bins Number of bins for the ELP histograms
		 */
		sann(const std::vector<alglib::decisionforest>& dfs, const std::vector<std::string>& seqs, 
				const std::vector<int>& max_link, int iter_max=500, int n_c=100, int n_bins=64, 
					int n_best=8);
		
		/**
		 * Destructor.
		 */
		~sann(){};
	
		/**
		 * Function to calculate the energy of a sequence configuration. This is basically the inverse probability to be
		 * a CPP and have high efficiency.
		 * 
		 * @param link The linker configuration.
		 * @return The "energy", meaning E = 1-Pr(CPP,efficient)
		 */
		double energy(const linker_set& link);
		
		/** 
		 * Update list of best solutions
		 * 
		 * @param links The linker set.
		 * @param energy The energy.
		 * @return 1 if there was an update, 0 if not.
		 */
		int update_best(const linker_set& links, double energy);
		
		/**
		 * Gets the best found solutions.
		 * 
		 * @param error Classification error to incorporate.
		 * @param full_seq Whether the full sequence or only linkers are printed.
		 * @return Formatted text containing the best solutions.
		 */
		std::string get_best(double error=0.0, int full_seq=1);  
		
		/**
		 * Function to sample a new linker configuration from a previous one.
		 * 
		 * @param link The linker configuration.
		 * @return A new linker configuration.
		 */
		void sample_linker(cand& c, const linker_set& old_link);
		 
		 /**
		  * Mutate an aminoacid for another one.
		  * 
		  * @param aa Integer code of the current amino acid.
		  * æreturn Integer code of new aminoacid.
		  */
		 int mutate(int aa);
		 
		 /**
		  * Remove a random amino acid from a linker
		  * 
		  * @param links The linker set to be modified.
		  * @param l_idx The linker index to be modified.
		  * @param pos Position of the removed amino acid
		  */
		 void cut(linker_set& links, int l_idx, int pos);
		 
		 /**
		  * A single annealing run which consists of generating several solution candidates in parallel, than sampling a new
		  * solution from all the generated + the current solution and updating the top solution list.
		  */
		 void anneal();
		 
		 /**
		  * Read BLOSUM qij odds-ratios and build up a matrix of substitution probabilities.
		  * 
		  * @param file The filename of the *.qij file as returned by the blosum program.
		  */
		 void read_blosum(std::string file);
		 
		 /**
		  * Formats the substitution matrix as a string for nice outputting.
		  * 
		  * @return Formatted matrix as string.
		  */
		 std::string sub_to_string();
		 
		 /**
		  * Assembles the complete peptide sequence in int coding.
		  * 
		  * @param link The linker set to be used.
		  * @return The assembled sequence as integer array.
		  */
		 std::vector<int> assemble_seq(const linker_set& link);
		 
		/**
		 * Converts integer sequence to amino acid sequence.
		 * 
		 * @param int_seq The integer sequence.
		 * @return The amino acid sequence as string.
		 */
		std::string get_seq(const std::vector<int>& int_seq);
		
		/**
		 * Gets some basic diagnostics of the optimization.
		 * 
		 * @return Diagnostics in string format.
		 */
		std::string get_state();
		
		#ifdef CURSES_HAVE_CURSES_H
		/**
		 * Initializes the curses screen
		 */ 
		int init_curses();
		
		/**
		 * Gets detailed diagnostics in a nice curses output
		 */
		void update_curses();
		
		/**
		 * Destroys the curses screen
		 */
		 void end_curses();
		#endif  
		
		/**
		 * Gets current energy value
		 */
		double current_energy(){ return current.energy; };
		
		/**
		 * Gets best found energy value
		 */
		 double get_best_energy(){ return best_energy; };
		 
		 /**
		  * Gets number of accepted solutions
		  */
		 int get_accepted(){ return accepted; };
		
		//Operators
		
		/**
		 * Output operator.
		 */
		friend std::ostream& operator<<(std::ostream& out, sann& opt);
};

#endif /* __DESIGN_H__ */
 
