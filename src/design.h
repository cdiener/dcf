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
#include <queue>

//Whether we can use curses for nicer output
#ifdef CURSES_HAVE_CURSES_H
#include <curses.h>
#endif

/**
General prototypes and definitions.
*/

const double R = 0.6; // Fraction of data used to construct individual trees

const int n_var = 26; //Number of variables used for classification
const int entropy_add = 8;
const double E_c = 0.005; // The largest competing energy difference that we still want to accept with P = 1 - min_acc
const double ELP = 1.0; // Factor how strong the EP should influence the energy values (0-none, 1-as much as the energy values)  

/** 
Helper functions and types
*/

/**
 * Represents a linker set
 */
typedef int* linker_set;

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

/* Classification */

/**
Function that reads positive and negative examples for a classification into the 
correct data structure.

@param pos Path to positive file.
@param neg Path to negative file.
@return The mean hydrophobicity of the sequence
*/
alglib::real_2d_array read_vars(std::string pos, std::string neg);
 
/**
Represents a data set for classification consisting of a training and test set.
*/
class class_set 
{
	public:
		alglib::real_2d_array train, test;

		/**
		Function that generates training and test set for classification

		@param data Full classification data.
		@param train_frac Fraction of full data set used for training (between 0 and 1).
		*/
		void generate(alglib::real_2d_array data, double train_frac);
		
		/**
		 * Function that generates a bootstrap sample from the data
		 * 
		 * @param data Full classification data.
		 * @param n Size of bootstrap sample.
		 */
		 void bootstrap(alglib::real_2d_array data, int n);
		 
		
		/**
		Constructors. Only initializes random number seed.
		*/	
		class_set();
		
		int n_train, n_test;
};

/**
 * Represents a classification trainer that also approximates the error using bootstrap or cross validation.
 */
class trainer
{
	public:
		std::vector<alglib::decisionforest> models;
		std::vector<double> errors;
		int n_tree, reset;
	
		trainer(int n_tree):n_tree(n_tree){ reset=1; };
	
		/**
		 * Estimate using bootstrap .632
		 * 
		 * @param data The full data set to be used.
		 * @param rep How often to repeat the bootstrap.
		 */
		void bootstrap(alglib::real_2d_array data, int rep);
		
		/**
		 * Estimate using cross validation
		 * 
		 * @param data The full data set to be used.
		 * @param folds How many folds to use.
		 */
		void cv(alglib::real_2d_array data, int folds);
		
		/**
		 * Estimate using repeated cross validation
		 * 
		 * @param data The full data set to be used.
		 * @param rep How many times to repeat cross validation.
		 * @param folds How many folds to use.
		 */
		void rep_cv(alglib::real_2d_array data, int rep, int folds);
		
		/**
		 * Return error mean.
		 *
		 * @return mean error of classification.
		 */
		double mean_error();
		
		/**
		 * Return error standard deviation.
		 *
		 * @return sd of classification error.
		 */
		double sd_error();
		
		/**
		 * Return minimum error.
		 *
		 * @return minimum error.
		 */
		double min_error();  
		
		/**
		 * Return the best model found.
		 * 
		 * @return Best found model.
		 */
		alglib::decisionforest best_model();  
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
		std::vector<int*> seqs;
		int n_link;
		int max_link;
		int n_candidates;
		int n_best;
		int accepted, iter, iter_max;
		std::vector<alglib::decisionforest> dfs;
		std::vector<linker_set> candidates;
		std::vector<double> energies;
		
		// Current solution
		solution current;
		
		// Best solutions queue
		std::priority_queue<solution, std::vector<solution>, solution_compare> best;
		double best_energy;
		
		// Boltzmann sample histograms
		int* link_idx;		//index distributions - this is not used for sampling 
		int action_idx[3];	//action distributions - actions are [ 0=add, 1=change, 2=delete ]
		
		// BLOSUM-based substitution
		std::vector<double*> subs;
		
		// ELP histogram
		double E_low; //lower energy value
		double E_up;  //upper energy value
		double step;  //step size for bins
		int n_bins;	  //number of bins
		int* hist;	  //actual histogram
		
		// Helper functions
		/**
		 * Function to sample an index from a weighted list of indices.
		 * 
		 * @param cweights Array of ints which denote the cumulative weights of the entries.
		 * @param n Number of indices.
		 * @return An index ranging from 
		 */
		template <class number_type> int sample_idx(number_type* cweights, int n);
		
		/**
		 * Function to get the index of an energy value in ELP histogram
		 * 
		 * @param energy The energy value
		 * @return The index
		 */
		int ELP_idx(double energy);
		
	public:
		/**
		 * Constructor
		 * 
		 * @param seq_file List of sequences to be linked into a CPP
		 * @param max_link Maximum length of the individual linkers
		 * @param n_c Number of candidates to be generated in each annealing step
		 * @param n_bins Number of bins for the ELP histograms
		 */
		sann(std::vector<alglib::decisionforest> dfs, std::vector<std::string> seqs, int iter_max=500, int max_link=8, 
				int n_c=100, int n_bins=64, int n_best=8);
		
		/**
		 * Destructor.
		 */
		~sann();
	
		/**
		 * Function to initialize and broadcast basic variables and classifiers
		 * 
		 * @return 0 - failure, 1 - success
		 */
		int init();
	
		/**
		 * Function to calculate the energy of a sequence configuration. This is basically the inverse probability to be
		 * a CPP and have high efficiency.
		 * 
		 * @param link The linker configuration.
		 * @return The "energy", meaning E = 1-Pr(CPP,efficient)
		 */
		double energy(linker_set link);
		
		/** 
		 * Update list of best solutions
		 * 
		 * @param links The linker set.
		 * @param energy The energy.
		 * @return 1 if there was an update, 0 if not.
		 */
		int update_best(linker_set links, double energy);
		
		/**
		 * Gets the best found solutions.
		 * 
		 * @param error Classification error to incorporate.
		 * @param full_seq Whether the full sequence or only linkers are printed.
		 * @return Formatted text containing the best solutions.
		 */
		std::string get_best(double error=0.0, int full_seq=1);  
		
		/**
		 * Get the length of a linker in a linker set 
		 * 
		 * @param links The linker set.
		 * @param idx index of linker.
		 * @return length of linker.
		 */
		int link_size(linker_set links, int idx);
		
		/**
		 * Function to sample a new linker configuration from a previous one.
		 * 
		 * @param link The linker configuration.
		 * @return A new linker configuration.
		 */
		void sample_linker(linker_set new_link, linker_set old_link);
		 
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
		 void cut(linker_set links, int l_idx, int pos);
		 
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
		 int* assemble_seq(linker_set link);
		 
		/**
		 * Converts integer sequence to amino acid sequence.
		 * 
		 * @param int_seq The integer sequence.
		 * @return The amino acid sequence as string.
		 */
		std::string get_seq(int* int_seq);
		
		/**
		 * Gets some basic diagnostics of the optimization.
		 * 
		 * @return Diagnostics in string format.
		 */
		std::string get_state();
		
		#ifdef CURSES_HAVE_CURSES_H
		/**
		 * Gets detailed diagnostics in a nice curses output
		 */
		void get_curses();
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
 
