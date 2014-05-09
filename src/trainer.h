/*
 * trainer.h
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

#ifndef __TRAINER_H__
#define __TRAINER_H__

#include "alglib/ap.h"
#include "alglib/dataanalysis.h"
#include "sys_config.h"
#include <string>

const double R = 0.66;	// Fraction of data used to construct individual trees

/* Classification */

/**
Function that reads positive and negative examples for a classification into the 
correct data structure.

@param pos Path to positive file.
@param neg Path to negative file.
@return The mean hydrophobicity of the sequence
*/
alglib::real_2d_array read_vars(const std::string pos, const std::string neg);
 
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
		void generate(const alglib::real_2d_array data, double train_frac);
		
		/**
		 * Function that generates a bootstrap sample from the data
		 * 
		 * @param data Full classification data.
		 * @param n Size of bootstrap sample.
		 */
		 void generate_bootstrap(const alglib::real_2d_array data, int n);
		 
		
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
	
		trainer(int n_tree);
	
		/**
		 * Estimate using bootstrap .632
		 * 
		 * @param data The full data set to be used.
		 * @param rep How often to repeat the bootstrap.
		 */
		void bootstrap(const alglib::real_2d_array data, int rep);
		
		/**
		 * Estimate using cross validation
		 * 
		 * @param data The full data set to be used.
		 * @param folds How many folds to use.
		 */
		void cv(const alglib::real_2d_array data, int folds);
		
		/**
		 * Estimate using repeated cross validation
		 * 
		 * @param data The full data set to be used.
		 * @param rep How many times to repeat cross validation.
		 * @param folds How many folds to use.
		 */
		void rep_cv(const alglib::real_2d_array data, int rep, int folds);
		
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
		alglib::decisionforest best_model() const;  
};
 

#endif
