/*
 * iztli.h
 * This file is part of Iztli
 *
 * Copyright (C) 2013 - Christian Diener
 *
 * Iztli is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Iztli is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with optim. If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef __IZTLI_H__
#define __IZTLI_H__

#include <string>
#include <map>
#include <vector>
#include <cmath>

/**
General prototypes and definitions.

Note that we encode the amoinoacids internally as integer vectors as this allows the sequence to be used as
a hash and makes most access operations O(1).
*/

//The amino acid coding
const std::string AA = "ILWFVMYAPTSCGNDQEHKR";
const std::map<char, int> AAmap = { {'I',0},{'L',1},{'W',2},{'F',3},{'V',4},{'M',5},{'Y',6},{'A',7},
		{'P',8},{'T',9},{'S',10}, {'C',11},{'G',12},{'N',13},{'D',14},{'Q',15},{'E',16},{'H',17},
		{'K',18},{'R',19} };

/**
Some hydrophobicity scales. Eisenberg is taken from Eisenberg. et. al. and
CCS is taken from Tossi. et. al.
*/
const double eisenberg[] = { 0.73, 0.53, 0.37, 0.61, 0.54, 0.26, 0.02, 0.25, -0.07, -0.18, -0.26,
								0.04, 0.16, -0.64, -0.72, -0.69, -0.62, -0.40, -1.10, -1.76 };

const double CCS[] = { 8.7, 9.7, 9.7, 10.0, 4.1, 4.6, 2.5, -1.1, -0.2, -3.8, -4.3, -2.3, -2.4, 
						-7.1, -8.3, -6.0, -8.3, -3.8, -9.9, -10.0};

#define SCALE CCS //define the used scale

/**
Linear model coefficients to calculate the partition coefficient logP.
Taken from Tao et. al. 1999
*/
const double logP_coef[] = {-3.25, 0.7, 0.8, 1.47, 1.16, 0.32, 0.51, 0.55, -0.27, 0.15, -0.26, -0.45, 
								0.83, -0.22, -0.98, -0.28, -1.0, -0.34, -0.31, 0.17, -0.79};

/**
 * Propensities of the amino acids to participate in an alpha-helix. Taken from
 * averaging the propensities of Jiang (1996) over alpha, alpha+beta and alpha/beta
 * proteins.
 */
const double alpha_prop[] = {1.02, 1.34, 0.95, 0.93, 1.01, 1.33, 0.9, 1.42, 0.4, 0.79, 0.75, 0.78, 0.44, 0.8, 0.91, 1.31, 1.38, 0.9, 1.04, 1.29}; 

const double alpha_cut = 1.0;
 
/**
 Parameters for the pI calculation. For speed we precompute the powers of 10 which are used repeatedly which leaves us
with a single power calculation per pH value (10^pH). pKa values are taken from 
http://en.wikipedia.org/wiki/List_of_standard_amino_acids .
*/
const double pKC = pow(10.0, 8.18);
const double pKD = pow(10.0, 3.9);
const double pKE = pow(10.0, 4.07);
const double pKY = pow(10.0, 10.46);
const double pKCOOH = pow(10.0, 3.65);
const double pKH = pow(10.0, 6.04);
const double pKK = pow(10.0, 10.54);
const double pKR = pow(10.0, 12.48);
const double pKNH2 = pow(10.0, 8.2);

//Constants
const int W = 7; //window size
const double TOLPI = 0.1; //Bisection tolerance for pI
const double ANGINC = 1.0/180.0*M_PI; //angle increment
const double ANGMIN = 95.0/180.0*M_PI; //minimum angle
const double ANGMAX = 108.0/180.0*M_PI; //maximum angle

//Function prototypes

/**
Calculates the mean hydrophobicity of a peptide over windows of size W.

@param seq A peptide sequence in integer encoding
@return The mean hydrophobicity of the sequence
*/
double Hm(const std::vector<int>& seq);

/**
Calculates the charge of a peptide

@param seq A peptide sequence in integer encoding
@return The average charge of the peptide. 
*/
double charge(const std::vector<int>& seq);

/**
Calculates the charge variation of a peptide

@param seq A peptide sequence in integer encoding
@return The charge variation of the peptide. 
*/
double charge_var(const std::vector<int>& seq);


/**
Calculates the isoelectric point of the peptide via the Henderson-Hasselbalch
equation. The strategy is to begin with some rough bisection steps and than
identify the zero point by linear interpolation.

@param seq A peptide sequence in integer encoding
@return The isoelectric point of the peptide
*/
double pI(const std::vector<int>& seq);

/**
Calculates the range of the maximum hydrophobic moment over windows of size W.

@param seq A peptide sequence in integer encoding
@return The maximum hydrophobic moment
*/
double rangeHM(const std::vector<int>& seq);

/**
Calculates the number of individual aminoacids

@param seq A peptide sequence in integer encoding
@return Array of amino acid counts
*/
std::vector<int> countAA(const std::vector<int>& seq);

/**
Approximates the water-octanol partition coefficient logP. 

@param seq A peptide sequence in integer encoding
@param counts The amino acid counts
@return logP value
*/
double logP(const std::vector<int>& counts);

/**
 * Approximates the proportion of alpha-helices in the peptide
 * using the Chou-Fasmann method with updated propensities.
 * 
 * @param seq The sequence in integer coding
 * @return proportion of alpha-helices.
 */ 
double alpha(const std::vector<int>& seq);

/**
Helper function to convert a protein sequence to an integer array

@param seq A sequence in string format
@return The sequence in integer format
*/
std::vector<int> string_to_array(const std::string& seq);

/**
Helper function to read a list of proteins from a text file.
Automatically detects whether it is a fasta or normal text file.

@param path The location of the text file
@return Vector of sequences.
*/
std::vector<std::string> read_seq(const std::string& path); 

#endif /* __IZTLI_H__ */

 
