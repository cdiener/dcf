/*
 * pfam.h
 * This file is part of pfam_splitter
 *
 * Copyright (C) 2013 - Christian Diener
 *
 * Pfam_splitter is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Pfam_splitter is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pfam_splitter. If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef __PFAM_H__
#define __PFAM_H__

#include <vector>

/**
General prototypes and definitions.

Note that we encode the amoinoacids internally as integer vectors as this allows the sequence to be used as
a hash and makes most access operations O(1). For better usage the integer array begins with the number of elements
in the array and is immediately followed by the actual integer sequence.
*/

const char FAM_SEP[] = ";";
const int REPINT = 100; 


//Function prototypes

/**
Helper function to split a huge PFAM Fasta file into individual per-family files.

@param path The location of the fasta file
@return Vector of family identifiers
*/
std::vector<int> pfam_split(const char* path); 

#endif /* __PFAM_H__ */

 
