/*
 * main.cxx
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

#include <chrono>
#include <iostream>
#include "pfam.h"

using namespace std;

int main (int argc, char const* argv[])
{
	if(argc!=2)
	{	
		cout<<"Usage: ./pfam_splitter fasta_file"<<endl;
		return 1;
	}	

	cout<<"Parsing";
	auto start = chrono::high_resolution_clock::now();

	vector<int> fams = pfam_split(argv[1]);
	auto end = chrono::high_resolution_clock::now();	
	double time_m = chrono::duration<double, milli>(end-start).count()/60000.0;
	cout<<endl;

	for(unsigned int i=1; i<fams.size(); i++) fams[0] += fams[i];

	cout<<"Parsed "<<fams[0]<<" peptides in "<<fams.size()<<" families in "<<time_m<<" minutes."<<endl;


	return 0;
}
