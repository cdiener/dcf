/*
 * pfam.cxx
 * This file is part of pfam_splitter
 *
 * Copyright (C) 2013 - Christian Diener
 *
 * pfam_splitter is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * pfam_splitter is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pfam_splitter. If not, see <http://www.gnu.org/licenses/>.
 */
 
#include <iostream>
#include <fstream>
#include <string>
#include "pfam.h"

std::vector<int> pfam_split(const char* path)
{
	std::vector<int> n_seqs;
	std::ifstream seqfile(path);

	int isfasta = (seqfile.peek() == '>');
	int fam_start, fam_stop;

	std::string tmp, family_new, family_old;
	family_new = family_old = "";
	std::ofstream out;

	if(!isfasta)
	{
		std::cerr<<"Not a valid fasta file!"<<std::endl;
		return n_seqs;
	}

	while(!seqfile.eof())
	{
		tmp = "";
		//Next sequence
		if(seqfile.peek() == '>')
		{
			std::getline(seqfile, tmp);
			fam_start = tmp.find(FAM_SEP, 0);
			fam_stop = tmp.find(FAM_SEP, fam_start+1);
			
			if(fam_start != fam_stop) //indeed both family separators found
			{
				family_new = tmp.substr(fam_start+1, fam_stop-fam_start-1);
				if(family_new != family_old)
				{
					if(out.is_open()) out.close();
					out.open(family_new+".fasta");
					family_old = family_new;
					n_seqs.push_back(1);
					if(n_seqs.size() % REPINT == 0) std::cout<<".";
					std::cout.flush();
				}
				
				out<<tmp<<std::endl;
			}
			n_seqs[n_seqs.size()-1]++;
		}
		else
		{
			std::getline(seqfile, tmp);
			out<<tmp<<std::endl;
		}
	}

	seqfile.close();
	out.close();

	return n_seqs;
}
