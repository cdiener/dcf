/*
 * predict.cxx
 * This file is part of dcf
 *
 * Copyright (C) 2016 - Christian Diener
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
 
#include <iostream>
#include <fstream>
#include <omp.h>
#include "iztli.h"

/**
prop program

@section DESCRIPTION

Calculates the physiochemical properties for several peptide sequences.
 
Usage:
To get the properties for a set of txt or fasta files in csv use
@code{.sh}
> prop [FILES]
@endcode 
*/

using namespace std;

const string usage = "Usage: prop [FILES]\n"
    "Calculate some basic physio-chemical properties for peptide sequences.\n"
    "Calculation will take in parallel and, thus, will be fast.\n"
    "You can pass as many files as you want in TXT (one sequence per line)\n"
    "or FASTA format. Output will be in CSV format.";
const int n_var = 27;

int main (int argc, char* argv[])
{	
	// Check command line arguments
	if(argc<2) {		
		cerr<<"Missing arguments! Need at least one file!"<<endl;
		cout<<usage<<endl;
		return 1;
	}
    string farg = argv[1];
    
    if(farg == "-h" || farg == "--help") {
        cout<<usage<<endl;
        return 0;
    } 
    
    double start, end, time_s;
	start = omp_get_wtime();
	
    int nseq = 0;
    
    #pragma omp parallel for
    for(int i=1; i<argc; ++i) {
        vector<string> seqs = read_seq(argv[i]);
        
        #pragma omp critical 
        {
            nseq += seqs.size();
        }
        
        string out_name = argv[i];
        string::size_type const p(out_name.find_last_of('.'));
        out_name = out_name.substr(0, p);
        out_name += "_props.csv";
        ofstream out(out_name);
        out<<"sequence,charge,hydrophobicity,pI,hydrophobic_moment_range,charge_range,logP,fraction_alpha"; 
		for(unsigned int k=0; k<20; k++) out<<","<<AA[k];
        out<<endl;
        
        #pragma parallel for
        for(unsigned int k=0; k<seqs.size(); k++) {
            if(seqs[k].size()==0) continue;
            
            vector<int> iseq = string_to_array(seqs[k]);
            vector<double> props(n_var);
            
            props[0] = charge(iseq);
            props[1] = Hm(iseq);
            props[2] = pI(iseq);
            props[3] = rangeHM(iseq);
            props[4] = charge_var(iseq);
            vector<int> counts = countAA(iseq);
            props[5] = logP(counts);
            props[6] = alpha(iseq);
            
            out<<seqs[k];
            for(int pi=0; pi<7; ++pi) out<<","<<props[pi];
            for(int ci=0; ci<counts.size(); ++ci) 
                out<<","<<1.0*counts[ci]/iseq.size();
            out<<endl;
        }
        
        out.close();
        
    }
    
	end = omp_get_wtime();
	time_s = end-start;
	cout<<"Needed "<<time_s*1.0e3<<" ms (";
    cout<<1.0*time_s/nseq*1e6<<" \xC2\xB5s per peptide)."<<endl;


	return 0;
}

  
 
