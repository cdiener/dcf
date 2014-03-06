/*
 * iztli.cxx
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
 * along with Iztli. If not, see <http://www.gnu.org/licenses/>.
 */

#include "iztli.h"
#include <iostream>
#include <fstream>

// Function implementations 

double Hm(const std::vector<int>& seq)
{
	double hwin, hm=0.0;
	int w = std::min( (int)seq.size() ,W);
	
	for(unsigned int i=0; i<=seq.size()-w; i++)
	{
		hwin = 0.0;

		//#pragma unroll 7
		for(unsigned int j=i; j<=i+w-1; j++)
		{ 	
			hwin += SCALE[ seq[j] ];
			if(SCALE[ seq[j] ] > 100.0) std::cout<<j<<std::endl;
		}
		hm += hwin/w; 
	}

	return hm/(seq.size()-w+1);
}

double charge(const std::vector<int>& seq)
{
	int K,R,D,E;
	K=R=D=E=0;

	for(unsigned int i=0; i<seq.size(); i++)
	{
		switch(seq[i])
		{
			case 18: ++K; break;
			case 19: ++R; break;
			case 14: ++D; break;
			case 16: ++E; break;
		}
	}

	return ( 1.0*(K+R-D-E) )/( seq.size() ); //1.0 for double conversion
}

double charge_var(const std::vector<int>& seq)
{
	int K,R,D,E;
	K=R=D=E=0;

	double cwin, cmin, cmax;
	cmin = 1.0e32;
	cmax = 0.0; 

	int w = std::min( (int)seq.size(),W);

	for(unsigned int i=0; i<=seq.size()-w; i++)
	{
		K=R=D=E=0;		
		
		//#pragma unroll 7
		for(unsigned int j=i; j<=i+w-1; j++)
		{ 	
			switch(seq[j])
			{	
				case 18: ++K; break;
				case 19: ++R; break;
				case 14: ++D; break;
				case 16: ++E; break;
			}
		}

		cwin = 1.0*(K+R-D-E)/w;
		cmin = std::min(cmin, cwin);
		cmax = std::max(cmax, cwin); 		
	}

	return cmax-cmin;
}

double pI(const std::vector<int>& seq)
{
	int C,D,E,Y,H,K,R;
	C=D=E=Y=H=K=R=0;

	for(unsigned int i=0; i<seq.size(); i++)
	{
		switch(seq[i])
		{
			case 11: ++C; break;
			case 14: ++D; break;
			case 16: ++E; break;
			case  6: ++Y; break;
			case 17: ++H; break;
			case 18: ++K; break;
			case 19: ++R; break;
		}
	}

	double pH = 6.5; //average pI
	double pHlow = 0.0;
	double pHup = 14.0;
	double eq = 0.0;
	double pH10 = 0.0;
	double eqlow = 0.0;
	double equp = 0.0;

	//bisection search for pI
	while( (pHup-pHlow)>TOLPI )
	{
		pH10 = pow(10.0, pH);
		eq = -1.0/(1.0+(pKCOOH/pH10)) - 1.0*C/(1.0+(pKC/pH10)) - 1.0*D/(1.0+(pKD/pH10)) - 1.0*E/(1.0+(pKE/pH10)) \
			- 1.0*Y/(1.0+(pKY/pH10)) + 1.0*H/(1.0+(pH10/pKH)) + 1.0*K/(1.0+(pH10/pKK)) + 1.0*R/(1.0+(pH10/pKR)) \
			+ 1.0/(1.0+(pH10/pKNH2));

		if(eq < 0) //pH too high
		{
			pHup = pH;
			pH = 0.5*(pH+pHlow);
			equp = eq;
		}
		else //pH too low
		{ 
			pHlow = pH;
			pH = 0.5*(pH+pHup);
			eqlow = eq;
		}
	}

	//Linear interpolation if we ever updated both boundaries
	if( (eqlow-equp)>0.0 ) return pHlow + eqlow*(pHup-pHlow)/(eqlow-equp);
	else return pH;
}

double rangeHM(const std::vector<int>& seq)
{
	int w = std::min( (int)seq.size() ,W);
	double Hcos, Hsin, Hup, Hlow;
	std::vector<double> Hmax = std::vector<double>(seq.size()-w+1, 0.0);
	Hcos = Hsin = 0.0;
	std::vector<double> cos_t(seq.size()); 
	std::vector<double> sin_t(seq.size());

	Hup = 0.0;
	Hlow = 1.0e32;	

	for(double angle=ANGMIN; angle<=ANGMAX; angle+=ANGINC)
	{
		Hcos = Hsin = 0.0;

		// We build up sine and cosine tables for faster calulation
		for(unsigned int i=0; i<seq.size(); i++)
		{
			cos_t[i] = cos(i*angle);

			//sin is calculated via the relation sin^2 + cos^2 = 1
			sin_t[i] = sqrt(1.0 - cos_t[i]*cos_t[i]);
			if( (int)(i*angle/M_PI) % 2 == 0) sin_t[i] = -sin_t[i];
		}
		
		//Calculate first window
		for(unsigned int i=0; i<w; i++)
		{
			Hcos += SCALE[ seq[i] ]*cos_t[i];
			Hsin += SCALE[ seq[i] ]*sin_t[i];
		}
		Hmax[0] = std::max(Hmax[0], sqrt(Hcos*Hcos + Hsin*Hsin)/w);

		//Update only non-overlapping window parts
		for(unsigned int i=1; i<=seq.size()-w; i++)
		{
			Hcos += SCALE[ seq[i+w-1] ]*cos_t[i+w-1] - SCALE[ seq[i-1] ]*cos_t[i-1];
			Hsin += SCALE[ seq[i+w-1] ]*sin_t[i+w-1] - SCALE[ seq[i-1] ]*sin_t[i-1];

			Hmax[i] = std::max(Hmax[i], sqrt(Hcos*Hcos + Hsin*Hsin)/w);
		}
	}

	// Calculate range 

	double delta;
	for(unsigned int i=0; i<seq.size()-w+1; i++)
	{
		Hup = std::max(Hup, Hmax[i]);
		Hlow = std::min(Hlow, Hmax[i]);
	}

	return Hup-Hlow;
}

std::vector<int> countAA(const std::vector<int>& seq)
{
	std::vector<int> counts(20, 0);
	//for(unsigned int i=0; i<20; i++) counts[i] = 0;

	for(unsigned int i=0; i<seq.size(); i++) counts[ seq[i] ]++;

	return counts;
}

double logP(const std::vector<int>& seq, const std::vector<int>& counts)
{
	double logP = logP_coef[0];

	for(unsigned int i=1; i<=20; i++) logP += logP_coef[i]*counts[i-1];

	return logP;
}

double alpha(const std::vector<int>& seq)
{
	unsigned int n = seq.size();
	if(n<7) return 0.0;
	
	int w = 5;
	std::vector<double> tetra_pep(n-3, 0.0);
	std::vector<int> alpha(n-3, 0);
	int alpha_sum; 
	
	// Get tetra-peptide propensities
	tetra_pep[0] = (alpha_prop[ seq[0] ] + alpha_prop[ seq[1] ] + alpha_prop[ seq[2] ] +
					alpha_prop[ seq[3] ])/4.0;
	for(unsigned int i=1; i<n-3; i++)
	{
		tetra_pep[i] = tetra_pep[i-1] - alpha_prop[ seq[i-1] ]/4.0 + alpha_prop[ seq[i+3] ]/4.0;
	}
	
	w = std::min(w, (int)tetra_pep.size() );
	//Check for alpha-fragments (4 out of 6 are alpha-prone)
	alpha_sum = 0;
	for(unsigned int j=0; j<w; j++) alpha_sum += (int)(tetra_pep[j] > alpha_cut);
	if(alpha_sum>4) for(unsigned int k=0; k<=w; k++) alpha[k] = 1;
	
	for(unsigned int i=1; i<tetra_pep.size()-w; i++)
	{
		alpha_sum = alpha_sum - (int)(tetra_pep[i-1] > alpha_cut) + (int)(tetra_pep[i+w] > alpha_cut);
		if(alpha_sum>4) for(unsigned int k=i; k<=i+w; k++) alpha[k] = 1;
	}
	
	// Extend alpha-sections
	int j=0;
	for(unsigned int i=1; i<alpha.size()-1; i++)
	{
		// Extend left
		if(alpha[i-1] == 0 && alpha[i] == 1)
		{
			j=i-1;
			while( j>=0 && alpha[j]==0 && (tetra_pep[j] > alpha_cut) )
			{
				alpha[j] = 1;
				j--;
			}
		}
		//Extend right
		else if( alpha[i] == 0 && alpha[i+1] == 1 )
		{
			j=i+1;
			while( j<alpha.size() && alpha[j]==0 && (tetra_pep[j] > alpha_cut) )
			{
				alpha[j] = 1;
				j++;
			}
		}
	}
	
	// Finally count the alpha-proportion
	alpha_sum = 0;
	for(unsigned int i=0; i<alpha.size(); i++) alpha_sum += alpha[i];
	
	return 1.0*alpha_sum/alpha.size();
}

std::vector<int> string_to_array(const std::string& seq)
{
	unsigned int n = seq.size();
	std::vector<int> out(n);

	for(unsigned int i=0; i<n; i++) out[i] = AAmap.find( seq[i] )->second;

	return out; 
}

std::vector<std::string> read_seq(const std::string& path)
{
	std::vector<std::string> seqs;	
	std::ifstream seqfile( path.c_str() );
	
	if(!seqfile.good())
	{
		std::cerr<<"Sequence file \""<<path<<"\" could not be openend!"<<std::endl;
		throw "error";
		return seqs;
	}

	int isfasta = (seqfile.peek() == '>');
	int i=-1;
	std::string tmp;

	if(isfasta)
	{
		while(!seqfile.eof())
		{
			tmp = "";
			//Next sequence
			if(seqfile.peek() == '>')
			{
				seqs.push_back( std::string() );
				std::getline(seqfile, tmp);				
				i++;
			}
			else
			{
				std::getline(seqfile, tmp);
				seqs[i] += tmp;
			}
		}
	}
	else
	{
		while(!seqfile.eof())
		{
			i++;
			std::getline(seqfile, tmp);
			if(tmp != "") seqs.push_back( std::string(tmp) );
		}
	}

	seqfile.close();

	return seqs;
}

 
 
