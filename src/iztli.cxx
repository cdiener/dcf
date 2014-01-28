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

double Hm(int* seq)
{
	double hwin, hm=0.0;
	int w = std::min(seq[0],W);
	
	for(unsigned int i=1; i<=seq[0]-w+1; i++)
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

	return hm/(seq[0]-w+1);
}

double charge(int* seq)
{
	int K,R,D,E;
	K=R=D=E=0;

	for(unsigned int i=1; i<=seq[0]; i++)
	{
		switch(seq[i])
		{
			case 18: ++K; break;
			case 19: ++R; break;
			case 14: ++D; break;
			case 16: ++E; break;
		}
	}

	return ( 1.0*(K+R-D-E) )/( seq[0] ); //1.0 for double conversion
}

double charge_var(int* seq)
{
	int K,R,D,E;
	K=R=D=E=0;

	double cwin, cmin, cmax;
	cmin = 1.0e32;
	cmax = 0.0; 

	int w = std::min(seq[0],W);

	for(unsigned int i=1; i<=seq[0]-w+1; i++)
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

double pI(int* seq)
{
	int C,D,E,Y,H,K,R;
	C=D=E=Y=H=K=R=0;

	for(unsigned int i=1; i<=seq[0]; i++)
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

double maxHM(int* seq)
{
	int w = std::min(seq[0],W);
	double Hcos, Hsin, Hup, Hlow;
	double* Hmax = new double[seq[0]-w+1];
	for(unsigned int i=0; i<seq[0]-w+1; i++) Hmax[i] = 0.0;
	Hcos = Hsin = 0.0;
	double* cos_t = new double[seq[0]]; 
	double* sin_t = new double[seq[0]];

	Hup = 0.0;
	Hlow = 1.0e32;	

	for(double angle=ANGMIN; angle<=ANGMAX; angle+=ANGINC)
	{
		Hcos = Hsin = 0.0;

		// We build up sine and cosine tables for faster calulation
		for(unsigned int i=0; i<seq[0]; i++)
		{
			cos_t[i] = cos(i*angle);

			//sin is calculated via the relation sin^2 + cos^2 = 1
			sin_t[i] = sqrt(1.0 - cos_t[i]*cos_t[i]);
			if( (int)(i*angle/M_PI) % 2 == 0) sin_t[i] = -sin_t[i];
		}
		
		//Calculate first window
		for(unsigned int i=1; i<=w; i++)
		{
			Hcos += SCALE[ seq[i] ]*cos_t[i-1];
			Hsin += SCALE[ seq[i] ]*sin_t[i-1];
		}
		Hmax[0] = std::max(Hmax[0], sqrt(Hcos*Hcos + Hsin*Hsin)/w);

		//Update only non-overlapping window parts
		for(unsigned int i=2; i<=seq[0]-w+1; i++)
		{
			Hcos += SCALE[ seq[i+w-1] ]*cos_t[i+w-2] - SCALE[ seq[i-1] ]*cos_t[i-2];
			Hsin += SCALE[ seq[i+w-1] ]*sin_t[i+w-2] - SCALE[ seq[i-1] ]*sin_t[i-2];

			Hmax[i-1] = std::max(Hmax[i-1], sqrt(Hcos*Hcos + Hsin*Hsin)/w);
		}
	}

	// Calculate mean and standard deviation 

	double delta;
	for(unsigned int i=0; i<seq[0]-w+1; i++)
	{
		Hup = std::max(Hup, Hmax[i]);
		Hlow = std::min(Hlow, Hmax[i]);
	}

	delete[] cos_t;
	delete[] sin_t;
	delete[] Hmax;

	return Hup-Hlow;
}

int* countAA(int* seq)
{
	int* counts = new int[20];
	for(unsigned int i=0; i<20; i++) counts[i] = 0;

	for(unsigned int i=1; i<=seq[0]; i++) counts[ seq[i] ]++;

	return counts;
}

double logP(int* seq, int* counts)
{
	double logP = logP_coef[0];

	for(unsigned int i=1; i<=20; i++) logP += logP_coef[i]*counts[i-1];

	return logP;
}

int* string_to_array(std::string seq)
{
	unsigned int n = seq.size();
	int* out = new int[n+1];
	out[0] = n;

	for(unsigned int i=1; i<=n; i++) out[i] = AAmap.find( seq[i-1] )->second;

	return out; 
}

std::vector<std::string> read_seq(std::string path)
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

 
 
