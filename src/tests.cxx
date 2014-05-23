/*
 * tests.cxx
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
 
#define CATCH_CONFIG_MAIN
#include "catch.h"

/**
 * Here we will run all the unit tests for individual parts of the program.
 */

#include "iztli.h"
#include "trainer.h"
#include "design.h"
#include "alglib/ap.h"

// c++ includes
#include <string>
#include <vector>

using namespace std;

TEST_CASE("iztli file reader works", "[iztli]")
{
	vector<string> seqs = read_seq("examples/pos_CPP.txt");
	
	unsigned int good = 0;
	for(unsigned int i=0; i<seqs.size(); i++) good += (seqs[i].size() > 0);
	REQUIRE( seqs.size() > 0 );
	REQUIRE( good == seqs.size() ); 
}

TEST_CASE("iztli works", "[iztli]")
{
	const string seq = "IIIRRR";
	const vector<int> real_iseq = {0,0,0,19,19,19};
	vector<int> real_counts = vector<int>(20, 0);
	real_counts[0] = 3;
	real_counts[19] = 3; 
	
	const vector<int> iseq = string_to_array(seq);
	REQUIRE( iseq == real_iseq );	
	
	
	// value tests
	REQUIRE( Hm(iseq)>= -10.0 ); 
	REQUIRE( Hm(iseq)<=10.0 );
	REQUIRE( abs(charge(iseq) - 0.5) < 1.0e-6 );
	REQUIRE( charge_var(iseq)>=0.0 ); 
	REQUIRE( charge_var(iseq)<=1.0 );
	REQUIRE( pI(iseq)>=0.0 );
	REQUIRE( pI(iseq)<=14.0 ); 
	REQUIRE( rangeHM(iseq)>=0.0 );
	vector<int> counts = countAA(iseq);
	REQUIRE( counts == real_counts );
	REQUIRE( abs(logP(counts)) > 0.0 );
	REQUIRE( alpha(iseq) < 1.0e-6 );
}


