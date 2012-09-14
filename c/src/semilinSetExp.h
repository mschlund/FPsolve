/*
 * semilinSetExp.h
 *
 *  Created on: 14.09.2012
 *      Author: maxi
 */

#ifndef SEMILINSETEXP_H_
#define SEMILINSETEXP_H_

#include "var.h"

typedef std::pair<std::vector<unsigned int>, std::set<std::vector<unsigned int> > > LinSet;

// TODO: bigints instead of finite width ints ?

// use the identities
// (1) (x*)* = x*
// (2) (x+y)* = x*y*
// (3) (xy*)* = 1 + xx*y*   to push stars inwards
// use distributive law when multiplying


// represents the counting-SR over an alphabet of size k as (fully expanded) semilinear sets,
// i.e. a finite union of (shifted) (pointed) integer cones in \nat^k (=linear sets).
// each linear set is represented as a pair (v_0,{v_1,...,v_n}) with vectors v_i of natural numbers
// v_0 is usually called "offset", v_1,...v_n are the "generators" of the cone
// addition of two semilinear sets is just set-union A + B := A \cup B
// multiplication is componentwise addition of vectors A*B = {a+b| a\in A, b\in B} and implemented like this:
// {(v_00,{v_10,...,v_n0}),..., (v_0p,{v_1p,...,v_np})} * {(w_00,{w_10,...,w_n0}),..., (w_0q,{w_1q,...,w_nq})} =
// {(v_0i+w_0j,{v_1i,...,v_ni}\cup{w_1j,...,w_nj} ) | i=0..p, j=0..q}

// star{(v_00,{v_10,...,v_n0}),..., (v_0p,{v_1p,...,v_np})} =

class SemilinSetExp : public CountingSemiring
{
private:
	std::set<LinSet> Val;
	unsigned int alphabetSize;
	static std::set<LinSet>* elem_null; // null = {} (empty set)
	static std::set<LinSet>* elem_one; // one = {(0,0,...0)}

public:
	SemilinSetExp(int k) {
		alphabetSize = k;
	}

	SemilinSetExp() {

	}

	SemilinSetExp null() {
		if(elem_null)
			return elem_null;
		else {
			elem_null = new std::set<LinSet>();
		}
	}

	SemilinSetExp one() {

	}

	SemilinSetExp operator + (const SemilinSetExp& sl)
	{

	};

	SemilinSetExp operator * (const SemilinSetExp& sl)
	{

	};

	SemilinSetExp star()
	{

	};

	std::string string() const
	{

	};


};

#endif /* SEMILINSETEXP_H_ */
