/*
 * semilinSetExp.h
 *
 *  Created on: 14.09.2012
 *      Author: maxi
 */

#ifndef SEMILINSETEXP_H_
#define SEMILINSETEXP_H_

#include "var.h"
#include "counting-semiring.h"

#include <cassert>

#include <algorithm>
#include <set>
#include <map>
#include <list>

// TODO: bigints instead of finite precision ints ?

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

// star{(v_00,{v_10,...,v_n0}),..., (v_0p,{v_1p,...,v_np})} = complicated...
// star(l_1,l_2,...,l_n) = \prod_{i=1}^n star(l_i)
// where star(l_i) gives a semilinear set

typedef std::map<VarPtr,unsigned int> VecSparse;
//typedef std::pair<VecSparse, std::set<VecSparse>  > LinSet;
typedef std::list<VecSparse> LinSet; // new implementation..

// FIXME: avoid deep-copying (+,*,..)
// TODO: refactoring ?? ... decouple semirings and semiring-elements! ???

class SemilinSetExp : public Semiring<SemilinSetExp> {
private:
	std::set<LinSet> val;

public:
	static SemilinSetExp* elem_null; // null = {} (empty set)
	static SemilinSetExp* elem_one; // one = {(0,0,...0)}

	SemilinSetExp();
	SemilinSetExp(std::set<LinSet> val);
	SemilinSetExp(VarPtr v);
	virtual ~SemilinSetExp();
	static SemilinSetExp null();
	static SemilinSetExp one();
	SemilinSetExp operator + (const SemilinSetExp& sl) const;
	SemilinSetExp operator * (const SemilinSetExp& sl) const;
	static std::set<LinSet> star(LinSet ls);
	SemilinSetExp star() const;
	std::string string() const;
	bool operator == (const SemilinSetExp& sl) const;
	std::set<LinSet> getVal() const;
	std::ostream& operator<<(std::ostream& os) const;

	static bool is_idempotent;
	static bool is_commutative;
};






#endif /* SEMILINSETEXP_H_ */
