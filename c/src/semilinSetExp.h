/*
 * semilinSetExp.h
 *
 *  Created on: 14.09.2012
 *      Author: maxi
 */

#ifndef SEMILINSETEXP_H_
#define SEMILINSETEXP_H_

#include "var.h"


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

// star{(v_00,{v_10,...,v_n0}),..., (v_0p,{v_1p,...,v_np})} = complicated...
// star(l_1,l_2,...,l_n) = \prod_{i=1}^n star(l_i)
// where star(l_i) gives a semilinear set

typedef std::pair<std::vector<unsigned int>, std::set<std::vector<unsigned int> > > LinSet;

//TODO: nothing is tested, nothing is done yet...


// adding two vectors componentwise... could be put in a util-class ?
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}


class SemilinSetExp : public CountingSemiring {
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

	SemilinSetExp star(const LinSet& ls)
	{
		// star of a linear set is a semilinear set:
		// (w_0.w_1*.w_2*...w_n*)* = 1 + \sum_{i=1}^{n-1} w_0^i.w_1*.w_2*...w_n* + w_0^n.w_0*.w_1*...w_n*

	};

	static LinSet operator * (const LinSet& ls1, const LinSet& ls2)
	{
		std::set<std::vector<unsigned int> > generators();
		std::set_union(ls1.second.begin(),ls1.second.end(),ls2.second.begin(),ls2.second.end(), generators.begin());
		// add the offsets, union on the generators
		return std::make_pair(ls1.first+ls2.first, generators);
	};




};






#endif /* SEMILINSETEXP_H_ */
