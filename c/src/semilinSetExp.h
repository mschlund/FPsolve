/*
 * semilinSetExp.h
 *
 *  Created on: 14.09.2012
 *      Author: maxi
 */

#ifndef SEMILINSETEXP_H_
#define SEMILINSETEXP_H_

#include "var.h"
#include <algorithm>

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


typedef std::pair<std::map<Var,unsigned int>, std::set<std::map<Var,unsigned int> > > LinSet;

// TODO: nothing is tested, nothing is done yet...
// FIXME: avoid deep-copying (+,*,..)
// TODO: refactoring needed!!...decouple semirings and semiring-elements! ???


// adding two maps componentwise... could be put in a util-class ?
std::map<Var,unsigned int> operator+(const std::map<Var, unsigned int>& a, const std::map<Var, unsigned int>& b)
{
	std::map<Var,unsigned int> result = a;
	for(std::map<Var, unsigned int>::const_iterator it = b.begin(); it!=b.end(); ++it) {
		if(result.find(it->first) != 0)
			result[it->first] += it->second;
		else
			result[it->first] = it->second;
	}

    return result;
}

std::ostream& operator<<(std::ostream& os, const std::map<Var, unsigned int>& v)
{
	os << "<";
    for (typename std::map<Var, unsigned int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        os << " " << it->second;
    }
    os << ">";
    return os;
}

// TODO: replace fixed-length vectors by a map: Var->int (multiplicity) ...

class SemilinSetExp : public CountingSemiring {
private:
	std::set<LinSet> val;
	static unsigned int alphabetSize;
	static std::map<Var, unsigned int> variables; //maps a variable to an index
	static unsigned int num_variables = 0;
	static std::set<LinSet>* elem_null = 0; // null = {} (empty set)
	static std::set<LinSet>* elem_one = 0; // one = {(0,0,...0)}


public:

	SemilinSetExp() {
		// not used normally :)
	}

	SemilinSetExp(const Var& var, int k) {
		if(variables.count(var) == 0) {
			variables.insert(std::make_pair(var,num_variables));
			num_variables++;
		}
		std::vector<unsigned int> offset = std::vector<unsigned int>();
	}

	SemilinSetExp(std::set<LinSet> val) {
		this->val = val;
	}

	static SemilinSetExp null() {
		if(!elem_null) {
			elem_null = new std::set<LinSet>();
		}
		return *elem_null;
	}

	static SemilinSetExp one() {
		if(!elem_one) {
			elem_one = new std::set<LinSet>();
			zeros = std::vector<unsigned int>(alphabetSize, 0);
			std::set<std::vector<unsigned int> > emptyset();
			elem_one->insert(std::make_pair(zeros, emptyset));
		}
		return *elem_one;
	}

	SemilinSetExp operator + (const SemilinSetExp& sl) {
		std::set<LinSet> result();
		std::set_union(this->val.begin(),this->val.end(),sl.val.begin(),sl.val.end(),result.begin());
		return SemilinSetExp(result);
	};

	SemilinSetExp operator * (const SemilinSetExp& sl) {
		std::set<LinSet> result;
		for(std::set<LinSet>::const_iterator it_arg = sl.val.begin(); it_arg != sl.val.end(); ++it_arg) {
			for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
						result.insert((*it_arg) * (*it_m));
			}
		}
	};

	SemilinSetExp star() {
		std::set<LinSet> result;
		for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
								result.insert(star(*it));
		}
		return SemilinearSetExp(result);
	};

	std::string string() const {
		std::stringstream ss;
		ss << "[" << std::endl;
		for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
			ss << it_m->first << " + ";
			for(std::set<std::vector<unsigned int> >::const_iterator it_gen = it_m->second.begin(); it_gen != it_m->second.end(); ++it_gen) {
				ss << *it_gen << " , ";
			}

			ss << std::endl;
		}
		ss << "]" << std::endl;
		return ss.str();
	};

	std::ostream& operator<<(std::ostream& os, const SemilinSetExp& slSet) {
		return os << slSet.string();
	}


	SemilinSetExp star(const LinSet& ls) {
		std::set<LinSet> result;
		// star of a linear set is a semilinear set:
		// (w_0.w_1*.w_2*...w_n*)* = 1 + \sum_{i=1}^{n-1} w_0^i.w_1*.w_2*...w_n* + w_0^n.w_0*.w_1*...w_n*
		result.insert(this->one());
		gens = ls.second;
		std::vector<unsigned int> tmp_offset = ls.first;

		int n = ls.second.size();
		for(int k=1; k<n; ++k) {
			result.insert(std::make_pair(tmp_offset,gens));
			for(int i =0; i<offset.size(); ++i) {
				tmp_offset[i] = tmp_offset[i] + ls.first[i];
			}
		}
		gens.insert(ls.first);
		result.insert(std::make_pair(tmp_offset,gens));

		return SemilinearSetExp(result);
	};

	static LinSet operator * (const LinSet& ls1, const LinSet& ls2) {
		std::set<std::vector<unsigned int> > generators();
		std::set_union(ls1.second.begin(),ls1.second.end(),ls2.second.begin(),ls2.second.end(), generators.begin());
		// add the offsets, union on the generators
		return std::make_pair(ls1.first+ls2.first, generators);
	};



};






#endif /* SEMILINSETEXP_H_ */
