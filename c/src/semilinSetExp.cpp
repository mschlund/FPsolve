/*
 * semilinSetExp.cpp
 *
 *  Created on: 20.09.2012
 *      Author: maxi
 */

#include "semilinSetExp.h"



// adding two Var-maps componentwise... could be put in a util-class ?
std::map<Var,unsigned int> operator+(const std::map<Var, unsigned int>& a, const std::map<Var, unsigned int>& b)
{
	std::map<Var,unsigned int> result = a;
	for(VecSparse::const_iterator it = b.begin(); it!=b.end(); ++it) {
		if(result.count(it->first) > 0)
			result[it->first] += it->second;
		else
			result[it->first] = it->second;
	}

    return result;
}

LinSet operator * (LinSet ls1, LinSet ls2) {
	std::set<VecSparse> generators;

	//std::set_union(ls1.second.begin(),ls1.second.end(),ls2.second.begin(),ls2.second.end(), generators.begin());
	for(std::set<VecSparse>::const_iterator it_1=ls1.second.begin(); it_1 != ls1.second.end(); ++it_1) {
		generators.insert(*it_1);
	}
	for(std::set<VecSparse>::const_iterator it_2=ls2.second.begin(); it_2 != ls2.second.end(); ++it_2) {
		generators.insert(*it_2);
	}

	// add the offsets, union on the generators
	return std::make_pair(ls1.first+ls2.first, generators);
}

std::ostream& operator<<(std::ostream& os, const VecSparse& v) {
	os << "<";
    for (typename std::map<Var, unsigned int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        os << it->first.string() << ":" << it->second << ", ";
    }
    os << ">";
    return os;
}



SemilinSetExp::SemilinSetExp() {
	this->val = std::set<LinSet>();
}

SemilinSetExp::SemilinSetExp(Var var) {
	this->val = std::set<LinSet>();
	std::map<Var,unsigned int> offset;
	offset.insert(std::make_pair(var, 1));
	std::set<std::map<Var,unsigned int> > gens;
	LinSet l = std::make_pair(offset,gens);
	this->val.insert(l);
}

SemilinSetExp::SemilinSetExp(std::set<LinSet> val) {
	this->val = val;
}

SemilinSetExp::~SemilinSetExp() {
	if(SemilinSetExp::elem_null)
		delete SemilinSetExp::elem_null;
	if(SemilinSetExp::elem_one)
		delete SemilinSetExp::elem_one;
}


SemilinSetExp SemilinSetExp::null() {
	if(!SemilinSetExp::elem_null)
		SemilinSetExp::elem_null = new SemilinSetExp(std::set<LinSet>());
	return *SemilinSetExp::elem_null;
}

SemilinSetExp SemilinSetExp::one() {
	if(!SemilinSetExp::elem_one) {
		std::set<LinSet> elone;
		VecSparse zeros;
		std::set<VecSparse> emptyset;
		elone.insert(std::make_pair(zeros, emptyset));
		SemilinSetExp::elem_one = new SemilinSetExp(elone);
	}
	return *SemilinSetExp::elem_one;
}

SemilinSetExp SemilinSetExp::operator + (const SemilinSetExp& sl) const {
	std::set<LinSet> result;

	for(std::set<LinSet>::const_iterator it_arg = sl.getVal().begin(); it_arg != sl.val.end(); ++it_arg) {
			result.insert(*it_arg);
	}
	for(std::set<LinSet>::const_iterator it_this = val.begin(); it_this != val.end(); ++it_this) {
			result.insert(*it_this);
	}

	//std::set_union(this->val.begin(),this->val.end(),sl.getVal().begin(),sl.getVal().end(),result.begin());

	return SemilinSetExp(result);
}

SemilinSetExp SemilinSetExp::operator * (const SemilinSetExp& sl) const {
	std::set<LinSet> result;
	for(std::set<LinSet>::const_iterator it_arg = sl.getVal().begin(); it_arg != sl.val.end(); ++it_arg) {
		for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
					result.insert((*it_arg) * (*it_m));
		}
	}
	return SemilinSetExp(result);
}

std::set<LinSet> SemilinSetExp::getVal() const {
	return val;
}

bool SemilinSetExp::operator == (const SemilinSetExp& sl) const {
	return (this->val == sl.getVal());
}

SemilinSetExp SemilinSetExp::star(LinSet ls) {
	std::set<LinSet> result = one().getVal();
	// star of a linear set is a semilinear set:
	// (w_0.w_1*.w_2*...w_n*)* = 1 + \sum_{i=1}^{n-1} w_0^i.w_1*.w_2*...w_n* + w_0^n.w_0*.w_1*...w_n*
	std::set<VecSparse> gens = ls.second;
	VecSparse tmp_offset = ls.first;

	//int n = ls.second.size();

	// for all elements of gens:
	//   for all keys of tmp_offset
	//      tmp_offset[key] = tmp_offset[key] + ls.first[key]
	for(std::set<VecSparse>::const_iterator it_gen = gens.begin(); it_gen!=gens.end(); ++it_gen) {
		result.insert(std::make_pair(tmp_offset,gens));
		// calculate w_0^i
		for(VecSparse::iterator it_tmp_offset = tmp_offset.begin(); it_tmp_offset != tmp_offset.end(); ++it_tmp_offset) {
			it_tmp_offset->second = it_tmp_offset->second + ls.first[it_tmp_offset->first];
		}
	}
	//the last summand is w_0^n.w_0*.w_1*...w_n*
	gens.insert(ls.first);
	result.insert(std::make_pair(tmp_offset,gens));

	return SemilinSetExp(result);
}

SemilinSetExp SemilinSetExp::star() const {
	std::set<LinSet> result;
	for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
		std::set<LinSet> star_linset = star(*it_m).getVal();
		result.insert(star_linset.begin(), star_linset.end());
	}
	return SemilinSetExp(result);
}

std::string SemilinSetExp::string() const {
	std::stringstream ss;
	ss << "[" << std::endl;
	for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
		ss << it_m->first << " + ";
		for(std::set<std::map<Var, unsigned int> >::const_iterator it_gen = it_m->second.begin(); it_gen != it_m->second.end(); ++it_gen) {
			ss << *it_gen << " , ";
		}

		ss << std::endl;
	}
	ss << "]" << std::endl;
	return ss.str();
};

std::ostream& SemilinSetExp::operator<<(std::ostream& os) const {
	return os << this->string();
}


bool SemilinSetExp::is_idempotent = true;
bool SemilinSetExp::is_commutative = true;
SemilinSetExp* SemilinSetExp::elem_null = 0;
SemilinSetExp* SemilinSetExp::elem_one = 0;
