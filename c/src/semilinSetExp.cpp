/*
 * semilinSetExp.cpp
 *
 *  Created on: 20.09.2012
 *      Author: maxi
 */

#include "semilinSetExp.h"

// adding two Var-maps componentwise... could be put in a util-class ?
VecSparse operator+(const VecSparse& a, const VecSparse& b) {
	VecSparse result = a;
	for(VecSparse::const_iterator it = b.begin(); it!=b.end(); ++it) {
		if(result.count(it->first) > 0)
			result[it->first] += it->second;
		else
			result[it->first] = it->second;
	}

    return result;
}

LinSet operator * (const LinSet& ls1, const LinSet& ls2) {
	//assert(!ls1.empty() && !ls2.empty());

	//union on the generators
	std::set<VecSparse> gens;
	std::insert_iterator<std::set<VecSparse> > it(gens, gens.begin());
	set_union(ls1.second.begin(), ls1.second.end(), ls2.second.begin(), ls2.second.end(), it);

	LinSet result;

	//add the offsets
	result.first = ls1.first + ls2.first;
	result.second = gens;
	return result;
}

std::ostream& operator<<(std::ostream& os, const VecSparse& v) {
	os << "<";
    for (typename std::map<VarPtr, unsigned int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        os << it->first << ":" << it->second << ", ";
    }
    os << ">";
    return os;
}

std::ostream& operator<<(std::ostream& os, const LinSet& ls) {
	os << ls.first;

	for(std::set<VecSparse>::const_iterator it_veclist = ls.second.begin(); it_veclist != ls.second.end(); ++it_veclist) {
		os << "+" << *it_veclist;
	}
    return os;
}

SemilinSetExp::SemilinSetExp() {
	this->val = std::set<LinSet>();
}

SemilinSetExp::SemilinSetExp(VarPtr var) {
	this->val = std::set<LinSet>();
	LinSet ls = LinSet();
	VecSparse offset;
	offset.insert(std::make_pair(var, 1));
	ls.first = offset;
	ls.second = std::set<VecSparse>();
	this->val.insert(ls);
}

SemilinSetExp::SemilinSetExp(std::set<LinSet> val) {
	this->val = std::set<LinSet>(val.begin(),val.end());
}

SemilinSetExp::~SemilinSetExp() {
	// do NOT delete static pointers!!!
}


SemilinSetExp SemilinSetExp::null() {
	if(!SemilinSetExp::elem_null)
		SemilinSetExp::elem_null = std::shared_ptr<SemilinSetExp>(new SemilinSetExp(std::set<LinSet>()));
	return *SemilinSetExp::elem_null;
}

SemilinSetExp SemilinSetExp::one() {
	if(!SemilinSetExp::elem_one) {
		std::set<LinSet> elone = std::set<LinSet>();
		LinSet ls = LinSet();
		ls.first = VecSparse();
		ls.second = std::set<VecSparse>();
		elone.insert(ls);
		SemilinSetExp::elem_one = std::shared_ptr<SemilinSetExp>(new SemilinSetExp(elone));
	}
	return *SemilinSetExp::elem_one;
}

SemilinSetExp SemilinSetExp::operator += (const SemilinSetExp& sl) {
	std::set<LinSet> result;

	std::insert_iterator<std::set<LinSet> > it(result, result.begin());
	std::set_union(this->val.begin(),this->val.end(),sl.getVal().begin(),sl.getVal().end(),it);

	*this = SemilinSetExp(result);
	return *this;
}

SemilinSetExp SemilinSetExp::operator *= (const SemilinSetExp& sl) {
	std::set<LinSet> result;
	for(std::set<LinSet>::const_iterator it_arg = sl.val.begin(); it_arg != sl.val.end(); ++it_arg) {
		for(std::set<LinSet>::const_iterator it_m = this->val.begin(); it_m != this->val.end(); ++it_m) {
					result.insert((*it_arg) * (*it_m));
		}
	}
	*this = SemilinSetExp(result);
	return *this;
}

std::set<LinSet> SemilinSetExp::getVal() const {
	return val;
}

//TODO: semantic equivalence check or at least some more sophisticated check
bool SemilinSetExp::operator == (const SemilinSetExp& sl) const {
	return (this->val == sl.getVal());
}

// TODO: multiple points of return might kill RVO ... try returning a pointer instead?
std::set<LinSet> SemilinSetExp::star(LinSet ls) {

	// if we do not have any generators, i.e. ls = w for some word w, just return w* (instead of 1 + ww*)
	if(ls.second.empty()) {
		LinSet r = ls;
		// if w is not the one-element, move w to the generators
		if(ls.first != VecSparse()) {
			r.second.insert(ls.first);
			r.first = VecSparse();
		}
		std::set<LinSet> res;
		res.insert(r);
		return res;
	}

	SemilinSetExp tmp_one = one();
	std::set<LinSet> v = tmp_one.getVal();
	std::set<LinSet> result = std::set<LinSet>(v.begin(),v.end());

	// star of a linear set is a semilinear set:
	// (w_0.w_1*.w_2*...w_n*)* = 1 +  (w_0.w_0*.w_1*.w_2*...w_n*)

	VecSparse new_gen = ls.first;
	LinSet ls_tmp = ls;
	ls_tmp.second.insert(new_gen);
	result.insert(ls_tmp);

	return result;
}

SemilinSetExp SemilinSetExp::star() const {
	SemilinSetExp result = SemilinSetExp::one();
	for(std::set<LinSet>::const_iterator it_m = val.begin(); it_m != val.end(); ++it_m) {
		std::set<LinSet> star_ls = star(*it_m);
		result = result * SemilinSetExp(star_ls);
	}
	return result;
}



std::string SemilinSetExp::string() const {
	std::stringstream ss;
	ss << "[" << '\n';
	for(std::set<LinSet>::const_iterator it_ls = this->val.begin(); it_ls != this->val.end(); ++it_ls) {
		ss << *it_ls;
		ss << '\n';
	}
	ss << "]" << '\n';
	return ss.str();
}


std::ostream& SemilinSetExp::operator<<(std::ostream& os) const {
	return os << this->string();
}


bool SemilinSetExp::is_idempotent = true;
bool SemilinSetExp::is_commutative = true;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_null;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_one;
