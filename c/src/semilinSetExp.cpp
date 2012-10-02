/*
 * semilinSetExp.cpp
 *
 *  Created on: 20.09.2012
 *      Author: maxi
 */

#include "semilinSetExp.h"



// adding two Var-maps componentwise... could be put in a util-class ?
VecSparse operator+(const VecSparse& a, const VecSparse& b)
{
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
	assert(!ls1.empty() && !ls2.empty());

	//union on the generators
	std::list<VecSparse> g1 = ls1;
	std::list<VecSparse> g2 = ls2;

	//remove the offsets
	g1.pop_front();
	g2.pop_front();

	g1.sort();
	g2.sort();

	std::list<VecSparse> result;
	std::insert_iterator<std::list<VecSparse> > it(result, result.begin());
	set_union(g1.begin(), g1.end(), g2.begin(), g2.end(), it);

	//add the offsets
	result.push_front(ls1.front()+ls2.front());
	return result;

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

std::ostream& operator<<(std::ostream& os, const LinSet& ls) {
	for(std::list<VecSparse>::const_iterator it_veclist = ls.begin(); it_veclist != ls.end(); ++it_veclist) {
		os << *it_veclist << " + ";
	}
    return os;
}

SemilinSetExp::SemilinSetExp() {
	this->val = std::set<LinSet>();
}

SemilinSetExp::SemilinSetExp(Var var) {
	this->val = std::set<LinSet>();
	LinSet ls = LinSet();
	VecSparse offset;
	offset.insert(std::make_pair(var, 1));
	ls.push_front(offset);
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
		SemilinSetExp::elem_null = new SemilinSetExp(std::set<LinSet>());
	return *SemilinSetExp::elem_null;
}

SemilinSetExp SemilinSetExp::one() {
	if(!SemilinSetExp::elem_one) {
		std::set<LinSet> elone = std::set<LinSet>();
		LinSet ls = LinSet();
		ls.push_front(VecSparse());
		elone.insert(ls);
		SemilinSetExp::elem_one = new SemilinSetExp(elone);
	}
	return *SemilinSetExp::elem_one;
}

SemilinSetExp SemilinSetExp::operator + (const SemilinSetExp& sl) const {
	std::set<LinSet> result;

	for(std::set<LinSet>::const_iterator it_arg = sl.val.begin(); it_arg != sl.val.end(); ++it_arg) {
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
	for(std::set<LinSet>::const_iterator it_arg = sl.val.begin(); it_arg != sl.val.end(); ++it_arg) {
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

// TODO: multiple points of return might kill RVO ... try returning a pointer instead!
std::set<LinSet> SemilinSetExp::star(LinSet ls) {

	// if we do not have generators, i.e. ls = w for some word w, just return w* (instead of 1 + ww*)
	if(ls.size() == 1) {
		LinSet r = ls;
		r.push_front(VecSparse());

		std::set<LinSet> res;
		res.insert(r);
		return res;
	}

	SemilinSetExp tmp_one = one();
	std::set<LinSet> v = tmp_one.getVal();
	std::set<LinSet> result = std::set<LinSet>(v.begin(),v.end());

	// star of a linear set is a semilinear set:
	// (w_0.w_1*.w_2*...w_n*)* = 1 + \sum_{i=1}^{n-1} w_0^i.w_1*.w_2*...w_n* + w_0^n.w_0*.w_1*...w_n*


	// for all elements of gens:
	//   for all keys of tmp_offset
	//      tmp_offset[key] = tmp_offset[key] + ls.first[key]

	LinSet ls_tmp = ls;
	VecSparse& tmp_offset = ls_tmp.front();


	int n = ls.size() - 1;

	//for(std::list<VecSparse>::iterator it_gens = gen_ls; it_gens!=ls.end(); ++it_gens) {
	for(int i=1; i<n; i++) {
		result.insert(ls_tmp);
		// calculate w_0^i
		for(VecSparse::iterator it_tmp_offset = tmp_offset.begin(); it_tmp_offset != tmp_offset.end(); ++it_tmp_offset) {
			it_tmp_offset->second = it_tmp_offset->second + ls.front()[it_tmp_offset->first];
		}
	}

	//the last summand is w_0^n.w_0*.w_1*...w_n*
	VecSparse new_gen = ls.front();
	std::list<VecSparse>::iterator gen_ls = ls_tmp.begin();
	++gen_ls;
	ls_tmp.insert(gen_ls, new_gen);
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
};


std::ostream& SemilinSetExp::operator<<(std::ostream& os) const {
	return os << this->string();
}


bool SemilinSetExp::is_idempotent = true;
bool SemilinSetExp::is_commutative = true;
SemilinSetExp* SemilinSetExp::elem_null = 0;
SemilinSetExp* SemilinSetExp::elem_one = 0;
