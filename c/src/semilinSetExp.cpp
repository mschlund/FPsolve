/*
 * semilinSetExp.cpp
 *
 *  Created on: 20.09.2012
 *      Author: maxi
 */

#include <algorithm>
#include <cassert>
#include <iostream>

#include "semilinSetExp.h"

// adding two Var-maps componentwise... could be put in a util-class ?
VecSparse operator+(const VecSparse &a, const VecSparse &b) {
  VecSparse result{a};
  for(auto &pair : b) {
    auto iter_bool = result.insert(pair);
    if (!iter_bool.second) {
      iter_bool.first->second += pair.second;
    }
  }
  return result;
}

/*
 * check if b = k*a for a natural number k.
 */
bool SemilinSetExp::divides(const VecSparse &a, const VecSparse &b) {
  if(a.size() != b.size())
    return false;

  unsigned int k=0;

  for(auto &pair : b) {
    if(a.count(pair.first) > 0) {
      // have we already obtained a candidate for a multiple?
      if(0 != k) {
        if(pair.second != k*a.at(pair.first) )
          return false;
      }
      // no candidate for k yet
      else if(pair.second % a.at(pair.first) != 0)
        return false;
      else
        // division leaves no remainder -> candidate k found
        k = pair.second / a.at(pair.first);
    }
    else {
      // Domains are different! => a does not divide b
      return false;
    }
  }
  return true;
}
/*
// check for inclusion between linsets via solving an ILP
bool SemilinSetExp::is_included(LinSet &ls1, LinSet &ls2) {

}

// simple inclusion check: just check if offset of ls1 is reachable in ls2 and if generators of ls1
// are included (set inclusion) in those of ls2
// even simpler: check if offsets are the SAME and if generators are subsets
bool SemilinSetExp::is_included_simple(LinSet &ls1, LinSet &ls2) {

}

// compute the (unique, minimal) Hilbert-basis for the linear set
void SemilinSetExp::min_hilbert(LinSet &ls) {

}

// check for pairwise inclusions between the linsets,
// if inclusion found: remove the subset
void SemilinSetExp::clean_slset() {

}
*/

/*
 * Perform simple optimization:
 * if g = k*h for some generators g!=h and natural number k, then remove g from generators
 */
void SemilinSetExp::clean_generators(LinSet &ls) {
  // check for all pairs of generators (g,h) if divides(h,g)
  for(auto &g : ls.second) {
    for(auto &h : ls.second) {
      if(g==h)
        continue;
      if(divides(h,g))
        ls.second.erase(g);
    }
  }
}




LinSet operator*(const LinSet &ls1, const LinSet &ls2) {

  LinSet result;

  /* Add the offsets... */
  result.first = ls1.first + ls2.first;

  /* ... and union on the generators */
  std::insert_iterator< std::set<VecSparse> > it(result.second,
                                                 result.second.begin());
  set_union(ls1.second.begin(), ls1.second.end(),
            ls2.second.begin(), ls2.second.end(), it);

  SemilinSetExp::clean_generators(result);
  return result;
}

std::ostream& operator<<(std::ostream &os, const VecSparse &v) {
  os << "<";
  for (auto &pair : v) {
    os << pair.first << ":" << pair.second << ", ";
  }
  os << ">";
  return os;
}

std::ostream& operator<<(std::ostream &os, const LinSet &ls) {
  os << ls.first;
  for (auto &veclist : ls.second) {
    os << "+" << veclist;
  }
  return os;
}

SemilinSetExp::SemilinSetExp() : val() { }

SemilinSetExp::SemilinSetExp(VarPtr var) : val() {
  VecSparse offset = { std::make_pair(var, 1) };
  LinSet ls{};
  ls.first = offset;
  val.insert(std::move(ls));
}

SemilinSetExp::SemilinSetExp(VarPtr var, unsigned int cnt) : val() {
  if(0 != cnt) {
    VecSparse offset = { std::make_pair(var, cnt) };
    LinSet ls{};
    ls.first = offset;
    val.insert(std::move(ls));
  }
  else {
    std::cerr << "[INFO] SL-set: tried to generate slset having a variable with count zero.. ignoring it." << std::endl;
  }
}

SemilinSetExp::SemilinSetExp(const std::set<LinSet> &v) {
  val = v;
}

SemilinSetExp::~SemilinSetExp() {
  // do NOT delete static pointers!!!
}


SemilinSetExp SemilinSetExp::null() {
  if(!SemilinSetExp::elem_null) {
    SemilinSetExp::elem_null =
      std::make_shared<SemilinSetExp>(std::set<LinSet>());
  }
  return *SemilinSetExp::elem_null;
}

SemilinSetExp SemilinSetExp::one() {
  if (!SemilinSetExp::elem_one) {
    std::set<LinSet> elone = { LinSet{} };
    SemilinSetExp::elem_one = std::make_shared<SemilinSetExp>(elone);
  }
  return *SemilinSetExp::elem_one;
}

// TODO: check for obvious inclusions and remove them
SemilinSetExp SemilinSetExp::operator+=(const SemilinSetExp &sl) {
  std::set<LinSet> result;
  std::insert_iterator< std::set<LinSet> > it(result, result.begin());
  std::set_union(val.begin(), val.end(), sl.val.begin(), sl.val.end(), it);
  val = std::move(result);
  return *this;
}

SemilinSetExp SemilinSetExp::operator*=(const SemilinSetExp &sl) {
  std::set<LinSet> result;
  for(auto &lin_set_rhs : sl.val) {
    for(auto &lin_set_lhs : val) {
      result.insert(lin_set_rhs * lin_set_lhs);
    }
  }
  val = std::move(result);
  return *this;
}

// TODO: semantic equivalence check or at least some more sophisticated check
bool SemilinSetExp::operator == (const SemilinSetExp &sl) const {
  return (val == sl.val);
}

std::set<LinSet> SemilinSetExp::star(const LinSet &ls) {

  /* If we do not have any generators, i.e.,
   *   ls = w  (for some word w)
   * just return
   *   w*
   * instead of 1 + ww*.  */
  if (ls.second.empty()) {
    LinSet r = ls;
    /* If w is not the one-element, move w to the generators. */
    if (ls.first != VecSparse()) {
      r.second.insert(ls.first);
      r.first = VecSparse();
    }
    std::set<LinSet> res = { r };
    return res;
  }

  SemilinSetExp tmp_one{one()};
  std::set<LinSet> v{tmp_one.val};
  std::set<LinSet> result{v};

  /* Star of a linear set is a semilinear set:
   * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*)
   */

  LinSet ls_tmp = ls;
  ls_tmp.second.insert(ls.first);
  result.insert(std::move(ls_tmp));

  return result;
}

SemilinSetExp SemilinSetExp::star() const {
  SemilinSetExp result = SemilinSetExp::one();
  for (auto &ls : val) {
    result *= SemilinSetExp(star(ls));
  }
  return result;
}

std::string SemilinSetExp::string() const {
  std::stringstream ss;
  ss << "[";
  for (auto &ls : val) {
    ss << ls << "; ";
  }
  ss << "]";
  return ss.str();
}

std::ostream& SemilinSetExp::operator<<(std::ostream &os) const {
  return os << string();
}


const bool SemilinSetExp::is_idempotent = true;
const bool SemilinSetExp::is_commutative = true;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_null;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_one;
