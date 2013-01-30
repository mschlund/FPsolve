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

LinSet operator*(const LinSet &ls1, const LinSet &ls2) {
  //assert(!ls1.empty() && !ls2.empty());

  LinSet result;

  /* Add the offsets... */
  result.first = ls1.first + ls2.first;

  /* ... and union on the generators */
  std::insert_iterator< std::set<VecSparse> > it(result.second,
                                                 result.second.begin());
  set_union(ls1.second.begin(), ls1.second.end(),
            ls2.second.begin(), ls2.second.end(), it);

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
  VecSparse offset = { std::make_pair(var, cnt) };
  LinSet ls{};
  ls.first = offset;
  val.insert(std::move(ls));
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

std::set<LinSet> SemilinSetExp::getVal() const {
  return val;
}

// TODO: semantic equivalence check or at least some more sophisticated check
bool SemilinSetExp::operator == (const SemilinSetExp &sl) const {
  return (val == sl.getVal());
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
  std::set<LinSet> v{tmp_one.getVal()};
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
    result = result * SemilinSetExp(star(ls));
  }
  return result;
}

std::string SemilinSetExp::string() const {
  std::stringstream ss;
  ss << "[" << std::endl;
  for (auto &ls : val) {
    ss << ls << std::endl;
  }
  ss << "]" << std::endl;
  return ss.str();
}

std::ostream& SemilinSetExp::operator<<(std::ostream &os) const {
  return os << string();
}


const bool SemilinSetExp::is_idempotent = true;
const bool SemilinSetExp::is_commutative = true;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_null;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_one;
