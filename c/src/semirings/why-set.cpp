/*
 * why-set.cpp
 *
 *  Created on: 15.05.2014
 *      Author: maxi
 */

#include <algorithm>
#include "why-set.h"
#include "../utils/profiling-macros.h"
#include "../utils/string_util.h"

// std::map in polynomial.h wants this constructor...
WhySemiring::WhySemiring()
{
  this->val = WhySet();
}

WhySemiring::WhySemiring(VarId v)
{
  this->val = WhySet({VarSet({v})});
}

WhySemiring::WhySemiring(std::string str_val)
{
  if(0 == str_val.compare("1")) {
    this->val = WhySet();
    val.insert(VarSet());
  }
  else if (0 == str_val.compare("0")) {
    this->val = WhySet();
  }
  else {
    this->val = WhySet();
    VarSet tmp = VarSet();
    tmp.insert(Var::GetVarId(str_val));
    val.emplace(tmp);
  }
}

WhySemiring::~WhySemiring()
{

}

WhySemiring WhySemiring::operator+=(const WhySemiring& elem)
{
  OPADD;
  this->val.insert(elem.val.begin(),elem.val.end());
  return *this;
}

WhySemiring WhySemiring::operator*=(const WhySemiring& elem)
{
  OPMULT;
  if(this->val.size()==0 || elem.val.size()==0) {
    *this = null();
  }
  else {
    WhySet res;
    for (auto vs : this->val) {
      for (auto es : elem.val) {
        auto tmp = vs;
        tmp.insert(es.begin(),es.end());
        res.insert(tmp);
        }
    }
    this->val = res;
  }
  return *this;
}

bool WhySemiring::operator==(const WhySemiring& elem) const
{
  return this->val == elem.val;
}

WhySemiring WhySemiring::star() const
{
  OPSTAR;
  WhySemiring tmp = *this;
  std::uint_fast16_t N = this->getVars().size();
  return (one() + pow(tmp,N));
}

WhySemiring WhySemiring::null()
{
  if(!WhySemiring::elem_null)
    WhySemiring::elem_null = std::shared_ptr<WhySemiring>(new WhySemiring("0"));
  return *WhySemiring::elem_null;
}

WhySemiring WhySemiring::one()
{
  if(!WhySemiring::elem_one)
    WhySemiring::elem_one = std::shared_ptr<WhySemiring>(new WhySemiring("1"));
  return *WhySemiring::elem_one;
}

std::string WhySemiring::string() const
{
  std::stringstream ss;
  ss << (this->val);
  return ss.str();
}


VarSet WhySemiring::getVars() const {
  VarSet res;
  for (auto ws : val) {
    res.insert(ws.begin(), ws.end());
  }
  return res;
}

std::shared_ptr<WhySemiring> WhySemiring::elem_null;
std::shared_ptr<WhySemiring> WhySemiring::elem_one;
