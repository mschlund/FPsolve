/*
 * why-set.h
 *
 *  Created on: 15.05.2014
 *      Author: schlund
 */

/*
 * Naive implementation of the why-semiring using sets
 */

#ifndef WHY_SET_H_
#define WHY_SET_H_

#include "semiring.h"
#include "../datastructs/var.h"
#include <set>

typedef std::set<VarId> VarSet;
typedef std::set<VarSet> WhySet;

class WhySemiring : public StarableSemiring<WhySemiring, Commutativity::Commutative, Idempotence::Idempotent>{

private:
  WhySet val;
  static std::shared_ptr<WhySemiring> elem_null;
  static std::shared_ptr<WhySemiring> elem_one;
public:
  WhySemiring();
  WhySemiring(VarId v);
  WhySemiring(std::string str_val);
  virtual ~WhySemiring();
  WhySemiring operator += (const WhySemiring& elem);
  WhySemiring operator *= (const WhySemiring& elem);
  bool operator == (const WhySemiring& elem) const;
  WhySemiring star () const;
  static WhySemiring null();
  static WhySemiring one();
  VarSet getVars() const;
  std::string string() const;
};



#endif /* WHY_SET_H_ */
