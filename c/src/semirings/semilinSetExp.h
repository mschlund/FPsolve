/*
 * semilinSetExp.h
 *
 *  Created on: 14.09.2012
 *      Author: maxi
 */

#ifndef SEMILINSETEXP_H_
#define SEMILINSETEXP_H_

#include <map>
#include <memory>
#include <set>

#include "../datastructs/var.h"
#include "semiring.h"

/* TODO:
 * - Bigints instead of finite precision ints?  Maybe should be configurable
 *   (e.g. templates)?
 * - Use the identities:
 *   (1) (x*)* = x*
 *   (2) (x+y)* = x*y*
 *   (3) (xy*)* = 1 + xx*y*  (to push stars inwards)
 * - Use distributive law when multiplying.
 */


// represents the counting-SR over an alphabet of size k as (fully expanded) semilinear sets,
// i.e. a finite union of (shifted) (pointed) integer cones in \nat^k (=linear sets).
// each linear set is represented as a pair (v_0,{v_1,...,v_n}) with vectors v_i of natural numbers
// v_0 is usually called "offset", v_1,...v_n are the "generators" of the cone
// addition of two semilinear sets is just set-union A + B := A \cup B
// multiplication is componentwise addition of vectors A*B = {a+b| a\in A, b\in B} and implemented like this:
// {(v_00,{v_10,...,v_n0}),..., (v_0p,{v_1p,...,v_np})} * {(w_00,{w_10,...,w_n0}),..., (w_0q,{w_1q,...,w_nq})} =
// {(v_0i+w_0j,{v_1i,...,v_ni}\cup{w_1j,...,w_nj} ) | i=0..p, j=0..q}

// star{(v_00,{v_10,...,v_n0}),..., (v_0p,{v_1p,...,v_np})}:
// star(l_1,l_2,...,l_n) = \prod_{i=1}^n star(l_i)
// where star((v_00,{v_10,...,v_n0})) = 1 +  (v_00,{v_00,v_10,...,v_n0})


typedef std::map<VarId, unsigned int> VecSparse;

/* At some point we used
 *   typedef std::list<VecSparse> LinSet;
 * but the current one turns out to be more convenient.  Reasons: set takes care
 * of uniqueness of generators, sorting sets is trivial...
 */
typedef std::pair< VecSparse, std::set<VecSparse> > LinSet;


class SemilinSetExp : public StarableSemiring<SemilinSetExp, Commutativity::Commutative, Idempotence::Idempotent> {
  private:
    std::set<LinSet> val;
    /* null = {} (empty set) */
    static std::shared_ptr<SemilinSetExp> elem_null;
    /* one = (0, {(0,0,...0)}) */
    static std::shared_ptr<SemilinSetExp> elem_one;

  public:
    SemilinSetExp();
    SemilinSetExp(const std::set<LinSet> &val);
    SemilinSetExp(VarId v);
    SemilinSetExp(VarId var, unsigned int cnt);

    virtual ~SemilinSetExp();

    static SemilinSetExp null();
    static SemilinSetExp one();
    static std::set<LinSet> star(const LinSet &ls);


    SemilinSetExp star() const;
    std::string string() const;

    SemilinSetExp operator += (const SemilinSetExp& sl);
    SemilinSetExp operator *= (const SemilinSetExp& sl);

    bool operator == (const SemilinSetExp& sl) const;
    std::ostream& operator<<(std::ostream& os) const;

    void clean_slset();

    static const bool is_idempotent;
    static const bool is_commutative;
};





#endif /* SEMILINSETEXP_H_ */
