/*
 * commutativeRExp.h
 *
 *  Created on: 30.08.2012
 *      Author: maxi
 */

#ifndef COMMUTATIVEREXP_H_
#define COMMUTATIVEREXP_H_

// Stores commutative regular expressions fully expanded as a list (or set?) (->finite union) of linear sets
// each linear set is just a pair of a set (offset) and a set (->idempotent sum) of multisets (-> commutative products) of variables (generators)

// Ok... storing them fully expanded might be a really bad idea (??) (huge blowup expected just because of applying the distributive law..?)
// we will worry about that later :)

// use the identities
// (1) (x*)* = x*
// (2) (x+y)* = x*y*
// (3) (xy*)* = 1 + xx*y*   to push stars inwards
// use distributive law when multiplying


class CommutativeRExp : public Semiring<CommutativeRExp>
{

};

#endif /* COMMUTATIVEREXP_H_ */
