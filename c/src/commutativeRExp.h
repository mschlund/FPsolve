/*
 * commutativeRExp.h
 *
 *  Created on: 30.08.2012
 *      Author: maxi
 */

#ifndef COMMUTATIVEREXP_H_
#define COMMUTATIVEREXP_H_
#include <initializer_list>
#include <set>
#include "counting-semiring.h"
#include "var.h"

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
public:
	typedef std::pair<std::multiset<Var>, std::set<Var> > linset;
private:
	std::set<linset> sets;
	static linset *elem_epsilon;
	CommutativeRExp(std::set<linset> sets);
	static std::set<linset> ExpUnion(std::set<linset> e1, std::set<linset> e2);
	static linset ExpConcat(linset e1, linset e2);
	static std::set<linset> ExpConcat(std::set<linset> e1, std::set<linset> e2);
	static std::set<linset> ExpStar(linset e);
	static linset epsilon();

public:
	static CommutativeRExp *elem_null;
	static CommutativeRExp *elem_one;

	CommutativeRExp();
	CommutativeRExp(Var var);
	CommutativeRExp(const CommutativeRExp& expr);
	virtual ~CommutativeRExp();
	CommutativeRExp operator + (const CommutativeRExp& term) const;
	CommutativeRExp operator * (const CommutativeRExp& term) const;
	bool operator == (const CommutativeRExp& term) const;
	CommutativeRExp star () const;
	static CommutativeRExp null();
	static CommutativeRExp one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

#endif /* COMMUTATIVEREXP_H_ */
