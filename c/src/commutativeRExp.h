/*
 * commutativeRExp.h
 *
 *  Created on: 30.08.2012
 *      Author: maxi
 */

#ifndef COMMUTATIVEREXP_H_
#define COMMUTATIVEREXP_H_
#include <set>
#include <string>
#include <memory>
#include "var.h"
#include "semiring.h"

/* commutative RegExp (commutative, associative, idempotent (CAI) addition and CA multiplication)
 storing them fully expanded might be a really bad idea (??) (huge blowup expected just because of applying the distributive law..?)
 thus we will store them non-expanded but we will (try) to use idempotence to reduce their size
 Examples:
	(a+(b.c+c.b).(a.b + c + b.a)*) = a + (b.c).(ab + c)*
	(a.b+c) + (c + b.a) = (a.b + c)
	(c + b.a) . (a.b + c) = (a.b+c).(a.b+c)
	(c + b.a) + (a.b + c) = (a.b + c)
	(a.b + a.c) + (a . (c+b)) = (a.b + a.c) +  (a . (b+c))

	so:
	- "sort" every expression (inductively) e.g. lexicographically
 	- apply idempotence, if sorted expressions are syntactically equal
 	- do not apply the distributive-law
 	- no semantic equivalence check which is hard (NP-hard, in fact even Pi_2-complete or something)

*/

/* TODO:
 * mögliche Rewrite-Regeln zur Optimierung :)

- x + (x)* = (x)*
[- 1 + x(x)* = x*] (enthalten im Nächsten wenn die +-Optimierung gemacht wird)
- 1 + (x)+ = x* (aber auch: 1 + (x)+ + (y)+ = (x)* + (y)*
- (x)* + (y)+ = (x)* + (y)*

*/



class CommutativeRExp : public Semiring<CommutativeRExp, Commutativity::Commutative, Idempotence::Idempotent>
{
public:
	enum optype {Empty, Element, Addition, Multiplication, Star, Plus};
	enum optype type;
private:
	VarId elem;
	std::shared_ptr<std::set<CommutativeRExp>> seta;
	std::shared_ptr<std::multiset<CommutativeRExp>> setm;
	std::shared_ptr<CommutativeRExp> rexp;
	void optimize_starplus();
	void optimize(bool one);
	std::string generateString() const;

public:
	static std::shared_ptr<CommutativeRExp> elem_null;
	static std::shared_ptr<CommutativeRExp> elem_one;

	CommutativeRExp();
	CommutativeRExp(int zero);
	CommutativeRExp(VarId var);
	CommutativeRExp(enum optype type, std::shared_ptr<std::set<CommutativeRExp>> seta);
	CommutativeRExp(enum optype type, std::shared_ptr<std::multiset<CommutativeRExp>> setm);
	CommutativeRExp(enum optype type, std::shared_ptr<CommutativeRExp> term);
	CommutativeRExp(const CommutativeRExp& expr);
	virtual ~CommutativeRExp();
	CommutativeRExp operator += (const CommutativeRExp& term);
	CommutativeRExp operator *= (const CommutativeRExp& term);
	bool operator < (const CommutativeRExp& term) const;
	bool operator == (const CommutativeRExp& term) const;
	CommutativeRExp star () const;
	static CommutativeRExp null();
	static CommutativeRExp one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

#endif /* COMMUTATIVEREXP_H_ */
