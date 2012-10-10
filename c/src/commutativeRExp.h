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
#include "counting-semiring.h"
#include "var.h"

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


class CommutativeRExp : public Semiring<CommutativeRExp>
{
public:
	enum optype {Empty, Element, Addition, Multiplication, Star};
	enum optype type;
private:
	VarPtr elem;
	std::set<CommutativeRExp>* seta;
	std::multiset<CommutativeRExp>* setm;
	CommutativeRExp* rexp;
	std::string str;
	std::string generateString() const;

public:
	static CommutativeRExp *elem_null;
	static CommutativeRExp *elem_one;

	CommutativeRExp();
	CommutativeRExp(VarPtr var);
	CommutativeRExp(enum optype type, std::set<CommutativeRExp>* seta);
	CommutativeRExp(enum optype type, std::multiset<CommutativeRExp>* setm);
	CommutativeRExp(enum optype type, const CommutativeRExp& term);
	CommutativeRExp(const CommutativeRExp& expr);
	virtual ~CommutativeRExp();
	CommutativeRExp operator + (const CommutativeRExp& term) const;
	CommutativeRExp operator * (const CommutativeRExp& term) const;
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
