#ifndef FREE_SEMIRING_H
#define FREE_SEMIRING_H

#include <string>
#include <memory>
#include "semiring.h"
#include "var.h"

class FreeSemiring : public Semiring<FreeSemiring>
{
	enum optype {Element, Addition, Multiplication, Star};
private:
	Var* elem;
	std::shared_ptr<FreeSemiring> left_ptr;
	std::shared_ptr<FreeSemiring> right_ptr;
	enum optype type;
public:
	FreeSemiring(Var var);
	FreeSemiring(const FreeSemiring& term);
	FreeSemiring(optype type, FreeSemiring left);
	FreeSemiring(optype type, FreeSemiring left, FreeSemiring right);
	virtual ~FreeSemiring();
	FreeSemiring operator + (const FreeSemiring& term) const;
	FreeSemiring operator * (const FreeSemiring& term) const;
	bool operator == (const FreeSemiring& term) const;
	FreeSemiring star () const;
	static FreeSemiring null;
	static FreeSemiring one;
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

#endif
