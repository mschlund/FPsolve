#ifndef BOOL_SEMIRING_H
#define BOOL_SEMIRING_H

#include <string>
#include <memory>

#include "semiring.h"

class BoolSemiring : public Semiring<BoolSemiring, Commutativity::Commutative, Idempotence::Idempotent>
{
private:
	bool val;
	static std::shared_ptr<BoolSemiring> elem_null;
	static std::shared_ptr<BoolSemiring> elem_one;
public:
	BoolSemiring();
	BoolSemiring(bool val);
	virtual ~BoolSemiring();
	BoolSemiring operator += (const BoolSemiring& elem);
	BoolSemiring operator *= (const BoolSemiring& elem);
	bool operator == (const BoolSemiring& elem) const;
	BoolSemiring star () const;
	static BoolSemiring null();
	static BoolSemiring one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

#endif
