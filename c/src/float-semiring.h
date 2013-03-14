#ifndef FLOAT_SEMIRING_H
#define FLOAT_SEMIRING_H

#include <string>
#include <memory>

#include "semiring.h"

class FloatSemiring : public Semiring<FloatSemiring>
{
private:
	float val;
	static std::shared_ptr<FloatSemiring> elem_null;
	static std::shared_ptr<FloatSemiring> elem_one;
public:
	FloatSemiring();
	FloatSemiring(const float val);
	virtual ~FloatSemiring();
	FloatSemiring operator += (const FloatSemiring& elem);
	FloatSemiring operator *= (const FloatSemiring& elem);
	bool operator == (const FloatSemiring& elem) const;
	FloatSemiring star () const;
	static FloatSemiring null();
	static FloatSemiring one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

#endif
