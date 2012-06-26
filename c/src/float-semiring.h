#ifndef FLOAT_SEMIRING_H
#define FLOAT_SEMIRING_H

#include <string>
#include "semiring.h"

class FloatSemiring : public Semiring<FloatSemiring>
{
private:
	float val;
public:
	FloatSemiring();
	FloatSemiring(const float val);
	~FloatSemiring();
	FloatSemiring operator + (const FloatSemiring& elem) const;
	FloatSemiring operator * (const FloatSemiring& elem) const;
	bool operator == (const FloatSemiring elem) const;
	FloatSemiring star () const;
	static FloatSemiring null();
	static FloatSemiring one();
	std::string string() const;
	bool is_idempotent() const ;
	bool is_commutative() const;
};

#endif
