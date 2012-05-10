#ifndef FLOAT_SEMIRING_H
#define FLOAT_SEMIRING_H

#include <string>
#include "semiring.h"

class FloatSemiring : public Semiring<FloatSemiring>
{
private:
	float val;
public:
	FloatSemiring(const float val);
	~FloatSemiring();
	FloatSemiring operator + (const FloatSemiring& elem);
	FloatSemiring operator * (const FloatSemiring& elem);
	FloatSemiring star ();
	static FloatSemiring null();
	std::string string();
	bool is_idempotent();
	bool is_commutative();
};

#endif
