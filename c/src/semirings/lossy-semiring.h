#ifndef LOSSY_SEMIRING_H_
#define LOSSY_SEMIRING_H_

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "semiring.h"

template<typename SR>
class Evaluator;

class LossySemiring: public StarableSemiring<LossySemiring,
		Commutativity::NonCommutative, Idempotence::Idempotent> {
public:

	static LossySemiring null();
	static LossySemiring one();

	LossySemiring star() const;

	LossySemiring operator+(const LossySemiring &x);
	LossySemiring& operator+=(const LossySemiring &x);
	LossySemiring operator*(const LossySemiring &x);
	LossySemiring& operator*=(const LossySemiring &x);

	bool operator==(const LossySemiring &x) const;

	std::string string() const;

protected:

	friend struct std::hash<LossySemiring>;

};

#endif
