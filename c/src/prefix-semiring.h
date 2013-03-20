#ifndef PREFIX_SEMIRING_H
#define PREFIX_SEMIRING_H

#include <string>
#include <vector>
#include <set>
#include "semiring.h"
#include "var.h"

class PrefixSemiring : public Semiring<PrefixSemiring, Commutativity::NonCommutative, Idempotence::Idempotent>
{
private:
	std::set<std::vector<VarId>> val;
	static unsigned int max_length;
	static std::vector<VarId> concatenate(std::vector<VarId> l, std::vector<VarId> r);
	static std::shared_ptr<PrefixSemiring> elem_null;
	static std::shared_ptr<PrefixSemiring> elem_one;
public:
	PrefixSemiring();
	PrefixSemiring(const std::vector<VarId>& val);
	virtual ~PrefixSemiring();
	PrefixSemiring operator += (const PrefixSemiring& elem);
	PrefixSemiring operator *= (const PrefixSemiring& elem);
	bool operator == (const PrefixSemiring& elem) const;
	PrefixSemiring star () const;
	static PrefixSemiring null();
	static PrefixSemiring one();
	std::string string() const;
};

#endif
