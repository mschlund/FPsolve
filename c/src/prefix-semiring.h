#ifndef PREFIX_SEMIRING_H
#define PREFIX_SEMIRING_H

#include <string>
#include <vector>
#include <set>
#include "semiring.h"
#include "var.h"

class PrefixSemiring : public Semiring<PrefixSemiring>
{
private:
	std::set<std::vector<VarPtr>> val;
	static unsigned int max_length;
	static std::vector<VarPtr> concatenate(std::vector<VarPtr> l, std::vector<VarPtr> r);
	static std::shared_ptr<PrefixSemiring> elem_null;
	static std::shared_ptr<PrefixSemiring> elem_one;
public:
	PrefixSemiring();
	PrefixSemiring(const std::vector<VarPtr>& val);
	virtual ~PrefixSemiring();
	PrefixSemiring operator += (const PrefixSemiring& elem);
	PrefixSemiring operator *= (const PrefixSemiring& elem);
	bool operator == (const PrefixSemiring& elem) const;
	PrefixSemiring star () const;
	static PrefixSemiring null();
	static PrefixSemiring one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

#endif
