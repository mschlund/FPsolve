#ifndef FREE_SEMIRING_H
#define FREE_SEMIRING_H

#include <string>
#include <memory>
#include <unordered_map>
#include "semiring.h"
#include "var.h"
#include "matrix.h"

class FreeSemiring : public Semiring<FreeSemiring>
{
public:
	Var* elem;
	std::shared_ptr<FreeSemiring> left_ptr;
	std::shared_ptr<FreeSemiring> right_ptr;
	enum optype {Element, Addition, Multiplication, Star};
	enum optype type;
	static FreeSemiring* elem_null;
	static FreeSemiring* elem_one;

	FreeSemiring();
	FreeSemiring(Var var);
	FreeSemiring(const FreeSemiring& term);
	FreeSemiring(optype type, FreeSemiring left);
	FreeSemiring(optype type, FreeSemiring left, FreeSemiring right);
	virtual ~FreeSemiring();
	FreeSemiring operator + (const FreeSemiring& term) const;
	FreeSemiring operator * (const FreeSemiring& term) const;
	bool operator == (const FreeSemiring& term) const;
	FreeSemiring star () const;
	static FreeSemiring null();
	static FreeSemiring one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

template <typename SR>
SR FreeSemiring_eval(FreeSemiring elem, std::unordered_map<FreeSemiring, SR, FreeSemiring>* valuation)
{
	SR result;
	switch(elem.type)
	{
	case FreeSemiring::Element:
	{
		typename std::unordered_map<FreeSemiring,SR,FreeSemiring>::const_iterator tmp = valuation->find(elem);
		assert(tmp!=valuation->end());
		result = tmp->second;
	}
		break;
	case FreeSemiring::Addition:
		result = FreeSemiring_eval(*elem.left_ptr, valuation) + FreeSemiring_eval(*elem.right_ptr, valuation);
		break;
	case FreeSemiring::Multiplication:
		result = FreeSemiring_eval(*elem.left_ptr, valuation) * FreeSemiring_eval(*elem.right_ptr, valuation);
		break;
	case FreeSemiring::Star:
		result = FreeSemiring_eval(*elem.left_ptr, valuation).star();
		break;
	};
	return result;
};

template <typename SR>
Matrix<SR> FreeSemiring_eval(Matrix<FreeSemiring> matrix, std::unordered_map<FreeSemiring, SR, FreeSemiring>* valuation)
{
	std::vector<FreeSemiring> elements = matrix.getElements();
	std::vector<SR> ret;

	for(unsigned int i=0; i<elements.size(); i++)
	{
		ret.push_back(FreeSemiring_eval<SR>(elements.at(i),valuation));
	}

	return Matrix<SR>(matrix.getRows(), matrix.getColumns(), ret);
};

#endif
