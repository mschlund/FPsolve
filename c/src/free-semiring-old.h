#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "matrix.h"
#include "semiring.h"
#include "var.h"


class FreeSemiring : public Semiring<FreeSemiring>
{
public:
	VarPtr elem;
	std::shared_ptr<FreeSemiring> left_ptr;
	std::shared_ptr<FreeSemiring> right_ptr;
	enum optype {Element, Addition, Multiplication, Star, Dummy};
	enum optype type;
	static std::shared_ptr<FreeSemiring> elem_null;
	static std::shared_ptr<FreeSemiring> elem_one;

	FreeSemiring();
	FreeSemiring(int zero);
	FreeSemiring(VarPtr var);
	FreeSemiring(const FreeSemiring& term);
	FreeSemiring(optype type, FreeSemiring left);
	FreeSemiring(optype type, FreeSemiring left, FreeSemiring right);
	virtual ~FreeSemiring();
	FreeSemiring operator += (const FreeSemiring& term);
	FreeSemiring operator *= (const FreeSemiring& term);
	bool operator == (const FreeSemiring& term) const;
	FreeSemiring star () const;
	static FreeSemiring null();
	static FreeSemiring one();
	std::string string() const;
	static bool is_idempotent;
	static bool is_commutative;
};

template <typename SR>
SR FreeSemiring_eval(FreeSemiring elem, std::unordered_map<VarPtr, SR>* valuation)
{
	SR result;
	switch(elem.type)
	{
	case FreeSemiring::Element:
	{
          if (Var::getVar("Null") == elem.elem) {
            return SR::null();
          } else if (Var::getVar("One") == elem.elem) {
            return SR::one();
          } else {
		typename std::unordered_map<VarPtr,SR>::const_iterator tmp = valuation->find(elem.elem);
		assert(tmp!=valuation->end());
		result = tmp->second;
          }
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
	case FreeSemiring::Dummy:
		assert(false);
	};
	return result;
}

template <typename SR>
Matrix<SR> FreeSemiring_eval(Matrix<FreeSemiring> matrix, std::unordered_map<VarPtr, SR>* valuation)
{
	std::vector<FreeSemiring> elements = matrix.getElements();
	std::vector<SR> ret;

	for(unsigned int i=0; i<elements.size(); i++)
	{
		ret.push_back(FreeSemiring_eval<SR>(elements.at(i),valuation));
	}

	return Matrix<SR>(matrix.getRows(), matrix.getColumns(), ret);
}
