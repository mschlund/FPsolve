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

// returns a pointer to a map which changes the access direction
// you have to delete the map by yourself!
template <typename SR>
std::unordered_map<FreeSemiring, SR, FreeSemiring>* reverse_map(const std::unordered_map<SR, FreeSemiring, SR>& valuation)
{
	auto result = new std::unordered_map<FreeSemiring,SR,FreeSemiring>();
	for(auto v_it = valuation.begin(); v_it != valuation.end(); ++v_it)
	{
		result->insert(result->begin(), std::pair<FreeSemiring,SR>(v_it->second, v_it->first));
	}
	return result;
}

// add an FreeSemiring → SR pair to the map. Delete an existing mapping from free_elem to some old sr_elem
template <typename SR>
void add_valuation(FreeSemiring free_elem, SR sr_elem, std::unordered_map<FreeSemiring, SR, FreeSemiring>* valuation)
{
	valuation->erase(free_elem);
	valuation->insert(valuation->begin(), std::pair<FreeSemiring,SR>(free_elem,sr_elem));
}
#endif
