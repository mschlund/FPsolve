#ifndef FREE_SEMIRING_H
#define FREE_SEMIRING_H

#include <memory>
#include <string>
#include <unordered_map>

#include "matrix.h"
#include "semiring.h"
#include "var.h"

#include "free-structure.h"

template <typename SR>
class Evaluator;

class FreeSemiring2 : public Semiring<FreeSemiring2> {
  public:
    FreeSemiring2(const VarPtr var) {
      assert(var);
      InitFactory();
      assert(factory_);
      node_ = factory_->NewElement(var);
      assert(node_);
    }

    static FreeSemiring2 null() {
      InitFactory();
      assert(factory_);
      return FreeSemiring2{factory_->GetEmpty()};
    }

    static FreeSemiring2 one() {
      InitFactory();
      assert(factory_);
      return FreeSemiring2{factory_->GetEpsilon()};
    }

    FreeSemiring2 star() const {
      assert(factory_);
      return FreeSemiring2{factory_->NewStar(node_)};
    }

    FreeSemiring2 operator+(const FreeSemiring2 &x) {
      assert(factory_);
      return FreeSemiring2{factory_->NewAddition(node_, x.node_)};
    }

    FreeSemiring2 operator*(const FreeSemiring2 &x) {
      assert(factory_);
      return FreeSemiring2{factory_->NewMultiplication(node_, x.node_)};
    }

    bool operator==(const FreeSemiring2 &x) const {
      return node_ == x.node_;
    }

    std::string string() const {
      std::stringstream ss;
      ss << *node_;
      return ss.str();
    }

    template <typename SR>
    SR Eval(const std::unordered_map<VarPtr, SR> &valuation) const;

    template <typename SR>
    SR Eval(Evaluator<SR> &evaluator) const;

    void PrintDot(std::ostream &out) {
      assert(factory_);
      factory_->PrintDot(out);
    }

    void GC() {
      assert(factory_);
      factory_->GC();
    }

  private:
    FreeSemiring2(NodePtr n) : node_(n) { assert(factory_); }

    static void InitFactory() {
      if (factory_ == nullptr) {
        factory_ = std::move(std::unique_ptr<NodeFactory>(new NodeFactory));
      }
    }

    NodePtr node_;
    static std::unique_ptr<NodeFactory> factory_;

    friend struct std::hash<FreeSemiring2>;
};

namespace std {

template <>
struct hash<FreeSemiring2> {
  inline std::size_t operator()(const FreeSemiring2 &fs) const {
    std::hash<NodePtr> h;
    return h(fs.node_);
  }
};

}  /* namespace std */


/*
 * Evaluator
 *
 * We want to memoize the result of evaluating every subgraph, since we can
 * reach the same node (and thus the same subgraph) many times.  So only the
 * first one we will actually perform the computation, in all subsequent visits
 * we will reuse the memoized reslt.  Note that this is ok, because we never
 * modify the actual semiring values.  We use shared_ptrs to avoid any
 * unnecessary copying of temporary values, etc.
 */
template <typename SR>
class Evaluator : public NodeVisitor {
  typedef std::shared_ptr<SR> SRPtr_;
  public:
    Evaluator(const std::unordered_map<VarPtr, SR> &v)
        : val_(v), evaled_(), result_() {}

    void Visit(const Addition &a) {
      LookupEval(a.GetLhs());
      auto temp = std::move(result_);
      LookupEval(a.GetRhs());
      result_ = std::make_shared<SR>(*temp + *result_);
    }

    void Visit(const Multiplication &m) {
      LookupEval(m.GetLhs());
      auto temp = std::move(result_);
      LookupEval(m.GetRhs());
      result_ = std::make_shared<SR>(*temp * *result_);
    }

    void Visit(const Star &s) {
      LookupEval(s.GetNode());
      result_ = std::make_shared<SR>(result_->star());
    }

    void Visit(const Element &e) {
      auto iter = val_.find(e.GetVar());
      assert(iter != val_.end());
      result_ = std::make_shared<SR>(iter->second);
    }

    void Visit(const Epsilon &e) {
      result_ = std::make_shared<SR>(SR::one());
    }

    void Visit(const Empty &e) {
      result_ = std::make_shared<SR>(SR::null());
    }

    std::shared_ptr<SR> GetResult() {
      return result_;
    }

  private:
    void LookupEval(const NodePtr &node) {
      auto iter = evaled_.find(node);
      if (iter != evaled_.end()) {
        result_ = iter->second;
      } else {
        node->Accept(*this);  /* Sets the result_ correctly */
        evaled_.emplace(node, result_);
      }
    }

    const std::unordered_map<VarPtr, SR> &val_;
    std::unordered_map<NodePtr, SRPtr_> evaled_;
    SRPtr_ result_;
};

template <typename SR>
SR FreeSemiring2::Eval(const std::unordered_map<VarPtr, SR> &valuation) const {
  Evaluator<SR> evaluator{valuation};
  node_->Accept(evaluator);
  return std::move(*evaluator.GetResult());
}

template <typename SR>
SR FreeSemiring2::Eval(Evaluator<SR> &evaluator) const {
  node_->Accept(evaluator);
  return std::move(*evaluator.GetResult());
}

template <typename SR>
Matrix<SR> FreeSemiringMatrixEval(const Matrix<FreeSemiring2> &matrix,
    const std::unordered_map<VarPtr, SR> &valuation) {

  std::vector<FreeSemiring2> elements = matrix.getElements();
  std::vector<SR> result;

  /* We have a single Evaluator that is used for all evaluations of the elements
   * of the original matrix.  This way if different elements refer to the same
   * FreeSemiring subexpression, we memoize the result and reuse it. */
  Evaluator<SR> evaluator{valuation};

  for(auto &elem : elements) {
    result.emplace_back(elem.Eval(evaluator));
  }

  return Matrix<SR>(matrix.getRows(), matrix.getColumns(), std::move(result));
}


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
	case FreeSemiring::Dummy:
		assert(false);
	};
	return result;
}

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
}

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
