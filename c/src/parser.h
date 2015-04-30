#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "semirings/commutativeRExp.h"
#include "semirings/prefix-semiring.h"
#include "semirings/lossy-finite-automaton.h"

#ifdef USE_GENEPI
#include "semirings/semilinSetNdd.h"
#endif

#include "polynomials/commutative_polynomial.h"
#include "polynomials/non_commutative_polynomial.h"

template <typename SR>
class CommutativePolynomial;

class CommutativeRExp;

class Parser
{
private:
public:
  Parser();
  std::vector<std::pair<VarId, CommutativePolynomial<CommutativeRExp>>> rexp_parser(std::string input);
#ifdef USE_GENEPI
  std::vector<std::pair<VarId, CommutativePolynomial<SemilinSetNdd>>> slsetndd_parser(std::string input);
#endif
  std::vector<std::pair<VarId, NonCommutativePolynomial<FreeSemiring>>> free_parser(std::string input);
  std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> lossy_fa_parser(std::string input);
  std::vector<std::pair<VarId, NonCommutativePolynomial<PrefixSemiring>>> prefix_parser(std::string input, unsigned int length);
};

#endif
