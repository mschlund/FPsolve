#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "semirings/commutativeRExp.h"
#include "semirings/prefix-semiring.h"
#include "semirings/lossy-finite-automaton.h"

#include "datastructs/equations.h"

#ifdef USE_GENEPI
#include "semirings/semilinSetNdd.h"
#endif

#include "polynomials/commutative_polynomial.h"
#include "polynomials/non_commutative_polynomial.h"
#include "polynomials/lossy_non_commutative_polynomial.h"

template <typename SR>
class CommutativePolynomial;

class CommutativeRExp;

class Parser
{
private:
public:
  Parser();
  Equations<CommutativeRExp> rexp_parser(std::string input);
#ifdef USE_GENEPI
  Equations<SemilinSetNdd> slsetndd_parser(std::string input);
#endif

  NCEquations<FreeSemiring> free_parser(std::string input);
  NCEquations<LossyFiniteAutomaton> lossy_fa_parser(std::string input);
  NCEquations<PrefixSemiring> prefix_parser(std::string input, unsigned int length);
};

#endif
