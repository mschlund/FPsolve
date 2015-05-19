#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "semirings/commutativeRExp.h"
#include "semirings/prefix-semiring.h"

#ifdef USE_LIBFA
#include "semirings/lossy-finite-automaton.h"
#include "polynomials/lossy_non_commutative_polynomial.h"
#endif

#include "datastructs/equations.h"

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
  Equations<CommutativeRExp> rexp_parser(std::string input);
#ifdef USE_GENEPI
  Equations<SemilinSetNdd> slsetndd_parser(std::string input);
#endif

  NCEquations<FreeSemiring> free_parser(std::string input);
  NCEquations<PrefixSemiring> prefix_parser(std::string input, unsigned int length);

#ifdef USE_LIBFA
  NCEquations<LossyFiniteAutomaton> lossy_fa_parser(std::string input);
#endif


};

#endif
