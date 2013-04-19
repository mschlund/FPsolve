#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "semirings/float-semiring.h"
#include "semirings/commutativeRExp.h"
#include "polynomials/polynomial.h"

#ifdef OLD_SEMILINEAR_SET
#include "semirings/semilinSetExp.h"
#else
#include "semirings/semilinear_set.h"
#endif

template <typename SR>
class Polynomial;

class FloatSemiring;
class CommutativeRExp;


class Parser
{
private:
public:
	Parser();
	std::vector<std::pair<VarId, Polynomial<FloatSemiring>>> float_parser(std::string input);
	std::vector<std::pair<VarId, Polynomial<CommutativeRExp>>> rexp_parser(std::string input);
	std::vector<std::pair<VarId, Polynomial<SemilinSetExp>>> slset_parser(std::string input);
        std::vector<std::pair<VarId, Polynomial<FreeSemiring>>> free_parser(std::string input);
};

#endif
