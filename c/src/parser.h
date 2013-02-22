#ifndef PARSER_H
#define PARSER_H

#include <string>
#include "float-semiring.h"
#include "commutativeRExp.h"
#include "polynomial.h"

#ifdef OLD_SEMILINEAR_SET
#include "semilinSetExp.h"
#else
#include "semilinear_set.h"
#endif


class Parser
{
private:
public:
	Parser();
	std::vector<std::pair<VarPtr, Polynomial<FloatSemiring>>> float_parser(std::string input);
	std::vector<std::pair<VarPtr, Polynomial<CommutativeRExp>>> rexp_parser(std::string input);
	std::vector<std::pair<VarPtr, Polynomial<SemilinSetExp>>> slset_parser(std::string input);
};

#endif
