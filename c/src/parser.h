#ifndef PARSER_H
#define PARSER_H

#include <string>

#ifdef OLD_SEMILINEAR_SET
#include "semilinSetExp.h"
#else
#include "semilinear_set.h"
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
};

#endif
