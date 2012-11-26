#ifndef PARSER_H
#define PARSER_H

#include <string>
#include "float-semiring.h"
#include "free-semiring.h"
#include "commutativeRExp.h"
#include "polynomial.h"

class Parser
{
private:
public:
	Parser();
	FloatSemiring parse_float(std::string input);
	FreeSemiring parse_free(std::string input);
	CommutativeRExp parse_rexp(std::string input);
	Polynomial<CommutativeRExp> parse_polyrexp(std::string input);
};

#endif
