#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include "parser.h"

Parser::Parser()
{
}

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phx = boost::phoenix;

// this function object is used for star'ing elements from within the parser
// the templates achieve, that this can be used for every type (hopefully a semiring)
// TODO: enforce that T is subclass of Semiring!
struct star_impl
{
	template <typename T> struct result { typedef T type; }; // needed because this is how boost knows of what type the result will be

	template <typename T>
	T operator()(T& s) const // this function does the actual work
	{
		return s.star();
	}
};
const phx::function<star_impl> star;

struct rexp_var_impl
{
	template <typename T>
	struct result { typedef CommutativeRExp type; }; // this tells Boost the return type

	const CommutativeRExp operator()(std::string& s) const
	{
		// create an element with the given var
		return CommutativeRExp(Var::getVar(s));
	}
};
const phx::function<rexp_var_impl> rexp_var;

struct var_impl
{
	template <typename T>
	struct result { typedef VarPtr type; }; // this tells Boost the return type

	const VarPtr operator()(std::string& var) const
	{
		return Var::getVar(var);
	}
};
const phx::function<var_impl> variable;

// some boost::qi methods for pushing values through our parser
using qi::_val;
using qi::_1;
using qi::_2;
using qi::lit;
using qi::lexeme;
using qi::eps;
using phx::ref; // must be used as phx::ref, otherwise g++ is confused with std::ref...
using phx::insert;
using phx::end;
using phx::begin;

typedef std::string::const_iterator iterator_type;

// parser for a float semiring element
struct float_elem_parser : qi::grammar<iterator_type, FloatSemiring()>
{
	float_elem_parser() : float_elem_parser::base_type(elem)
	{
		elem %= qi::float_;
	}
	qi::rule<iterator_type, FloatSemiring()> elem;
};

// parser for a commutative regular expression semiring element
struct rexp_elem_parser : qi::grammar<iterator_type, CommutativeRExp()>
{
	rexp_elem_parser() : rexp_elem_parser::base_type(elem)
	{
		elem = '"' >> qi::as_string[lexeme[+(ascii::char_ -'"')]] [_val = rexp_var(_1)] >> '"';
	}
	qi::rule<iterator_type, CommutativeRExp()> elem;
};


// use the given parser and return a list of equations of SR polynomials
template <typename SR_Parser, typename SR>
struct equation_parser : qi::grammar<iterator_type, std::vector<std::pair<VarPtr, Polynomial<SR>>>(), qi::space_type>
{
	std::vector<std::pair<VarPtr, Polynomial<SR>>> new_rules; // accumulate generated rules in this vector
	equation_parser() : equation_parser::base_type(equations)
	{
		// the first rule combines the generated equations with the read equations and returns them in one vector
		equations = (*equation)[_val = _1] [insert(_val, end(_val), begin(phx::ref(new_rules)), end(phx::ref(new_rules)))];
		equation %= (var >> lexeme["::="] >> polynomial >> ';');
		polynomial = summand [_val = _1] >> *('|' >> summand [_val = _val + _1]); // addition
		summand = eps [_val = Polynomial<SR>::one()] >> // set _val to one-element of the semiring
				*(
					('(' >> polynomial >> ')') [_val = _val * _1] | // parenthesised term
					var [_val = _val * _1] | // variable factor
					sr_elem [_val = _val * _1] // sr-elem factor
				 );
		var = varidentifier[_val = variable(_1)];
		varidentifier = qi::as_string[lexeme['<' >> +(ascii::char_ - '>') >> '>']];
	}
	qi::rule<iterator_type, std::vector<std::pair<VarPtr, Polynomial<SR>>>(), qi::space_type> equations;
	qi::rule<iterator_type, std::pair<VarPtr, Polynomial<SR>>(), qi::space_type> equation;
	qi::rule<iterator_type, Polynomial<SR>(), qi::space_type> polynomial;
	qi::rule<iterator_type, Polynomial<SR>(), qi::space_type> summand;
	qi::rule<iterator_type, std::string()> varidentifier;
	qi::rule<iterator_type, VarPtr()> var;
	SR_Parser sr_elem;
};

// generic function for using the parser
template <typename SR_Parser, typename SR> 
std::vector<std::pair<VarPtr, Polynomial<SR>>> parser(std::string input)
{
	typedef equation_parser<SR_Parser, SR> equation_parser;
	equation_parser equation;

	iterator_type iter = input.begin();
	iterator_type end = input.end();

	std::vector<std::pair<VarPtr, Polynomial<SR>>> result;
	if(!(qi::phrase_parse(iter, end, equation, qi::space, result) && iter == end))
		std::cout << "bad input, failed at: " << std::string(iter, end) << std::endl;

	return result;
}

// wrapper function for float equations
std::vector<std::pair<VarPtr, Polynomial<FloatSemiring>>> Parser::float_parser(std::string input)
{
	return parser<float_elem_parser, FloatSemiring>(input);
}

// wrapper function for regular expression equations
std::vector<std::pair<VarPtr, Polynomial<CommutativeRExp>>> Parser::rexp_parser(std::string input)
{
	return parser<rexp_elem_parser, CommutativeRExp>(input);
}
