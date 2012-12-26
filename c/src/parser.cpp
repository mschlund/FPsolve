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

// this function object is used for star'ing from within the parser
// a new equation with the given new variable is created and returned
struct star_impl
{
	template <typename T, typename SR> struct result { typedef std::pair<VarPtr, SR> type; }; // needed because this is how boost knows of what type the result will be

	template <typename SR>
	std::pair<VarPtr, SR> operator()(VarPtr var, SR& s) const // var → 1 + var×s
	{
		return std::pair<VarPtr, SR>(var, (SR::one() + SR(var)*s));
	}
};
const phx::function<star_impl> star_equation;

struct option_impl
{
	template <typename T> struct result { typedef T type; }; // needed because this is how boost knows of what type the result will be

	template <typename SR>
	SR operator()(SR& s) const // A → 1 + A
	{
		return (SR::one() + s);
	}
};
const phx::function<option_impl> option;

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

struct new_var_impl
{
	template <typename T>
	struct result { typedef VarPtr type; }; // this tells Boost the return type

	template <typename T>
	const VarPtr operator()(T t) const // TODO: why has there to be an argument??? compiling fails without dummy argument... overloading???
	{
		return Var::getVar(); // return a fresh anonymous variable
	}
};
const phx::function<new_var_impl> new_var;

// some boost::qi methods for pushing values through our parser
using qi::_val;
using qi::_1;
using qi::_2;
using qi::_a;
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
					('(' >> polynomial >> ')') [_val = _val * _1] | // group term
					('{' >> kstar >> '}') [_val = _val * _1] | // | // kleene star: {A} → B && B = 1 + BA
					('[' >> polynomial >>']')  [_val = _val * option(_1)] | // optional term: [A] → 1 + A
					var [_val = _val * _1] | // variable factor
					sr_elem [_val = _val * _1] // sr-elem factor
				 );
		var = varidentifier[_val = variable(_1)];
		kstar = polynomial	[_a = new_var(0)] // create a new variable _a for the new equation
					[phx::push_back(phx::ref(new_rules),star_equation(_a,_1))] // create a new equation _a = <stared equation> and save it
					[_val = _a]; // continue with the new variable
		varidentifier = qi::as_string[lexeme['<' >> +(ascii::char_ - '>') >> '>']];
	}
	qi::rule<iterator_type, std::vector<std::pair<VarPtr, Polynomial<SR>>>(), qi::space_type> equations;
	qi::rule<iterator_type, std::pair<VarPtr, Polynomial<SR>>(), qi::space_type> equation;
	qi::rule<iterator_type, Polynomial<SR>(), qi::space_type> polynomial;
	qi::rule<iterator_type, Polynomial<SR>(), qi::space_type> summand;
	qi::rule<iterator_type, std::string()> varidentifier;
	qi::rule<iterator_type, VarPtr(), qi::locals<VarPtr>, qi::space_type> kstar;
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
