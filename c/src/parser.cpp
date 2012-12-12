#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
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

struct free_var_impl
{
	template <typename T>
	struct result { typedef FreeSemiring type; }; // needed because this is how boost knows of what type the result will be

	//template <typename T>
	const FreeSemiring operator()(std::string& s) const // this function does the actual work
	{
		// create an element with the given var
		return FreeSemiring(Var::getVar(s));
	}
};
const phx::function<free_var_impl> free_var;

struct rexp_var_impl
{
	template <typename T>
	struct result { typedef CommutativeRExp type; }; // needed because this is how boost knows of what type the result will be

	//template <typename T>
	const CommutativeRExp operator()(std::string& s) const // this function does the actual work
	{
		// create an element with the given var
		return CommutativeRExp(Var::getVar(s));
	}
};
const phx::function<rexp_var_impl> rexp_var;

struct polyrexp_impl
{
	template <typename T,typename U>
	struct result { typedef Polynomial<CommutativeRExp> type; }; // needed because this is how boost knows of what type the result will be

	const Polynomial<CommutativeRExp> operator()(CommutativeRExp& s, std::vector<std::string>& v) const // this function does the actual work
	{
		// create an element with the given var
		std::vector<VarPtr> vars;
		for(auto it = v.begin(); it != v.end(); ++it)
			vars.push_back(Var::getVar(*it));
		return Polynomial<CommutativeRExp>({Monomial<CommutativeRExp>(s,vars)});
	}
};
const phx::function<polyrexp_impl> polyrexp;

struct polyrexp2_impl
{
	template <typename T>
	struct result { typedef Polynomial<CommutativeRExp> type; }; // needed because this is how boost knows of what type the result will be

	const Polynomial<CommutativeRExp> operator()(std::string& s) const // this function does the actual work
	{
		// create an element without a variable
		return Polynomial<CommutativeRExp>({Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar(s)),{})});
	}
};
const phx::function<polyrexp2_impl> polyrexp2;

struct polyvar_impl
{
	template <typename T>
	struct result { typedef Polynomial<CommutativeRExp> type; }; // needed because this is how boost knows of what type the result will be

	const Polynomial<CommutativeRExp> operator()(std::string& var) const // this function does the actual work
	{
		// create an element with the given var
		if(var == "1") // special case for the one-element
			return Polynomial<CommutativeRExp>::one();
		else
			return Polynomial<CommutativeRExp>({Monomial<CommutativeRExp>(CommutativeRExp::one(),{Var::getVar(var)})});
	}
};
const phx::function<polyvar_impl> polyvar;

// some boost::qi methods for pushing values through our parser
using qi::_val;
using qi::_1;
using qi::_2;
using qi::lit;
using qi::lexeme;
using qi::eps;

template <typename Iterator, typename Skipper>
struct float_parser : qi::grammar<Iterator, FloatSemiring(), Skipper>
{
	// specify the parser rules, term is the start rule (and in this case the only rule)
	float_parser() : float_parser::base_type(expression)
	{
		// what does a float expression look like?
		expression = term [_val  = _1] >> *( '+' >> term [_val = _val + _1] );
		term = starfactor [_val = _1] >> *( 'x' >> starfactor [_val = _val * _1] );
		starfactor = factor [_val = _1] >> -(lit('*'))[_val = star(_val)];
		factor =
			qi::float_ [_val = _1] |
			'(' >> expression [_val = _1] >> ')';
	}

	// this declare the available rules and the respecting return types of our parser
	qi::rule<Iterator, FloatSemiring(), Skipper> expression;
	qi::rule<Iterator, FloatSemiring(), Skipper> term;
	qi::rule<Iterator, FloatSemiring(), Skipper> starfactor;
	qi::rule<Iterator, FloatSemiring(), Skipper> factor;

};

template <typename Iterator, typename Skipper>
struct free_parser : qi::grammar<Iterator, FreeSemiring(), Skipper>
{
	// specify the parser rules, term is the start rule (and in this case the only rule)
	free_parser() : free_parser::base_type(expression)
	{
		// what does a free semiring term look like?
		expression = term [_val  = _1] >> *( '+' >> term [_val = _val + _1] );
		term = starfactor [_val = _1] >> *( '.' >> starfactor [_val = _val * _1] );
		starfactor = factor [_val = _1] >> -(lit('*'))[_val = star(_val)];
		factor =
			qi::as_string[lexeme[+qi::alpha]] [_val = free_var(_1)] |
			'(' >> expression [_val = _1] >> ')';
	}

	// this declare the available rules and the respecting return types of our parser
	qi::rule<Iterator, FreeSemiring(), Skipper> expression;
	qi::rule<Iterator, FreeSemiring(), Skipper> term;
	qi::rule<Iterator, FreeSemiring(), Skipper> starfactor;
	qi::rule<Iterator, FreeSemiring(), Skipper> factor;
};

template <typename Iterator, typename Skipper>
struct rexp_parser : qi::grammar<Iterator, CommutativeRExp(), Skipper>
{
	// specify the parser rules, term is the start rule (and in this case the only rule)
	rexp_parser() : rexp_parser::base_type(expression)
	{
		// what does a commutative regular expression term look like?
		expression = term [_val  = _1] >> *( '+' >> term [_val = _val + _1] );
		term = starfactor [_val = _1] >> *( '.' >> starfactor [_val = _val * _1] );
		starfactor = factor [_val = _1] >> -(lit('*'))[_val = star(_val)];
		factor =
			qi::as_string[lexeme[+qi::alpha]] [_val = rexp_var(_1)] |
			'(' >> expression [_val = _1] >> ')';
	}

	// this declare the available rules and the respecting return types of our parser
	qi::rule<Iterator, CommutativeRExp(), Skipper> expression;
	qi::rule<Iterator, CommutativeRExp(), Skipper> term;
	qi::rule<Iterator, CommutativeRExp(), Skipper> starfactor;
	qi::rule<Iterator, CommutativeRExp(), Skipper> factor;
};

template <typename Iterator, typename Skipper>
struct polyrexp_parser : qi::grammar<Iterator, Polynomial<CommutativeRExp>(), Skipper>
{
	// specify the parser rules, term is the start rule (and in this case the only rule)
	polyrexp_parser() : polyrexp_parser::base_type(expression)
	{
		// what does a polynomial in the commutative regular expression semiring look like?
		expression = term [_val  = _1] >> *( '+' >> term [_val = _val + _1] );
		term = ( rexp >> *( '.' >> var) )[_val = polyrexp(_1, _2)];
		var = +qi::as_string[lexeme[+qi::upper]] [_val = _1];
		//var = +qi::as_string[lexeme['<' >> +(ascii::char_ - '>') >> '>']] [_val = _1];
	}

	// this declare the available rules and the respecting return types of our parser
	qi::rule<Iterator, Polynomial<CommutativeRExp>(), Skipper> expression;
	qi::rule<Iterator, Polynomial<CommutativeRExp>(), Skipper> term;
	qi::rule<Iterator, std::string(), Skipper> var;
	rexp_parser<Iterator, Skipper> rexp;
};

//FIXME: empty lines in the input-file lead to a segfault :)
//TODO: allow rules going over multiple lines

template <typename Iterator, typename Skipper>
struct grammar_parser : qi::grammar<Iterator, std::pair<std::string, Polynomial<CommutativeRExp>>(), Skipper>
{
	// specify the parser rules, term is the start rule (and in this case the only rule)
	grammar_parser() : grammar_parser::base_type(rule)
	{
		// what does a polynomial in the commutative regular expression semiring look like?
		rule %= (varidentifier >> lexeme["::="] >> expression);
		expression = term [_val = _1] >> *('|' >> expression [_val = _val + _1]); // addition
		term = 	eps [_val = Polynomial<CommutativeRExp>::one()] >> // set _val to be the one-element of the semiring
			*(
				('(' >> expression >> ')') [_val = _val * _1] |
				literal	[_val = _val * _1] |
				var	[_val = _val * _1]  ); // multiplication
		//literal = '"' >> literalidentifier [_val = polyrexp2(_1)] >> '"';
		literal = sidentifier [_val = polyrexp2(_1)] |
				'"' >> literalidentifier [_val = polyrexp2(_1)] >> '"';
		var = varidentifier[_val = polyvar(_1)];
		literalidentifier = qi::as_string[lexeme[+(ascii::char_ - '"' - '|' -'('- ')' - ascii::space - '<' - '>')]];
		sidentifier = qi::as_string[lexeme[(ascii::char_ - '"' - '|' -'('- ')' - ascii::space - '<' - '>')]];
		//varidentifier = qi::as_string[lexeme[+(ascii::char_ - ':' - '=' - '|' - '(' - ')' - ascii::space)]];
		varidentifier = qi::as_string[lexeme['<' >> +(ascii::char_ - '>') >> '>']];
	}

	// this declares the available rules and the respecting return types of our parser
	qi::rule<Iterator, std::pair<std::string, Polynomial<CommutativeRExp>>(), Skipper> rule;
	qi::rule<Iterator, Polynomial<CommutativeRExp>(), Skipper> expression;
	qi::rule<Iterator, Polynomial<CommutativeRExp>(), Skipper> term;
	qi::rule<Iterator, Polynomial<CommutativeRExp>(), Skipper> var;
	qi::rule<Iterator, Polynomial<CommutativeRExp>(), Skipper> literal;
	qi::rule<Iterator, std::string(), Skipper> varidentifier;
	qi::rule<Iterator, std::string(), Skipper> literalidentifier;
	qi::rule<Iterator, std::string(), Skipper> sidentifier;
};

typedef std::string::const_iterator iterator_type;

FloatSemiring Parser::parse_float(std::string input)
{
	typedef float_parser<iterator_type, qi::space_type> float_parser;
	float_parser floater;

	iterator_type iter = input.begin();
	iterator_type end = input.end();

	FloatSemiring result;
	if(!(qi::phrase_parse(iter, end, floater, qi::space, result) && iter == end))
		std::cout << "bad input, failed at: " << std::string(iter, end) << std::endl;

	return result;
}


FreeSemiring Parser::parse_free(std::string input)
{
	typedef free_parser<iterator_type, qi::space_type> free_parser;
	free_parser freeer;

	iterator_type iter = input.begin();
	iterator_type end = input.end();

	FreeSemiring result;
	if(!(qi::phrase_parse(iter, end, freeer, qi::space, result) && iter == end))
		std::cout << "bad input, failed at: " << std::string(iter, end) << std::endl;

	return result;
}

CommutativeRExp Parser::parse_rexp(std::string input)
{
	typedef rexp_parser<iterator_type, qi::space_type> rexp_parser;
	rexp_parser rexper;

	iterator_type iter = input.begin();
	iterator_type end = input.end();

	CommutativeRExp result;
	if(!(qi::phrase_parse(iter, end, rexper, qi::space, result) && iter == end))
		std::cout << "bad input, failed at: " << std::string(iter, end) << std::endl;

	return result;
}

Polynomial<CommutativeRExp> Parser::parse_polyrexp(std::string input)
{
	typedef polyrexp_parser<iterator_type, qi::space_type> polyrexp_parser;
	polyrexp_parser polyrexper;

	iterator_type iter = input.begin();
	iterator_type end = input.end();

	Polynomial<CommutativeRExp> result;
	if(!(qi::phrase_parse(iter, end, polyrexper, qi::space, result) && iter == end))
		std::cout << "bad input, failed at: " << std::string(iter, end) << std::endl;


	return result;
}

std::pair<VarPtr, Polynomial<CommutativeRExp>> Parser::parse_grammar(std::string input)
{
	typedef grammar_parser<iterator_type, qi::space_type> grammar_parser;
	grammar_parser grammarer;

	iterator_type iter = input.begin();
	iterator_type end = input.end();

	std::pair<std::string, Polynomial<CommutativeRExp>> result;
	if(!(qi::phrase_parse(iter, end, grammarer, qi::space, result) && iter == end))
		std::cout << "bad input, failed at: " << std::string(iter, end) << std::endl;

	return std::pair<VarPtr, Polynomial<CommutativeRExp>>(Var::getVar(result.first), result.second);
}
