#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <map>
#include <set>
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include <assert.h>
#include "semiring.h"

template <typename SR>
class Polynomial : public Semiring<Polynomial<SR> >
{
private:
	int degree;
	std::set<char> variables;
	typedef std::map<std::string, SR> Tcoeff;
	Tcoeff coeff;

	// non-commutative version
	std::list<std::string> monomialDerivative(const char& var, const std::string& vars) const
	{
		std::list<std::string> ret;
		// (uvwx) = u'vwx + uv'wx + uvw'x + uvwx'
		// but uv'wx = 0 for v != var because uv'wx = u0wx = 0
		// so just use all terms uv'wx for v = var -> uv'wx = uwx
		for (int pos = 0; pos < vars.length(); ++pos)
		{
			if(vars[pos] == var)
			{
				std::string tmp = vars;
				tmp.erase(pos,1); // uv'wx -> uwx
				ret.push_back(tmp);
			}
		}
		return ret;
	};

	// insert a monomial into coeffs
	static void insertMonomial(const std::string& vars, const SR& coeff, Tcoeff* coeffs)
	{
		typename Tcoeff::const_iterator elem = coeffs->find(vars);
		if(elem == coeffs->end()) // not in the map yet
			(*coeffs)[vars] = coeff;
		else // there is already a coefficient
			(*coeffs)[vars] = (*coeffs)[vars] + coeff;
	};
public:

	// empty polynomial
	Polynomial()
	{
	};

	Polynomial(const Polynomial& polynomial)
	{
		this->variables = polynomial.variables;
		this->coeff = polynomial.coeff;
		this->degree = polynomial.degree;
	}

	Polynomial& operator=(const Polynomial& polynomial)
	{
		this->variables = polynomial.variables;
		this->coeff = polynomial.coeff;
		this->degree = polynomial.degree;
		return (*this);
	}

	Polynomial(std::set<char> variables, Tcoeff coeff)
	{
		this->variables = variables;
		this->coeff = coeff;
		// calculate the degree. will not change because polynomial is read only
		this->degree = 0;
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			this->degree = it->first.length() > this->degree ? it->first.length() : this->degree;
		}

	};

	Polynomial<SR> operator+(const Polynomial<SR>& poly) const
	{
		Tcoeff ret = this->coeff; // fill coefficients with all elements of the first operand
		for (typename Tcoeff::const_iterator it = poly.coeff.begin(); it != poly.coeff.end(); ++it)
			insertMonomial(it->first, it->second, &ret); // add all elements of the second operand

		std::set<char> new_vars = this->variables;
		new_vars.insert(poly.variables.begin(), poly.variables.end()); // concat variable set
		return Polynomial(new_vars, ret);
	}
	
	Polynomial<SR> operator*(const Polynomial<SR>& poly) const
	{
		Tcoeff ret;
		for (typename Tcoeff::const_iterator it1 = this->coeff.begin(); it1 != this->coeff.end(); ++it1)
		{
			for (typename Tcoeff::const_iterator it2 = poly.coeff.begin(); it2 != poly.coeff.end(); ++it2)
			{
				std::stringstream ss;
				ss << it1->first << it2->first; // non commutative multiplication of variables
				std::string tmp = ss.str();
				insertMonomial(tmp, it1->second * it2->second, &ret); // use semiring multiplication
			}
		}
		std::set<char> new_vars = this->variables;
		new_vars.insert(poly.variables.begin(), poly.variables.end()); // concat variable set
		return Polynomial(new_vars, ret);
	}

	Polynomial<SR> derivative(const char& var) const
	{
		Tcoeff ret;
		// derivative of a polynomial is the sum of all derivated monomials
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			// naive use of product rule in monomialDerivative
			std::string vars = it->first;
			std::list<std::string> derivs = monomialDerivative(var, vars);
			for(std::list<std::string>::const_iterator deriv = derivs.begin(); deriv != derivs.end(); ++deriv)
			{
				insertMonomial(*deriv, it->second, &ret);
			}
		}

		return Polynomial(this->variables, ret);
	}

	Polynomial<SR> derivative(const std::string& vars) const
	{
		Polynomial<SR> ret = *this;
		for(int pos = 0; pos < vars.length(); pos++)
		{
			ret = ret.derivative(vars[pos]);
		}
		return ret;
	}

	int get_degree()
	{
		return this->degree;
	}

	// some semiring functions
	Polynomial<SR> star() const
	{
		// TODO?
		return (*this);
	}

	static Polynomial<SR> null()
	{
		return Polynomial();
	}

	bool is_idempotent() const
	{
		return 0;
	}

	bool is_commutative() const
	{
		return 0;
	}

	std::string string() const
	{
		std::stringstream ss;
		for (typename Tcoeff::const_iterator it = coeff.begin(); it != coeff.end(); ++it) {
			if(it != coeff.begin())
				ss << " + ";
			SR tmp = it->second;
			ss << tmp << it->first;
		}
		return ss.str();
	}
};

#endif
