#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <map>
#include <set>
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include <assert.h>
#include "var.h"
#include "semiring.h"
#include "matrix.h"

template <typename SR>
class Polynomial : public Semiring<Polynomial<SR> >
{
private:
	int degree;
	std::multiset<Var> variables;
	typedef std::map<std::multiset<Var>, SR> Tcoeff;
	Tcoeff coeff;

	// non-commutative version
	std::list<std::multiset<Var> > monomialDerivative(const Var& var, const std::multiset<Var>& vars) const
	{
		std::list<std::multiset<Var> > ret;
		// (uvwx) = u'vwx + uv'wx + uvw'x + uvwx'
		// but uv'wx = 0 for v != var because uv'wx = u0wx = 0
		// so just use all terms uv'wx for v = var -> uv'wx = uwx
		for (std::multiset<Var>::iterator v = vars.begin(); v != vars.end(); ++v)
		{
			if( (*v) == var)
			{
				std::multiset<Var> tmp = vars;
				tmp.erase(v); // uv'wx -> uwx
				ret.push_back(tmp);
			}
		}

		return ret;
	};

	// insert a monomial into coeffs
	static void insertMonomial(const std::multiset<Var>& vars, const SR& coeff, Tcoeff* coeffs)
	{
		if(coeff == SR::null())
			return; // do not save 0 coefficients
		typename Tcoeff::const_iterator elem = coeffs->find(vars);
		if(elem == coeffs->end()) // not in the map yet
			(*coeffs)[vars] = coeff;
		else // there is already a coefficient
			(*coeffs)[vars] = coeffs->at(vars) + coeff;
	};
public:
	// empty polynomial
	Polynomial()
	{
		this->degree = 0;
		std::multiset<Var> vars = {Var("")};
		this->coeff[vars] = SR::null();
	};

	Polynomial(const Polynomial& polynomial)
	{
		this->variables = polynomial.variables;
		this->coeff = polynomial.coeff;
		this->degree = polynomial.degree;
	}

	// create a 'constant' polynomial
	Polynomial(const SR& elem)
	{
		this->degree = 0;
		std::multiset<Var> vars = {Var("")};
		this->coeff.at(vars) = elem;
	}

	Polynomial& operator=(const Polynomial& polynomial)
	{
		this->variables = polynomial.variables;
		this->coeff = polynomial.coeff;
		this->degree = polynomial.degree;
		return (*this);
	}

	Polynomial(std::multiset<Var> variables, Tcoeff coeff)
	{
		this->variables = variables;
		this->coeff = coeff;
		// calculate the degree. will not change because polynomial is read only
		this->degree = 0;
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			this->degree = it->first.size() > this->degree ? it->first.size() : this->degree;
		}

	};

	Polynomial<SR> operator+(const Polynomial<SR>& poly) const
	{
		Tcoeff ret = this->coeff; // fill coefficients with all elements of the first operand
		for (typename Tcoeff::const_iterator it = poly.coeff.begin(); it != poly.coeff.end(); ++it)
			insertMonomial(it->first, it->second, &ret); // add all elements of the second operand

		std::multiset<Var> new_vars = this->variables;
		new_vars.insert(poly.variables.begin(), poly.variables.end()); // concat variable multiset
		return Polynomial(new_vars, ret);
	}
	
	Polynomial<SR> operator*(const Polynomial<SR>& poly) const
	{
		Tcoeff ret;
		for (typename Tcoeff::const_iterator it1 = this->coeff.begin(); it1 != this->coeff.end(); ++it1)
		{
			for (typename Tcoeff::const_iterator it2 = poly.coeff.begin(); it2 != poly.coeff.end(); ++it2)
			{
				std::multiset<Var> tmp = it1->first;
				tmp.insert(it2->first.begin(), it2->first.end()); // commutative multiplication of variables
				insertMonomial(tmp, it1->second * it2->second, &ret); // use semiring multiplication
			}
		}
		std::multiset<Var> new_vars = this->variables;
		new_vars.insert(poly.variables.begin(), poly.variables.end()); // concat variable multiset
		return Polynomial(new_vars, ret);
	}

	friend Polynomial<SR> operator*(const SR& elem, const Polynomial<SR>& polynomial)
	{
		Tcoeff ret;
		for (typename Tcoeff::const_iterator it = polynomial.coeff.begin(); it != polynomial.coeff.end(); ++it)
		{
			insertMonomial(it->first, elem * it->second, &ret);
		}
		return Polynomial(polynomial.vars, ret);
	}

	// convert the given matrix to a matrix containing polynomials
	static Matrix<Polynomial<SR> > convert(const Matrix<SR>& mat)
	{
		std::vector<Polynomial<SR> > ret;
		for(int i=0; i<mat.getColumns() * mat.getRows(); ++i)
		{
			// create constant polynomials
			ret.push_back(Polynomial(mat.getElements().at(i)));
		}

		return Matrix<Polynomial<SR> >(mat.getColumns(), mat.getRows(), ret);
	}

	Polynomial<SR> derivative(const Var& var) const
	{
		Tcoeff ret;
		// derivative of a polynomial is the sum of all derivated monomials
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			// naive use of product rule in monomialDerivative
			std::multiset<Var> vars = it->first;
			std::list<std::multiset<Var> > derivs = monomialDerivative(var, vars);
			for(std::list<std::multiset<Var> >::const_iterator deriv = derivs.begin(); deriv != derivs.end(); ++deriv)
			{
				insertMonomial(*deriv, it->second, &ret);
			}
		}

		// if var was not in the polynomial return the null polynomial
		if(ret.empty())
			return Polynomial();

		return Polynomial(this->variables, ret);
	}

	Polynomial<SR> derivative(const std::multiset<Var>& vars) const
	{
		Polynomial<SR> ret = *this;
		for(std::multiset<Var>::const_iterator var = vars.begin(); var != vars.end(); ++var)
		{
			ret = ret.derivative(*var);
		}
		return ret;
	}

	static Matrix<Polynomial<SR> > jacobian(const std::vector<Polynomial<SR> >& polynomials, const std::multiset<Var>& variables)
	{
		std::vector<Polynomial<SR> > ret;
		for(typename std::vector<Polynomial<SR> >::const_iterator poly = polynomials.begin(); poly != polynomials.end(); ++poly)
		{
			for(std::multiset<Var>::const_iterator var = variables.begin(); var != variables.end(); ++var)
			{
				ret.push_back(poly->derivative(*var));
			}
		}
		return Matrix<Polynomial<SR> >(variables.size(), polynomials.size(), ret);
	};

	Matrix<Polynomial<SR> > hessian() const
	{
		std::vector<Polynomial<SR> > ret;
		for(std::multiset<Var>::const_iterator var2 = this->variables.begin(); var2 != this->variables.end(); ++var2)
		{
			Polynomial<SR> tmp = this->derivative(*var2);
			for(std::multiset<Var>::const_iterator var1 = this->variables.begin(); var1 != this->variables.end(); ++var1)
			{
				ret.push_back(tmp.derivative(*var1));
			}
		}
		return Matrix<Polynomial<SR> >(this->variables.size(), this->variables.size(), ret);
	}

	SR eval(std::map<Var,SR>& vars) const
	{
		SR ret = SR::null();
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			std::multiset<Var> variables = it->first;
			SR elem = it->second;
			for(std::multiset<Var>::const_iterator var = variables.begin(); var != variables.end(); ++var)
			{
				SR foo = vars.find((*var))->second;
				elem = elem * foo;
			}
			ret = ret + elem;
		}

		return ret;
	}

	Polynomial<SR> eval(std::map<Var,Polynomial<SR> >& vars) const
	{
		Polynomial<SR> ret = Polynomial<SR>::null();
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			std::multiset<Var> variables = it->first;
			Polynomial<SR> elem = it->second;
			for(std::multiset<Var>::const_iterator var = variables.begin(); var != variables.end(); ++var)
			{
				Polynomial<SR> foo = vars.find((*var))->second;
				elem = elem * foo;
			}
			ret = ret + elem;
		}

		return ret;
	}

	static Matrix<SR> eval(Matrix<Polynomial<SR> > polys, std::map<Var,SR> vars)
	{
		std::vector<Polynomial<SR> > polynomials = polys.getElements();
		std::vector<SR> ret;
		for(int i = 0; i < polys.getRows()*polys.getColumns(); i++)
		{
			ret.push_back(polynomials[i].eval(vars));
		}
		return Matrix<SR>(polys.getColumns(), polys.getRows(), ret);
	}

	static Matrix<Polynomial<SR> > eval(Matrix<Polynomial<SR> > polys, std::map<Var,Polynomial<SR> > vars)
	{
		std::vector<Polynomial<SR> > polynomials = polys.getElements();
		std::vector<Polynomial<SR> > ret;
		for(int i = 0; i < polys.getRows()*polys.getColumns(); i++)
		{
			ret.push_back(polynomials[i].eval(vars));
		}
		return Matrix<Polynomial<SR> >(polys.getColumns(), polys.getRows(), ret);
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

	static Polynomial<SR> one()
	{
		return Polynomial(SR::one());
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
			ss << tmp;
			//if(it->first != "")
			//	ss << "*";
			ss << it->first;
		}
		return ss.str();
	}
};

template <typename SR>
std::ostream& operator<<(std::ostream& os, const std::map<std::multiset<Var>, SR>& values)
{
	for(typename std::map<std::multiset<Var>, SR >::const_iterator value = values.begin(); value != values.end(); ++value)
	{
		os << value->first << "→" << value->second << ";";
	}
	return os;
}

#endif
