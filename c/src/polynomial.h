#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <map>
#include <unordered_map>
#include <set>
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include <initializer_list>
#include <assert.h>
#include "var.h"
#include "semiring.h"
#include "matrix.h"
#include "free-semiring.h"

//FIXME: do not generate and return new objects when performing operations like +,* ... inefficient :)


template <typename SR>
class Monomial
{
private:
	std::multiset<Var> variables;
	SR coeff;

	// private constructor to not leak the internal data structure
	Monomial(SR coeff, std::multiset<Var> variables)
	{
		this->coeff = coeff;
		if( !(coeff == SR::null())) // we only need to save the variables in this case
			this->variables = variables;
	}
public:
	// constant monomial coeff
	Monomial(SR coeff)
	{
		this->coeff = coeff;
	}

	Monomial(SR coeff, std::initializer_list<Var> variables)
	{
		this->coeff = coeff;
		if( !(coeff == SR::null())) // we only need to save the variables in this case
			this->variables = variables;
	}

	// add the coefficients of two monomials if there variables are equal
	Monomial operator+(const Monomial& monomial) const
	{
		assert(this->variables == monomial.variables);
		return Monomial(this->coeff + monomial.coeff, this->variables);
	}

	// multiply two monomials
	Monomial operator*(const Monomial& monomial) const
	{
		std::multiset<Var> variables = this->variables;
		// "add" the variables from one to the other monomial
		variables.insert(monomial.variables.begin(), monomial.variables.end());
		return Monomial(this->coeff * monomial.coeff, variables);
	}

	// multiply a monomial with a variable
	Monomial operator*(const Var& var) const
	{
		std::multiset<Var> variables = this->variables;
		// "add" the variables from one to the other monomial
		variables.insert(var);
		return Monomial(this->coeff, variables);
	}

	// commutative version of derivative
	Monomial derivative(const Var& var) const
	{
		// count number of occurences of var in variables
		int count = this->variables.count(var);

		// variable is not in variables, derivative is null
		if(count == 0)
			return Monomial(SR::null());

		// remove one of these by removing the first of them and then "multiply"
		// the coefficient with count
		std::vector<Var> variables(this->variables.begin(), this->variables.end());
		SR coeff = this->coeff;
		for(unsigned int i=0; i<variables.size(); ++i)
		{
			if(variables[i] == var)
			{
				variables.erase(variables.begin()+i);
				for(int j = 0; j<count-1; ++j)
					coeff = coeff + this->coeff;
				break;
			}
		}
		std::multiset<Var> result(variables.begin(), variables.end());
		return Monomial(coeff, result);
	}

	// evaluate the monomial at the position values
	SR eval(const std::map<Var, SR>& values) const
	{
		SR elem = this->coeff;

		for(std::multiset<Var>::const_iterator v_it = this->variables.begin(); v_it != this->variables.end(); ++v_it)
		{
			auto e = values.find(*v_it);
			assert(e != values.end()); // all variables should be in values
			SR foo = e->second;
			elem = elem * foo;
		}

		return elem;
	}

	// substitute variables with other variables
	Monomial<SR> subst(const std::map<Var, Var>& mapping) const
	{
		SR coeff = this->coeff;
		std::multiset<Var> variables = this->variables;

		for(std::map<Var, Var>::const_iterator m_it = mapping.begin(); m_it != mapping.end(); ++m_it)
		{
			int count = variables.count(m_it->first);
			variables.erase( (*m_it).first ); // erase all occurences
			for(int i = 0; i < count; ++i)
				variables.insert( (*m_it).second ); // and "replace" them
		}

		return Monomial(coeff, variables);
	}

	// convert this monomial to an element of the free semiring
	FreeSemiring make_free(std::unordered_map<SR,FreeSemiring,SR>* valuation) const
	{
		FreeSemiring result = FreeSemiring::one();
		for(std::multiset<Var>::const_iterator v_it = this->variables.begin(); v_it != this->variables.end(); ++v_it)
		{
			result = result * FreeSemiring(*v_it);
		}

		// change the SR element to a constant in the free semiring
		typename std::unordered_map<SR,FreeSemiring>::const_iterator elem = valuation->find(this->coeff);
		if(elem == valuation->end()) // this is a new SR element
		{
			// map 'zero' and 'one' element to respective free semiring element
			if(this->coeff == SR::null())
			{
				valuation->insert(valuation->begin(), std::pair<SR,FreeSemiring>(this->coeff,FreeSemiring::null()));
				result = FreeSemiring::null() * result;
			}
			else if(this->coeff == SR::one())
			{
				valuation->insert(valuation->begin(), std::pair<SR,FreeSemiring>(this->coeff,FreeSemiring::one()));
				result = FreeSemiring::one() * result;
			}
			else
			{
				// use a fresh constant - the constructor of Var() will do this
				Var tmp_var;
				FreeSemiring tmp(tmp_var);
				valuation->insert(valuation->begin(), std::pair<SR,FreeSemiring>(this->coeff,tmp));
				result = tmp * result;
			}
		}
		else // this is an already known element
		{
			result = (elem->second) * result;
		}
		return result;
	}

	// a monomial is smaller than another monomial if the variables are smaller
	bool operator<(const Monomial& monomial) const
	{
		return this->variables < monomial.variables;
	}

	// a monomial is equal to another monomial if the variables are equal
	// warning: the coefficient will not be regarded
	bool operator==(const Monomial& monomial) const
	{
		return this->variables == monomial.variables;
	}

	int get_degree() const
	{
		return this->variables.size();
	}

	std::string string() const
	{
		std::stringstream ss;
		ss << this->coeff;
		if( !(this->coeff == SR::null() || this->variables.empty()))
		{
			ss << "*" << this->variables;
		}
		return ss.str();
	}
};

template <typename SR>
class Polynomial : public Semiring<Polynomial<SR> >
{
private:
	int degree;
	std::set<Monomial<SR> > monomials;

	// private constructor to hide the internal data structure
	Polynomial(const std::set<Monomial<SR> >& monomials)
	{
		this->monomials = monomials;
		this->degree = 0;
		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			this->degree = (*m_it).get_degree() > this->degree ? (*m_it).get_degree() : this->degree;
		}
	}

	static Polynomial<SR>* elem_null;
	static Polynomial<SR>* elem_one;
public:
	// empty polynomial
	Polynomial()
	{
		this->degree = 0;
	};

	Polynomial(std::initializer_list<Monomial<SR> > monomials)
	{
		this->monomials = monomials;
		this->degree = 0;
		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			this->degree = (*m_it).get_degree() > this->degree ? (*m_it).get_degree() : this->degree;
		}
	}

	Polynomial(const Polynomial& polynomial)
	{
		this->monomials = polynomial.monomials;
		this->degree = polynomial.degree;
	}

	// create a 'constant' polynomial
	Polynomial(const SR& elem)
	{
		this->monomials = {Monomial<SR>(elem,{})};
		this->degree = 0;
	}

	Polynomial& operator=(const Polynomial& polynomial)
	{
		this->monomials = polynomial.monomials;
		this->degree = polynomial.degree;
		return (*this);
	}

	Polynomial<SR> operator+(const Polynomial<SR>& poly) const
	{
		// TODO: check this!
		std::set<Monomial<SR> > monomials = this->monomials;
		for(typename std::set<Monomial<SR> >::const_iterator m_it = poly.monomials.begin(); m_it != poly.monomials.end(); ++m_it)
		{
			// check if "same" monomial is already in the set
			typename std::set<Monomial<SR> >::const_iterator mon = monomials.find(*m_it);
			if(mon == monomials.end()) // this is not the case
				monomials.insert(*m_it); // just insert it
			else // monomial with the same variables found
				monomials.insert( (*mon) + (*m_it) ); // then add both of them and overwrite the old one
		}
		return Polynomial(monomials);
	}
	
	Polynomial<SR> operator*(const Polynomial<SR>& poly) const
	{
		std::set<Monomial<SR> > monomials;
		// iterate over both this and the poly polynomial
		for(typename std::set<Monomial<SR> >::const_iterator m_it1 = this->monomials.begin(); m_it1 != this->monomials.end(); ++m_it1)
		{
			for(typename std::set<Monomial<SR> >::const_iterator m_it2 = poly.monomials.begin(); m_it2 != poly.monomials.end(); ++m_it2)
			{
				monomials.insert( (*m_it1) * (*m_it2) ); // multiply them and insert them to the result set
			}
		}

		return Polynomial(monomials);
	}

	// multiplying a polynomial with a variable
	Polynomial<SR> operator*(const Var& var) const
	{
		std::set<Monomial<SR> > monomials;
		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			monomials.insert( (*m_it) * var );
		}

		return Polynomial(monomials);
	}

	friend Polynomial<SR> operator*(const SR& elem, const Polynomial<SR>& polynomial)
	{
		std::set<Monomial<SR> > monomials;

		for(typename std::set<Monomial<SR> >::const_iterator m_it = polynomial.monomials.begin(); m_it != polynomial.monomials.end(); ++m_it)
		{
			monomials.insert(elem * (*m_it));
		}

		return Polynomial(monomials);
	}

	// convert the given matrix to a matrix containing polynomials
	// TODO: needed?
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
		std::set<Monomial<SR> > monomials;

		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			if(SR::is_commutative) // TODO: check if compiler is optimizing this out
			{
				// take the derivative of m_it and add it to the result set
				Monomial<SR> derivative = (*m_it).derivative(var);
				typename std::set<Monomial<SR> >::const_iterator monomial = monomials.find(derivative);
				if(monomial == monomials.end()) // TODO: think about this and remove if not needed
				{
					monomials.insert(derivative);
				}
				else
				{
					Monomial<SR> tmp = (*monomial) + derivative;
					monomials.erase(monomial); // remove
					monomials.insert(tmp); // and insert the updated version
				}
			}
			else // non-commutative case
			{
				assert(false);	// TODO: not implemented yet
			}
		}

		if(monomials.empty()) // TODO: save variables explicit in this class and check if var is in vars
			return Polynomial();
		else
			return Polynomial(monomials);
	}

	Polynomial<SR> derivative(const std::vector<Var>& vars) const
	{
		Polynomial<SR> polynomial = *this; // copy the polynomial
		for(typename std::vector<Var>::const_iterator var = vars.begin(); var != vars.end(); ++var)
		{
			polynomial = polynomial.derivative(*var);
		}
		return polynomial;
	}

	static Matrix<Polynomial<SR> > jacobian(const std::vector<Polynomial<SR> >& polynomials, const std::vector<Var>& variables)
	{
		std::vector<Polynomial<SR> > ret;
		for(typename std::vector<Polynomial<SR> >::const_iterator poly = polynomials.begin(); poly != polynomials.end(); ++poly)
		{
			for(std::vector<Var>::const_iterator var = variables.begin(); var != variables.end(); ++var)
			{
				ret.push_back(poly->derivative(*var));
			}
		}
		return Matrix<Polynomial<SR> >(variables.size(), polynomials.size(), ret);
	};

	// TODO: i need the variables in this function!
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

	SR eval(const std::map<Var,SR>& values) const
	{
		SR result = SR::null();
		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			Monomial<SR> foo = (*m_it); // TODO: collapse the lines
			SR elem = foo.eval(values);
			result = result + elem;
		}
		return result;
	}

	// substitute variables with other variables
	Polynomial<SR> subst(const std::map<Var, Var>& mapping) const
	{
		std::set<Monomial<SR> > monomials;

		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
			monomials.insert((*m_it).subst(mapping));

		return Polynomial<SR>(monomials);
	}

	// TODO: is this method really needed?
	/*Polynomial<SR> eval(const std::map<Var,Polynomial<SR> >& values) const
	{
		Polynomial<SR> result = Polynomial<SR>::null;
		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			Monomial<SR> foo = (*m_it);
			Polynomial<SR> elem = foo.eval(values);
			result = result + elem;
		}

		return result;
	}*/

	static Matrix<SR> eval(const Matrix<Polynomial<SR> >& polys, const std::map<Var,SR>& values)
	{
		std::vector<Polynomial<SR> > polynomials = polys.getElements();
		std::vector<SR> ret;
		for(int i = 0; i < polys.getRows()*polys.getColumns(); i++)
		{
			ret.push_back(polynomials[i].eval(values));
		}
		return Matrix<SR>(polys.getColumns(), polys.getRows(), ret);
	}

	static Matrix<Polynomial<SR> > eval(Matrix<Polynomial<SR> > polys, std::map<Var,Polynomial<SR> > values)
	{
		std::vector<Polynomial<SR> > polynomials = polys.getElements();
		std::vector<Polynomial<SR> > ret;
		for(int i = 0; i < polys.getRows()*polys.getColumns(); i++)
		{
			ret.push_back(polynomials[i].eval(values));
		}
		return Matrix<Polynomial<SR> >(polys.getColumns(), polys.getRows(), ret);
	}

	// convert this polynomial to an element of the free semiring. regard the valuation map
	// which can already define a conversion from the SR element to a free SR constant
	// the valuation map is extended in this function
	FreeSemiring make_free(std::unordered_map<SR,FreeSemiring, SR>* valuation)
	{
		if(!valuation)
			valuation = new std::unordered_map<SR, FreeSemiring, SR>();

		FreeSemiring result = FreeSemiring::null();
		// convert this polynomial by adding all converted monomials
		for(typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			result = result + m_it->make_free(valuation);
		}

		return result;
	}

	// convert this matrix of polynomials to a matrix with elements of the free semiring
	static Matrix<FreeSemiring> make_free(const Matrix<Polynomial<SR> >& polys, std::unordered_map<SR, FreeSemiring, SR>* valuation)
	{
		std::vector<Polynomial<SR> > polynomials = polys.getElements();
		std::vector<FreeSemiring> ret;
		if(!valuation)
			valuation = new std::unordered_map<SR, FreeSemiring, SR>();

		for(int i = 0; i < polys.getRows()*polys.getColumns(); i++)
		{
			ret.push_back(polynomials[i].make_free(valuation));
		}
		return Matrix<FreeSemiring>(polys.getColumns(), polys.getRows(), ret);
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

	static Polynomial<SR> const null()
	{
		if(!Polynomial::elem_null)
			Polynomial::elem_null = new Polynomial(SR::null());
		return *Polynomial::elem_null;
	}

	static Polynomial<SR> const one()
	{
		if(!Polynomial::elem_one)
			Polynomial::elem_one = new Polynomial(SR::one());
		return *Polynomial::elem_one;
	}

	static bool is_idempotent;
	static bool is_commutative;

	std::string string() const
	{
		std::stringstream ss;
		for (typename std::set<Monomial<SR> >::const_iterator m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
		{
			if(m_it != this->monomials.begin())
				ss << " + ";
			ss << (*m_it);
		}

		return ss.str();
	}
};

template <typename SR> bool Polynomial<SR>::is_commutative = false;
template <typename SR> bool Polynomial<SR>::is_idempotent = false;
// initialize pointers
template <typename SR> Polynomial<SR>* Polynomial<SR>::elem_null = 0;
template <typename SR> Polynomial<SR>* Polynomial<SR>::elem_one = 0;

template <typename SR>
std::ostream& operator<<(std::ostream& os, const Monomial<SR>& monomial)
{
	return 	os << monomial.string();
}

template <typename SR>
std::ostream& operator<<(std::ostream& os, const std::map<Var, SR>& values)
{
	for(typename std::map<Var, SR>::const_iterator value = values.begin(); value != values.end(); ++value)
	{
		os << value->first << "→" << value->second << ";";
	}
	return os;
}
template <typename SR>
std::ostream& operator<<(std::ostream& os, const std::map<Var, Polynomial<SR> >& values)
{
	for(typename std::map<Var, Polynomial<SR> >::const_iterator value = values.begin(); value != values.end(); ++value)
	{
		os << value->first << "→" << value->second << ";";
	}
	return os;
}

#endif
