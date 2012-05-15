#ifndef POLYNOME_H
#define POLYNOME_H

#include <map>
#include <set>
#include <iostream>
#include <string>
#include <sstream>
#include <assert.h>

template <typename SR>
class Polynome
{
private:
	std::set<char> variables;
	typedef std::map<std::string, SR> Tcoeff;
	Tcoeff coeff;
public:
	// empty polynome
	Polynome()
	{
	};

	Polynome(std::set<char> variables, Tcoeff coeff)
	{
		this->variables = variables;
		this->coeff = coeff;
	};

	Polynome<SR> operator+(const Polynome<SR>& poly) const
	{
		Tcoeff ret;
		for (typename Tcoeff::const_iterator it = this->coeff.begin(); it != this->coeff.end(); ++it)
		{
			typename Tcoeff::const_iterator elem = poly.coeff.find(it->first);
			if(elem != poly.coeff.end()) // element found
			{
				ret[it->first] = it->second + elem->second;
			}
			else
			{
				ret[it->first] = it->second;
			}
		} // ret now contains all elements of this polynome excluding exclusive elements of the second operand

		// include them as well
		for (typename Tcoeff::const_iterator it = poly.coeff.begin(); it != poly.coeff.end(); ++it)
		{
			if(this->coeff.find(it->first) == this->coeff.end()) // this "poly.coeff" element is not in "this->coeff"
			{
				ret[it->first] = it->second;
			}
		}
		
		std::set<char> new_vars = this->variables;
		new_vars.insert(poly.variables.begin(), poly.variables.end()); // concat variable set
		return Polynome(new_vars, ret);
	}
	
	Polynome<SR> operator*(const Polynome<SR>& poly) const
	{
		Tcoeff ret;
		for (typename Tcoeff::const_iterator it1 = this->coeff.begin(); it1 != this->coeff.end(); ++it1)
		{
			for (typename Tcoeff::const_iterator it2 = poly.coeff.begin(); it2 != poly.coeff.end(); ++it2)
			{
				std::stringstream ss;
				ss << it1->first << it2->first; // non commutative multiplication of variables
				std::string tmp = ss.str();
				ret[tmp] = it1->second * it2->second; // semiring multiplication
			}
		}
		std::set<char> new_vars = this->variables;
		new_vars.insert(poly.variables.begin(), poly.variables.end()); // concat variable set
		return Polynome(new_vars, ret);
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

template <typename SR>
std::ostream& operator<<(std::ostream& os, Polynome<SR>& poly)
{
	return os << poly.string();
}

#endif
