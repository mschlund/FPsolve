#include <assert.h>
#include <string>
#include "commutativeRExp.h"
#include <iostream>

// this creates the empty set
CommutativeRExp::CommutativeRExp()
{
}

CommutativeRExp::CommutativeRExp(Var var)
{
	this->sets = {std::make_pair<std::multiset<Var>, std::set<Var> >({var},{})};
}

CommutativeRExp::CommutativeRExp(std::set<linset> sets)
{
	this->sets = sets;
}

CommutativeRExp::CommutativeRExp(const CommutativeRExp& expr)
{
	this->sets = expr.sets;
}

CommutativeRExp::~CommutativeRExp()
{
	if(elem_epsilon != 0)
	{
		delete elem_epsilon;
		elem_epsilon = 0;
	}
	/* this segfaults?!
	if(elem_null != 0)
	{
		delete elem_null;
		elem_null = 0;
	}*/
	if(elem_one != 0)
	{
		delete elem_one;
		elem_one = 0;
	}
}

// union operator
CommutativeRExp CommutativeRExp::operator +(const CommutativeRExp& expr) const
{
	return CommutativeRExp(ExpUnion(this->sets, expr.sets));
}

// concatenate all expressions from first set with all expressions of the second set
CommutativeRExp CommutativeRExp::operator *(const CommutativeRExp& expr) const
{
	return CommutativeRExp(ExpConcat(this->sets,expr.sets));
}


// create a union of two expression-sets
std::set<CommutativeRExp::linset> CommutativeRExp::ExpUnion(std::set<linset> e1, std::set<linset> e2)
{
	std::set<linset> ret = e1;
	ret.insert(e2.begin(),e2.end());
	return ret;
}

// concatenate two expressions
CommutativeRExp::linset CommutativeRExp::ExpConcat(linset e1, linset e2)
{
	linset ret;
	if( (e1.first.empty() && e1.second.empty()) || (e2.first.empty() && e2.second.empty()) ) // one of the expressions is {}
	{
		return ret; // return {}
	}

	if(e1 == epsilon()) // ε*x = x
	{
		ret = e2;
	}
	else if(e2 == epsilon()) // x*ε = x
	{
		ret = e1;
	}
	else // x*y = xy
	{
		ret = std::make_pair(e1.first,e1.second);
		ret.first.insert(e2.first.begin(),e2.first.end());
		ret.second.insert(e2.second.begin(),e2.second.end());
	}

	return ret;
}

// concatenate two expression sets
std::set<CommutativeRExp::linset> CommutativeRExp::ExpConcat(std::set<linset> e1, std::set<linset> e2)
{
	std::set<linset> ret;
	for(std::set<linset>::const_iterator set1 = e1.begin(); set1 != e1.end(); ++set1)
	{
		for(std::set<linset>::const_iterator set2 = e2.begin(); set2 != e2.end(); ++set2)
		{
			ret.insert(CommutativeRExp::ExpConcat(*set1,*set2));
		}
	}
	return ret;
}

// handles the starring of one regexp
std::set<CommutativeRExp::linset> CommutativeRExp::ExpStar(linset exp)
{
	// identities:
	// (1) (x*)* = x*
	// (2) (x+y)* = x*y* (not covered in this method)
	// (3) (xy*)* = ε + xx*y*

	std::set<linset> ret;
	if(exp == linset({},{})) // {}* = ε
	{
		ret.insert(epsilon());
		return ret;
	}

	if(exp.first.empty()) // (x*)* = x* (id)
	{
		ret.insert(exp);
		return ret;
	}

	if(exp.second.empty()) // (x)* = x*
	{
		std::set<Var> tmp(exp.first.begin(),exp.first.end());
		ret.insert(linset({},tmp));
		return ret;
	}

	// (xy*)* = ε + xx*y*
	ret.insert(epsilon());
	std::set<Var> tmp(exp.first.begin(),exp.first.end()); // x -> x*
	tmp.insert(exp.second.begin(),exp.second.end()); // concatenate x* with y* -> x*y*
	ret.insert(linset(exp.first,tmp)); // xx*y*;
	return ret;
}

// TODO
bool CommutativeRExp::operator ==(const CommutativeRExp& expr) const
{
	return this->sets == expr.sets;
}

// star the whole RExp
CommutativeRExp CommutativeRExp::star() const
{
	// (1) (x*)* = x*
	// (2) (x+y)* = x*y*
	// (3) (xy*)* = 1 + xx*y*
	// (?) (x+yz*)* = x*(yz*)* = x*(1+yy*z*) = x*+yx*y*z

	// {}* = ε ?
	if(*this == CommutativeRExp::null())
		return CommutativeRExp::one();

	std::set<linset> ret = {epsilon()};
	for(std::set<linset>::const_iterator set1 = this->sets.begin(); set1 != this->sets.end(); ++set1)
	{
		ret = ExpConcat(ret, ExpStar(*set1));
	}

	return CommutativeRExp(ret);
}

CommutativeRExp::linset CommutativeRExp::epsilon()
{
	if(!CommutativeRExp::elem_epsilon)
		CommutativeRExp::elem_epsilon = new linset(std::multiset<Var>({Var("ε")}),{});
	return *CommutativeRExp::elem_epsilon;
}

CommutativeRExp CommutativeRExp::null()
{
	if(!CommutativeRExp::elem_null)
		CommutativeRExp::elem_null = new CommutativeRExp();
	return *CommutativeRExp::elem_null;
}

CommutativeRExp CommutativeRExp::one()
{
	if(!CommutativeRExp::elem_one)
		CommutativeRExp::elem_one = new CommutativeRExp(Var("ε"));
	return *CommutativeRExp::elem_one;
}

std::string CommutativeRExp::string() const
{
	std::stringstream ss;
	ss << "(";
	for(std::set<linset>::const_iterator set = this->sets.begin(); set != this->sets.end(); ++set)
	{
		if(set != this->sets.begin())
			ss << "+";
		for(std::multiset<Var>::const_iterator var = set->first.begin(); var != set->first.end(); ++var)
			ss << (*var);
		for(std::multiset<Var>::const_iterator var = set->second.begin(); var != set->second.end(); ++var)
			ss << (*var) << "*";
	}
	ss << ")";
	return ss.str();
}

bool CommutativeRExp::is_idempotent = true;
bool CommutativeRExp::is_commutative = true;
CommutativeRExp::linset* CommutativeRExp::elem_epsilon = 0;
CommutativeRExp* CommutativeRExp::elem_null = 0;
CommutativeRExp* CommutativeRExp::elem_one = 0;
