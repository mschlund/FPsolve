#include <assert.h>
#include <string>
#include "commutativeRExp.h"

// this creates the empty set
CommutativeRExp::CommutativeRExp()
{
	this->type = Empty;
}

CommutativeRExp::CommutativeRExp(Var var)
{
	this->type = Element;
	this->elem = new Var(var);
	this->str = generateString();
}

CommutativeRExp::CommutativeRExp(enum optype type, std::set<CommutativeRExp>* seta)
{
	this->type = type;
	if(this->type == Addition)
		this->seta = seta;
	else
		assert(false); // should not be called with this constructor...
	this->str = generateString();
}

CommutativeRExp::CommutativeRExp(enum optype type, std::multiset<CommutativeRExp>* setm)
{
	this->type = type;
	if(this->type == Multiplication)
		this->setm = setm;
	else
		assert(false); // should not be called with this constructor...
	this->str = generateString();
}

CommutativeRExp::CommutativeRExp(enum optype type, const CommutativeRExp& rexp)
{
	this->type = type;
	if(this->type == Star)
		this->rexp = new CommutativeRExp(rexp);
	else
		assert(false); // should not be called with this constructor...
	this->str = generateString();
}

CommutativeRExp::CommutativeRExp(const CommutativeRExp& expr)
{
	this->type = expr.type;
	if(this->type == Element)
		this->elem = expr.elem;
	else if(this->type == Addition)
		this->seta = expr.seta;
	else if(this->type == Multiplication)
		this->setm = expr.setm;
	else if(this->type == Star)
		this->rexp = expr.rexp;
	else if(this->type == Empty)
		;// Do nothing
	else
		assert(false);
	this->str = generateString();
}

CommutativeRExp::~CommutativeRExp()
{
/*	if(elem_null != 0)
	{
		delete elem_null;
		elem_null = 0;
	}
	if(elem_one != 0)
	{
		delete elem_one;
		elem_one = 0;
	}
	if(this->type == Element)
	{
		delete elem;
	}
	if(this->type == Addition )
	{
		delete seta;
	}
	if(this->type == Multiplication)
	{
		delete setm;
	}
	if(this->type == Star)
	{
		delete rexp;
	}
*/
}

// union operator
CommutativeRExp CommutativeRExp::operator +(const CommutativeRExp& expr) const
{
	std::set<CommutativeRExp>* retset= new std::set<CommutativeRExp>();

	if(this->type == Element || this->type == Multiplication || this->type == Star)
		retset->insert(*this);
	else if(this->type == Addition)
		retset->insert(this->seta->begin(), this->seta->end());
	else if(this->type == Empty)
		return expr; // {} + x = x
	else
		assert(false); // this should not happen

	if(expr.type == Element || expr.type == Multiplication || expr.type == Star)
		retset->insert(expr);
	else if(expr.type == Addition)
		retset->insert(expr.seta->begin(), expr.seta->end());
	else if(expr.type == Empty)
		return *this; // x + {} = x
	else
		assert(false); // this should not happen

	// degenerated case, both sets have been equal
	if(retset->size() == 1)
		return *this; // so we can just return one of the elements // TODO: check this! quick test shows this is right
	return CommutativeRExp(Addition, retset);
}

// concatenate all expressions from first set with all expressions of the second set
CommutativeRExp CommutativeRExp::operator *(const CommutativeRExp& expr) const
{
	std::multiset<CommutativeRExp>* retset= new std::multiset<CommutativeRExp>();

	if(this->type == Element || this->type == Addition || this->type == Star)
	{
		if(this->type == Element && *this == one())
			return expr; // 1 * x = x
		retset->insert(*this);
	}
	else if(this->type == Multiplication)
		retset->insert(this->setm->begin(), this->setm->end());
	else if(this->type == Empty)
		return *this; // {} * x = {}
	else
		assert(false); // this should not happen

	if(expr.type == Element || expr.type == Addition || expr.type == Star)
	{
		if(expr.type == Element && expr == one())
			return *this; // x * 1 = x
		retset->insert(expr);
	}
	else if(expr.type == Multiplication)
		retset->insert(expr.setm->begin(), expr.setm->end());
	else if(expr.type == Empty)
		return expr; // x * {} = {}
	else
		assert(false); // this should not happen

	return CommutativeRExp(Multiplication, retset);
}

//bool operator <(const CommutativeRExp& lhs, const CommutativeRExp& rhs)
bool CommutativeRExp::operator <(const CommutativeRExp& rhs) const
{
	if(this->type < rhs.type)
		return true;
	else if(this->type > rhs.type)
		return false;
	else // same type
	{
		if(this->type == Element)
			return (*this->elem < *rhs.elem);
		else
			return this->str.compare(rhs.str) < 0;
	}
}

bool CommutativeRExp::operator ==(const CommutativeRExp& expr) const
{
	return this->str.compare(expr.str) == 0;
}

// star the whole RExp
CommutativeRExp CommutativeRExp::star() const
{
	if(this->type == Star)
		return *this; // absorb the star
	else if (this->type == Empty)
		return one(); // {}* = ε
	else
		return CommutativeRExp(Star, *this); // star the element
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

std::ostream& operator<<(std::ostream& os, const std::set<CommutativeRExp>& set)
{
	for(std::set<CommutativeRExp>::const_iterator s = set.begin(); s != set.end(); ++s)
	{
		if(s != set.begin())
			os << "+";
		os << *s;
	}
	return os;
}
std::ostream& operator<<(std::ostream& os, const std::multiset<CommutativeRExp>& set)
{
	for(std::multiset<CommutativeRExp>::const_iterator s = set.begin(); s != set.end(); ++s)
	{
		if(s != set.begin())
			os << ".";
		os << *s;
	}
	return os;
}

std::string CommutativeRExp::generateString() const
{
	std::stringstream ss;
	if(this->type == Element)
		ss << *this->elem;
	else if(this->type == Addition)
		ss << "(" << *this->seta << ")";
	else if(this->type == Multiplication)
		ss << "(" << *this->setm << ")";
	else if(this->type == Star)
		ss << "(" << *this->rexp << ")*";
	return ss.str();
}

std::string CommutativeRExp::string() const
{
	return this->str;
}

bool CommutativeRExp::is_idempotent = true;
bool CommutativeRExp::is_commutative = true;
CommutativeRExp* CommutativeRExp::elem_null = 0;
CommutativeRExp* CommutativeRExp::elem_one = 0;
