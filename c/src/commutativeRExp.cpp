#include <assert.h>
#include <string>
#include "commutativeRExp.h"

// this creates the empty set
CommutativeRExp::CommutativeRExp()
{
	this->type = Empty;
}

CommutativeRExp::CommutativeRExp(VarPtr var)
{
	this->type = Element;
	this->elem = Var::getVar(var);
}

CommutativeRExp::CommutativeRExp(enum optype type, std::shared_ptr<std::set<CommutativeRExp>> seta)
{
	this->type = type;
	if(this->type != Addition)
		assert(false); // should not be called with this constructor...
	this->seta = seta;
}

CommutativeRExp::CommutativeRExp(enum optype type, std::shared_ptr<std::multiset<CommutativeRExp>> setm)
{
	this->type = type;
	if(this->type != Multiplication)
		assert(false); // should not be called with this constructor...
	this->setm = setm;
}

CommutativeRExp::CommutativeRExp(enum optype type, std::shared_ptr<CommutativeRExp> rexp)
{
	this->type = type;
	if(this->type != Star)
		assert(false); // should not be called with this constructor...
	this->rexp = rexp;
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
}

CommutativeRExp::~CommutativeRExp()
{
}

// union operator
CommutativeRExp CommutativeRExp::operator +(const CommutativeRExp& expr) const
{
	std::shared_ptr<std::set<CommutativeRExp> > retset(new std::set<CommutativeRExp>());

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
	std::shared_ptr<std::multiset<CommutativeRExp> > retset(new std::multiset<CommutativeRExp>());

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
	// compare the expression type
	if(this->type < rhs.type)
		return true;
	else if(this->type > rhs.type)
		return false;
	else // same type, compare content
	{
		switch(this->type)
		{
			case Element: // reduce to element comparison
				return this->elem < rhs.elem;
			case Addition: // length- then lexicographic-order based on elements of seta
				if(this->seta->size() < rhs.seta->size())
					return true;

				for(auto l_it = this->seta->begin(), r_it = rhs.seta->begin();
					(l_it != this->seta->end() ) && (r_it != rhs.seta->end());
					++l_it, ++r_it )
				{
					if((*l_it) < (*r_it))
						return true;
					else
					{
						if( (*l_it) == (*r_it)) // elements are equal
							continue; // check the next element
						else
							return false; // 'this' is bigger
					}
				}
				return false; // all elements are equal
			case Multiplication: // length- then lexicographic-order based on elements of setm
				if(this->setm->size() < rhs.setm->size())
					return true;

				for(auto l_it = this->setm->begin(), r_it = rhs.setm->begin();
					(l_it != this->setm->end() ) && (r_it != rhs.setm->end());
					++l_it, ++r_it )
				{
					if((*l_it) < (*r_it))
						return true;
					else
					{
						if( (*l_it) == (*r_it)) // elements are equal
							continue; // check the next element
						else
							return false; // 'this' is bigger
					}
				}
				return false; // all elements are equal
			case Star: // reduce to element comparison
				return this->rexp < rhs.rexp;
			case Empty: // empty expressions are always equal
				return false;
		}
	}
}

bool CommutativeRExp::operator ==(const CommutativeRExp& expr) const
{
	if(this->type != expr.type)
		return false;
	else // same type, compare content
	{
		switch(this->type)
		{
			case Element:
				return this->elem == expr.elem;
			case Addition: // all elements have to be equal
				for(auto l_it = this->seta->begin(), r_it = expr.seta->begin();
					(l_it != this->seta->end() ) && (r_it != expr.seta->end());
					++l_it, ++r_it )
				{
					if( !((*l_it) == (*r_it)) ) // at least one element is different
						return false;
				}
				return true; // all elements are equal
			case Multiplication: // all elements have to be equal
				for(auto l_it = this->setm->begin(), r_it = expr.setm->begin();
					(l_it != this->setm->end() ) && (r_it != expr.setm->end());
					++l_it, ++r_it )
				{
					if( !((*l_it) == (*r_it) )) // at least one element is different
						return false;
				}
				return true; // all elements are equal
			case Star: // reduce to element comparison
				return *this->rexp == *expr.rexp;
			case Empty: // empty expressions are always equal
				return true;
		}
	}
}

// star the whole RExp
CommutativeRExp CommutativeRExp::star() const
{
	if(this->type == Star)
		return *this; // absorb the star
	else if (this->type == Empty)
		return one(); // {}* = ε
	else
		return CommutativeRExp(Star, std::shared_ptr<CommutativeRExp>(new CommutativeRExp(*this))); // star the element
}

CommutativeRExp CommutativeRExp::null()
{
	if(!CommutativeRExp::elem_null)
		CommutativeRExp::elem_null = std::shared_ptr<CommutativeRExp>(new CommutativeRExp());
	return *CommutativeRExp::elem_null;
}

CommutativeRExp CommutativeRExp::one()
{
	if(!CommutativeRExp::elem_one)
		CommutativeRExp::elem_one = std::shared_ptr<CommutativeRExp>(new CommutativeRExp(Var::getVar("ε")));
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
		ss << this->elem;
	else if(this->type == Addition)
		ss << "(" << *this->seta << ")";
	else if(this->type == Multiplication)
		ss << "(" << *this->setm << ")";
	else if(this->type == Star)
		ss << *this->rexp << "*";
	else if(this->type == Empty)
		ss << "{}";
	return ss.str();
}

std::string CommutativeRExp::string() const
{
	return this->generateString();
}

bool CommutativeRExp::is_idempotent = true;
bool CommutativeRExp::is_commutative = true;
std::shared_ptr<CommutativeRExp> CommutativeRExp::elem_null;
std::shared_ptr<CommutativeRExp> CommutativeRExp::elem_one;
