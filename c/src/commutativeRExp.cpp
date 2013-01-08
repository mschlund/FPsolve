#include <assert.h>
#include <string>
#include "commutativeRExp.h"

// this creates the empty set
CommutativeRExp::CommutativeRExp()
{
	this->type = Empty;
}

// used by boost::ublas with CommutativeRExp(0) to
// generate a zero Element
CommutativeRExp::CommutativeRExp(int zero)
{
	assert(zero == 0);

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
	if(this->type != Star && this->type != Plus)
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
	else if(this->type == Plus)
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
CommutativeRExp CommutativeRExp::operator +=(const CommutativeRExp& expr)
{
	std::shared_ptr<std::set<CommutativeRExp> > retset(new std::set<CommutativeRExp>());

	if(this->type == Element || this->type == Multiplication || this->type == Star || this->type == Plus)
		retset->insert(*this);
	else if(this->type == Addition)
		retset->insert(this->seta->begin(), this->seta->end());
	else if(this->type == Empty)
	{
		*this = expr;
		return *this; // {} + x = x
	}
	else
		assert(false); // this should not happen

	if(expr.type == Element || expr.type == Multiplication || expr.type == Star || expr.type == Plus)
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
	*this = CommutativeRExp(Addition, retset); // TODO: do not create a new object
	return *this;
}

// try to find a case of xx^* and convert it to x^+
std::shared_ptr<std::multiset<CommutativeRExp>> CommutativeRExp::optimize_starplus(std::shared_ptr<std::multiset<CommutativeRExp>>& set) const
{
	// we have to be careful, we will modify the set_copy!
	// loop over all elements in the set and find the stared elements
	std::shared_ptr<std::multiset<CommutativeRExp>> set_copy(new std::multiset<CommutativeRExp>());
	*set_copy = *set;
	bool changed = false;
	for(auto it = set->begin(); it != set->end(); ++it)
	{
		if(it->type == Star) // we found a stared element
		{
			// distinguish on the nested type
			switch(it->rexp->type)
			{
				case Element:  // easy, just try to find this element on the outer set, fall through case
				case Addition: // addition behaves like an element in the original multiplication set
				{
					auto elem = set_copy->find(*it->rexp);
					if(elem != set_copy->end())
					{
						// create x^+
						auto plus_elem = CommutativeRExp(Plus, it->rexp);
						// delete both x and x^*
						set_copy->erase(*elem);
						set_copy->erase(*it);
						// then insert x^+
						set_copy->insert(plus_elem);
						changed = true;
					}
					break;
				}
				case Multiplication: // more tricky case, because the inner stared elements are distributed in the original multiplication set
				{
					bool found_all_elements = true;
					// for every element in the starred multiplication
					for(auto it2 = it->rexp->setm->begin(); it2 != it->rexp->setm->end(); ++it2)
					{
						auto elem = set_copy->find(*it2);
						if(elem == set_copy->end())
						{
							found_all_elements = false;
							break;
						}
					}
					if(found_all_elements)
					{
						// create x^+
						auto plus_elem = CommutativeRExp(Plus, it->rexp);
						// delete both x and x^*
						for(auto it2 = it->rexp->setm->begin(); it2 != it->rexp->setm->end(); ++it2)
							set_copy->erase(*it2);
						set_copy->erase(*it);
						// then insert x^+
						set_copy->insert(plus_elem);
						changed = true;
					}
					break;
				}
				case Star:  // should not happen, fall through
				case Plus:  // no idea, what to do, fall through
				case Empty: // should not happen
					assert(false); // we should have been optimizing this, debug more!
					break;
			}
		}

	}
	if(changed)
		return set_copy;
	else
		return set;
}

// concatenate all expressions from first set with all expressions of the second set
CommutativeRExp CommutativeRExp::operator *=(const CommutativeRExp& expr)
{
	std::shared_ptr<std::multiset<CommutativeRExp> > retset(new std::multiset<CommutativeRExp>());

	if(this->type == Element || this->type == Addition || this->type == Star || this->type == Plus)
	{
		if(this->type == Element && *this == one())
		{
			*this = expr;
			return *this; // 1 * x = x
		}
		retset->insert(*this);
	}
	else if(this->type == Multiplication)
		retset->insert(this->setm->begin(), this->setm->end());
	else if(this->type == Empty)
		return *this; // {} * x = {}
	else
		assert(false); // this should not happen

	if(expr.type == Element || expr.type == Addition || expr.type == Star || expr.type == Plus)
	{
		if(expr.type == Element && expr == one())
			return *this; // x * 1 = x
		retset->insert(expr);
	}
	else if(expr.type == Multiplication)
		retset->insert(expr.setm->begin(), expr.setm->end());
	else if(expr.type == Empty)
	{
		*this = expr;
		return *this; // x * {} = {}
	}
	else
		assert(false); // this should not happen

	// x(x^*) = x^+
	// check if there was at least one star in the returned multiplication set
	bool star_found = false;
	for(auto it = retset->begin(); it != retset->end(); ++it) // TODO: maybe reverse search order!
	{
		if(it->type == Star)
		{
			star_found = true;
			break;
		}
	}

	if(star_found)
		retset = optimize_starplus(retset);

	*this = CommutativeRExp(Multiplication, retset);
	return *this;
}

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
			case Star: // reduce to element comparison, fall through case
			case Plus:
				return this->rexp < rhs.rexp; //FIXME: how are shared pointers compared??
			case Empty: // empty expressions are always equal
				return false;
		}
	}

	assert(false);
	return false;
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
			case Star: // reduce to element comparison, fall through case
			case Plus:
				return *this->rexp == *expr.rexp;
			case Empty: // empty expressions are always equal
				return true;
		}
	}
	assert(false); // should not happen
	return false;
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
	else if(this->type == Plus)
		ss << *this->rexp << "¤"; // TODO: fix output...
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
