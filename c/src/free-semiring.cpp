#include <sstream>
#include "free-semiring.h"

FreeSemiring::FreeSemiring()
{
	// Do not use this constructor. It is needed for a std::unordered_map...

	// dummy type, so valgrind does not complain about uninitialised value on
	// switch-statemenet in copy-constructor
	this->type = Dummy;
}

FreeSemiring::FreeSemiring(VarPtr var)
{
	this->type = Element;
	this->elem = Var::getVar(var);
}

FreeSemiring::FreeSemiring(const FreeSemiring& term)
{
	this->type = term.type;
	switch(this->type)
	{
	case Element:
		this->elem = term.elem;
		break;
	case Addition:
	case Multiplication:
		this->right_ptr = term.right_ptr;
	case Star:
		this->left_ptr = term.left_ptr;
		break;
	case Dummy:
		break; // do not copy anything
	}
}

FreeSemiring::FreeSemiring(optype type, FreeSemiring left, FreeSemiring right)
{
	this->type = type;
	this->left_ptr = std::make_shared<FreeSemiring>(left);
	this->right_ptr = std::make_shared<FreeSemiring>(right);
}

FreeSemiring::FreeSemiring(optype type, FreeSemiring left)
{
	this->type = type;
	this->left_ptr = std::make_shared<FreeSemiring>(left);
}

FreeSemiring::~FreeSemiring()
{
	//if(this->type == Element)
	//	delete this->elem;
}

FreeSemiring FreeSemiring::operator +=(const FreeSemiring& term)
{
	if(*this == FreeSemiring::null())
		*this = term;
	else if(term == FreeSemiring::null())
	; // do not do anything
	else
		*this = FreeSemiring(Addition, *this, term);
	return *this;
}

FreeSemiring FreeSemiring::operator *=(const FreeSemiring& term)
{
	if(*this == FreeSemiring::one())
		*this = term;
	else if(term == FreeSemiring::one())
		; // do not do anything
	else if(*this == FreeSemiring::null())
		; // do not do anything
	else if(term == FreeSemiring::null())
		*this = term;
	else
		*this = FreeSemiring(Multiplication, *this, term);
	return *this;
}

// TODO: at the moment == returns true iff the whole structure is identical
// it does not respect associativity
bool FreeSemiring::operator ==(const FreeSemiring& term) const
{
	switch(this->type)
	{
		case Element:
			if(term.type != Element)
				return false;
			else
				return this->elem == term.elem;
		case Addition:
		case Multiplication:
			return (this->type == term.type && *this->left_ptr == *term.left_ptr && *this->right_ptr == *term.right_ptr);
		case Star:
			return (this->type == term.type && *this->left_ptr == *term.left_ptr);
		case Dummy:
			assert(false);
	}
	return false;
}

FreeSemiring FreeSemiring::star() const
{
	if(*this == FreeSemiring::null())
		return FreeSemiring::one();
	return FreeSemiring(Star, *this);
}

std::string FreeSemiring::string() const
{
	std::stringstream ss;

	switch(this->type)
	{
	case Element:
		ss << this->elem->string();
		break;
	case Addition:
		ss << "(";
		ss << this->left_ptr->string() << " + " << this->right_ptr->string();
		ss << ")";
		break;
	case Multiplication:
		ss << "(";
		ss << this->left_ptr->string() << " . " << this->right_ptr->string();
		ss << ")";
		break;
	case Star:
		//ss << "s(" << this->left_ptr->string() << ")";
		ss << "(";
		ss << this->left_ptr->string() << "*";
		ss << ")";
		break;
	case Dummy:
		assert(false);
	}

	return ss.str();
}

FreeSemiring FreeSemiring::null()
{
	if(!FreeSemiring::elem_null)
		FreeSemiring::elem_null = std::shared_ptr<FreeSemiring>(new FreeSemiring(Var::getVar("Null")));
	return *FreeSemiring::elem_null;
}
FreeSemiring FreeSemiring::one()
{
	if(!FreeSemiring::elem_one)
		FreeSemiring::elem_one = std::shared_ptr<FreeSemiring>(new FreeSemiring(Var::getVar("One")));
	return *FreeSemiring::elem_one;
}

bool FreeSemiring::is_idempotent = false;
bool FreeSemiring::is_commutative = false;
std::shared_ptr<FreeSemiring> FreeSemiring::elem_null;
std::shared_ptr<FreeSemiring> FreeSemiring::elem_one;
