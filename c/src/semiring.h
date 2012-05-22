#ifndef SEMIRING_H
#define SEMIRING_H
#include <iostream>
#include <string>


template <typename SR>
class Semiring {
private:
protected:
	Semiring() {};
public:
	virtual SR operator * (const SR& elem) const = 0;
	virtual SR operator + (const SR& elem) const = 0;
	virtual SR star () const = 0;
	virtual bool is_idempotent() const = 0;
	virtual bool is_commutative() const = 0;
	SR null() { return SR::null(); }; // TODO: make this pure virtual
	virtual std::string string() const = 0;
};

template <typename SR>
std::ostream& operator<<(std::ostream& os, Semiring<SR>& elem)
{
	return os << elem.string();
}

#endif
