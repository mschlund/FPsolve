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
	virtual SR operator * (const SR& elem) = 0;
	virtual SR operator + (const SR& elem) = 0;
	virtual SR star () = 0;
	virtual bool is_idempotent() = 0;
	virtual bool is_commutative() = 0;
	SR null() { return -1; }; // TODO: make this pure virtual
	virtual std::string string() = 0;
};

template <typename SR>
std::ostream& operator<<(std::ostream& os, Semiring<SR>& elem)
{
	return os << elem.string();
}

#endif
