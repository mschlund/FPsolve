#ifndef SEMIRING_H
#define SEMIRING_H
#include <iostream>
#include <string>


template <typename SR>
class Semiring {
public:
	virtual SR operator * (const SR& elem) const = 0;
	virtual SR operator + (const SR& elem) const = 0;
	virtual SR star () const = 0;
	static bool is_idempotent;
	static bool is_commutative;
	static SR null();
	static SR one();
	virtual std::string string() const = 0;
};

template <typename SR>
std::ostream& operator<<(std::ostream& os, const Semiring<SR>& elem)
{
	return os << elem.string();
}

#endif
