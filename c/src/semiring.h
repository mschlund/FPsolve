#ifndef SEMIRING_H
#define SEMIRING_H
#include <iostream>
#include <string>
#include <functional> // for std::hash


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
	size_t operator()(const SR& sr) const
	{
		return std::hash<std::string>()(sr.string());
	}
};

template <typename SR>
std::ostream& operator<<(std::ostream& os, const Semiring<SR>& elem)
{
	return os << elem.string();
}

#endif
