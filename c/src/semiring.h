#ifndef SEMIRING_H
#define SEMIRING_H

#include <iosfwd>
#include <string>
#include <functional> // for std::hash


template <typename SR>
class Semiring {
public:
  virtual ~Semiring(){};

	friend SR operator * (const SR& lhs, const SR& rhs)
	{
		SR result = lhs;
		result *= rhs;
		return result;
	}
	friend SR operator + (const SR& lhs, const SR& rhs)
	{
		SR result = lhs;
		result += rhs;
		return result;
	}
	virtual SR star () const = 0;
	virtual bool operator ==(const SR& elem) const = 0;
	static bool is_idempotent;
	static bool is_commutative;
	static SR null();
	static SR one();
	virtual std::string string() const = 0;

        // FIXME: This might be _really_ inefficient.  Maybe we should just
        // require that the semirings provide the specialization for
        // std::hash...?
	size_t operator()(const SR& sr) const
	{
		return std::hash<std::string>()(sr.string());
	}
};

template <typename SR>
SR operator *= (SR& lhs, const SR& rhs);
template <typename SR>
SR operator += (SR& lhs, const SR& rhs);

template <typename SR>
std::ostream& operator<<(std::ostream& os, const Semiring<SR>& elem)
{
	return os << elem.string();
}

#endif
