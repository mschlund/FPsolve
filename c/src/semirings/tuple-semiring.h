#ifndef TUPLE_SEMIRING_H
#define TUPLE_SEMIRING_H

#include <string>
#include <tuple>
#include <memory>
#include <cassert>
#include <sstream>
#include "semiring.h"


// is it commutative or idempotent?
template <typename SR_A, typename SR_B>
class TupleSemiring : public Semiring<TupleSemiring<SR_A,SR_B>, Commutativity::NonCommutative, Idempotence::Idempotent>
{
private:
	std::tuple<SR_A, SR_B> val;
	static std::shared_ptr<TupleSemiring<SR_A, SR_B>> elem_null;
	static std::shared_ptr<TupleSemiring<SR_A, SR_B>> elem_one;
public:
	TupleSemiring();
	TupleSemiring(SR_A fst, SR_B snd);
	virtual ~TupleSemiring();
	TupleSemiring<SR_A, SR_B> operator += (const TupleSemiring<SR_A, SR_B>& elem);
	TupleSemiring<SR_A, SR_B> operator *= (const TupleSemiring<SR_A, SR_B>& elem);
	bool operator == (const TupleSemiring& elem) const;
	TupleSemiring<SR_A, SR_B> star () const;
	static TupleSemiring<SR_A, SR_B> null();
	static TupleSemiring<SR_A, SR_B> one();
	std::string string() const;
};

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B>::TupleSemiring()
{
  this->val = std::tuple<SR_A,SR_B>{SR_A::null(), SR_B::null()};
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B>::TupleSemiring(SR_A fst, SR_B snd)
{
  this->val = std::tuple<SR_A,SR_B>{fst, snd};
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B>::~TupleSemiring()
{
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B> TupleSemiring<SR_A,SR_B>::operator+=(const TupleSemiring<SR_A,SR_B>& elem)
{
  std::get<0>(this->val) += std::get<0>(elem.val);
  std::get<1>(this->val) += std::get<1>(elem.val);
  return *this;
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B> TupleSemiring<SR_A,SR_B>::operator*=(const TupleSemiring<SR_A,SR_B>& elem)
{
  std::get<0>(this->val) *= std::get<0>(elem.val);
  std::get<1>(this->val) *= std::get<1>(elem.val);
  return *this;
}

template <typename SR_A, typename SR_B>
bool TupleSemiring<SR_A,SR_B>::operator==(const TupleSemiring<SR_A,SR_B>& elem) const
{
  return
    (std::get<0>(this->val) == std::get<0>(elem.val) )&&
    (std::get<1>(this->val) == std::get<1>(elem.val) );
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B> TupleSemiring<SR_A,SR_B>::star() const
{
  TupleSemiring<SR_A,SR_B> result = {std::get<0>(this->val).star(), std::get<1>(this->val).star()};
  return result;
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B> TupleSemiring<SR_A,SR_B>::null()
{
  if(!TupleSemiring<SR_A,SR_B>::elem_null)
    TupleSemiring<SR_A,SR_B>::elem_null = std::shared_ptr<TupleSemiring<SR_A,SR_B>>(new TupleSemiring<SR_A,SR_B>());
  return *TupleSemiring<SR_A,SR_B>::elem_null;
}

template <typename SR_A, typename SR_B>
TupleSemiring<SR_A,SR_B> TupleSemiring<SR_A,SR_B>::one()
{
  if(!TupleSemiring<SR_A,SR_B>::elem_one)
    TupleSemiring<SR_A,SR_B>::elem_one = std::shared_ptr<TupleSemiring<SR_A,SR_B>>(new TupleSemiring<SR_A,SR_B>(SR_A::one(), SR_B::one()));
  return *TupleSemiring<SR_A,SR_B>::elem_one;
}

template <typename SR_A, typename SR_B>
std::string TupleSemiring<SR_A,SR_B>::string() const
{
  std::stringstream ss;
  ss << "<";
  ss << std::get<0>(this->val).string() << "," << std::get<1>(this->val);
  ss << ">";

  return ss.str();
}

template <typename SR_A, typename SR_B>
std::shared_ptr<TupleSemiring<SR_A,SR_B>> TupleSemiring<SR_A,SR_B>::elem_null;
template <typename SR_A, typename SR_B>
std::shared_ptr<TupleSemiring<SR_A,SR_B>> TupleSemiring<SR_A,SR_B>::elem_one;

#endif
