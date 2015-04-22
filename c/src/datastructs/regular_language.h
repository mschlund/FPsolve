#pragma once

#include <string>

template<typename RL>
class RegularLanguage {
public:
    virtual ~RegularLanguage(){};
	bool empty() const;
	RegularLanguage<RL> minimize() const;
	RegularLanguage<RL> complement() const;
	RegularLanguage<RL> intersectionWith(RegularLanguage<RL> other) const;
	RegularLanguage<RL> unionWith(RegularLanguage<RL> other) const;
    RegularLanguage<RL> concatenate(RegularLanguage<RL> other) const;
	RegularLanguage<RL> kleeneStar() const;

	std::string string() const;

	bool disjoint(RegularLanguage<RL> other) const;
	bool contains(RegularLanguage<RL> sub) const;
	bool equals(RegularLanguage<RL> other) const;
};
