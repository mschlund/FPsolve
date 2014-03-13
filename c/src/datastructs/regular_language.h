#pragma once

#include <string>
#include "finite_automaton.h"


template<typename RL>
class RegularLanguage {
public:
    virtual ~RegularLanguage(){};
	bool empty();
	RegularLanguage<RL> minimize();
	RegularLanguage<RL> complement();
	RegularLanguage<RL> intersectionWith(RegularLanguage<RL> other);
	RegularLanguage<RL> unionWith(RegularLanguage<RL> other);
    RegularLanguage<RL> concatenate(RegularLanguage<RL> other);
	RegularLanguage<RL> kleeneStar();

	std::string string();

	bool disjoint(RegularLanguage<RL> other);
	bool containedIn(RegularLanguage<RL> super);
	bool contains(RegularLanguage<RL> sub);
	bool equals(RegularLanguage<RL> other);
};
