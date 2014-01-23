#pragma once

#include <string>
#include "finite_automaton.h"


template<typename RL>
class RegularLanguage {
public:
	bool empty();
	RegularLanguage<RL> minimize();
	RegularLanguage<RL> complement();
	RegularLanguage<RL> intersectionWith(RegularLanguage<RL> other);
	RegularLanguage<RL> unionWith(RegularLanguage<RL> other);
    RegularLanguage<RL> concatenate(RegularLanguage<RL> other);
	RegularLanguage<RL> kleeneStar();

    /*
     * Returns a regular expression that describes the language of this automaton.
     */
	std::string string();

	bool disjoint(RegularLanguage<RL> other) {
		return intersectionWith(other).empty();
	}

	bool containedIn(RegularLanguage<RL> super) {
	    return intersectionWith(super.complement()).empty();
	}

	bool contains(RegularLanguage<RL> sub) {
	    return complement().intersectionWith(sub).empty();
	}

	bool equals(RegularLanguage<RL> other) {
	    return contains(other) && containedIn(other);
	}
};
