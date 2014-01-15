#pragma once

#include <string>
#include "libfa_fa.h"

enum representation {libfa};
class RegularLanguage {
public:
	bool empty();
	RegularLanguage minimize();
	RegularLanguage complement();
	RegularLanguage intersectionWith(RegularLanguage other);
	RegularLanguage unionWith(RegularLanguage other);
	std::string regularExpression();

	RegularLanguage toType(representation type) {
		std::string expression = regularExpression();

		switch(type) {
		    case libfa: return LibfaFA(expression);
		}
	}

	representation getType();

	static bool disjoint(RegularLanguage A, RegularLanguage B) {
		return A.intersectionWith(B).empty();
	}

	static bool isSubsetOf(RegularLanguage sub, RegularLanguage super) {
	    return super.complement().intersectionWith(sub).empty();
	}
};
