#pragma once

#include <string>

#include "../../libraries/augeas/src/fa.h"

enum representation {libfa};

class regularLanguage {
public:
	boolean empty();
	regularLanguage complement();
	regularLanguage intersectionWith(regularLanguage other);
	regularLanguage unionWith(regularLanguage other);
	regularLanguage minimize();
	std::string regularExpression();
	regularLanguage fromRegularExpression(std::string regularExpression);

	regularLanguage toType(representation type) {
		if(getType() != type) {
			std::string expression = regularExpression();
		}

		return this;
	}

	representation getType(){
		return type;
	}

	static boolean disjoint(regularLanguage A, regularLanguage B) {
		return A.intersectionWith(B).empty();
	}

	static boolean isSubsetOf(regularLanguage sub, regularLanguage super) {

	}

private:
	representation type;
};
