/*
 * counting-semiring.h
 *
 *  Created on: 30.08.2012
 *      Author: schlund
 */

#ifndef COUNTING_SEMIRING_H_
#define COUNTING_SEMIRING_H_

#include "semirings/semiring.h"

class CountingSemiring : public Semiring<CountingSemiring>
{

public:
	static bool is_idempotent;
	static bool is_commutative;
	virtual std::string string() const = 0;
};

//CountingSemiring::is_idempotent = true;
//ConutingSemiring::is_commutative = true

#endif /* COUNTING_SEMIRING_H_ */
