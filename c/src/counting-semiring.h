/*
 * counting-semiring.h
 *
 *  Created on: 30.08.2012
 *      Author: schlund
 */

#ifndef COUNTING_SEMIRING_H_
#define COUNTING_SEMIRING_H_

class CountingSemiring : public Semiring<CountingSemiring>
{

public:
	static bool is_idempotent = true;
	static bool is_commutative = true;
};



#endif /* COUNTING_SEMIRING_H_ */
