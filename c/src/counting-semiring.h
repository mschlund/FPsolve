/*
 * counting-semiring.h
 *
 *  Created on: 30.08.2012
 *      Author: schlund
 */

#ifndef COUNTING_SEMIRING_H_
#define COUNTING_SEMIRING_H_

template <typename CountingSRImpl>
class CountingSemiring : public Semiring<CountingSemiring<CountingSRImpl> >
{
private:
	CountingSRImpl val;
	static CountingSemiring<CountingSRImpl>* elem_null = 0;
	static CountingSemiring<CountingSRImpl>* elem_one = 0;

public:
	CountingSemiring<CountingSRImpl>(const CountingSRImpl val) {
		this->val = val;
	}

	void ~CountingSemiring<CountingSRImpl>() {
		if(elem_null) {
			delete(elem_null);
			elem_null =0;
		}
		if(elem_one) {
			delete(elem_one);
			elem_one = 0;
		}
	}

	CountingSemiring<CountingSRImpl> operator + (const CountingSemiring<CountingSRImpl>& elem) {
		return CountingSemiring<CountingSRImpl>(this->val + elem.val);
	}

	CountingSemiring<CountingSRImpl> operator * (const CountingSemiring<CountingSRImpl>& elem) {
		return CountingSemiring<CountingSRImpl>(this->val * elem.val);
	}

	CountingSemiring<CountingSRImpl> star () {
		return CountingSemiring<CountingSRImpl>(this->val.star());
	}

	static CountingSemiring<CountingSRImpl> null() {
		if (!elem_null)
			elem_null = new CountingSRImpl(0);
	}
	static CountingSemiring<CountingSRImpl> one() {
		if (!elem_one)
			elem_one = new CountingSRImpl(0);
	}

	std::string string() const {
		std::stringstream ss;
		ss << this->val;
		return ss.str();
	}

	static bool is_idempotent = true;
	static bool is_commutative = true;
};



#endif /* COUNTING_SEMIRING_H_ */
