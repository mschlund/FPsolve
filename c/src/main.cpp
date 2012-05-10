#include <iostream>
#include "float-semiring.h"
#include "matrix.h"

int main(int argc, char* argv[])
{
	std::cout << "hello world" << std::endl;
	FloatSemiring foo(0.5);
	FloatSemiring bar(0.2);
	FloatSemiring baz = foo * bar;
	std::cout << foo << " + " << bar << " = " << baz << std::endl;

	Matrix<FloatSemiring> bazinga(3,3,FloatSemiring(0.6));
	bazinga = bazinga + bazinga;
	std::cout << bazinga << std::endl;

	FloatSemiring elems[] = {
		FloatSemiring(0.5), FloatSemiring(0.9),
		FloatSemiring(0.2), FloatSemiring(0.1)};
	Matrix<FloatSemiring> foobar(2,2,std::vector<FloatSemiring>(elems, elems+4));
	foobar = foobar * foobar;
	std::cout << foobar << std::endl;
	return 0;
}
