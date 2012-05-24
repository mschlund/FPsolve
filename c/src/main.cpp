#include <iostream>
#include "float-semiring.h"
#include "matrix.h"
#include "polynomial.h"

int main(int argc, char* argv[])
{
	FloatSemiring foo(0.5);
	FloatSemiring bar(0.2);

	Matrix<FloatSemiring> bazinga(3,3,FloatSemiring(0.6));
	bazinga = bazinga + bazinga;

	FloatSemiring elems[] = {
		FloatSemiring(0.5), FloatSemiring(0.9),
		FloatSemiring(0.2), FloatSemiring(0.1)};
	Matrix<FloatSemiring> foobar(2,2,std::vector<FloatSemiring>(elems, elems+4));
	foobar = foobar * foobar;

	// define new polynomial 2xx+7y
	std::map<std::string, FloatSemiring> coeff;
	coeff["xx"] = FloatSemiring(2);
	coeff["z"] = FloatSemiring(7);
	std::set<char> variables;
	variables.insert('x');
	variables.insert('z');
	Polynomial<FloatSemiring> poly1(variables,coeff);

	// define new polynomial 5xx+3xy+6yy
	coeff.clear();
	coeff["xx"] = FloatSemiring(5);
	coeff["xy"] = FloatSemiring(3);
	coeff["yy"] = FloatSemiring(6);
	variables.clear();
	variables.insert('x');
	variables.insert('y');
	Polynomial<FloatSemiring> poly2(variables,coeff);

	Polynomial<FloatSemiring> tmpPoly;
	tmpPoly = poly1 + poly2;
	std::cout << "(" << poly1 << ") + (" << poly2 << ") = " << tmpPoly << std::endl;
	tmpPoly = poly2 * poly2;
	std::cout << "(" << poly2 << ")^2 = " << tmpPoly << std::endl;


	return 0;
}
