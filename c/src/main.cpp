#include <iostream>
#include "float-semiring.h"
#include "free-semiring.h"
#include "matrix.h"
#include "polynomial.h"
#include "newton.h"

void test_addition()
{
	FloatSemiring first(0.5);
	FloatSemiring second(0.3);
	FloatSemiring result = first + second;
	std::cout << "- float semiring addition: " << std::endl;
	std::cout << first << " + " << second << " = " << result << std::endl;
}

void test_multiplication()
{
	FloatSemiring first(0.5);
	FloatSemiring second(0.3);
	FloatSemiring result = first * second;
	std::cout << "- float semiring multiplication: " << std::endl;
	std::cout << first << " * " << second << " = " << result << std::endl;
}

void test_variables()
{
	Var var1;
	Var var2;
	std::cout << "- variables:" << std::endl;
	std::cout << var1 << "," << var2 << std::endl;
}

void test_monomials()
{
	std::cout << "- monomials:" << std::endl;
	Monomial<FloatSemiring> mon1(FloatSemiring(5),{Var("x"),Var("x"),Var("y")});
	Monomial<FloatSemiring> mon2(FloatSemiring(2),{Var("y"), Var("x")});
	Monomial<FloatSemiring> res1 = mon1 + mon1;
	Monomial<FloatSemiring> res2 = mon1 * mon2;
	std::cout << mon1 << " + " << mon1 << " = " << res1 << std::endl;
	std::cout << mon1 << " * " << mon2 << " = " << res2 << std::endl;
	std::cout << "derivative of " << mon1 << " = " << mon1.derivative(Var("x")) << std::endl;
}

Polynomial<FloatSemiring> get_first_polynomial()
{
	// define new polynomial 2xx+7y
	return Polynomial<FloatSemiring>(
			{	Monomial<FloatSemiring>(FloatSemiring(2),{Var("x"),Var("x")}),
				Monomial<FloatSemiring>(FloatSemiring(7),{Var("z")}) });
}

Polynomial<FloatSemiring> get_second_polynomial()
{
	// define new polynomial 5xx+3xy+6yy
	return Polynomial<FloatSemiring>(
			{	Monomial<FloatSemiring>(FloatSemiring(5),{Var("x"),Var("x")}),
				Monomial<FloatSemiring>(FloatSemiring(3),{Var("x"),Var("y")}),
				Monomial<FloatSemiring>(FloatSemiring(6),{Var("y"),Var("y")})});
}

void test_polynomial_addition()
{
	Polynomial<FloatSemiring> first = get_first_polynomial();
	Polynomial<FloatSemiring> second = get_second_polynomial();
	Polynomial<FloatSemiring> result = first + second;
	std::cout << "- polynomial addition:" << std::endl;
	std::cout << first << " + " << second << " = " << result << std::endl;
}

void test_polynomial_multiplication()
{
	Polynomial<FloatSemiring> first = get_first_polynomial();
	Polynomial<FloatSemiring> second = get_second_polynomial();
	Polynomial<FloatSemiring> result = first * second;
	std::cout << "- polynomial multiplication:" << std::endl;
	std::cout << first << " * " << second << " = " << result << std::endl;
}

std::map<Var,FloatSemiring> get_variable_values()
{
	std::map<Var,FloatSemiring> values;
	values.insert(values.begin(), std::pair<Var,FloatSemiring>(Var("x"),FloatSemiring(5)));
	values.insert(values.begin(), std::pair<Var,FloatSemiring>(Var("y"),FloatSemiring(4)));
	values.insert(values.begin(), std::pair<Var,FloatSemiring>(Var("z"),FloatSemiring(0)));
	return values;
}

void test_polynomial_evaluation()
{
	Polynomial<FloatSemiring> polynomial = get_second_polynomial();
	std::map<Var,FloatSemiring> values = get_variable_values();
	FloatSemiring result = polynomial.eval(values);
	std::cout << "- polynomial evaluation:" << std::endl;
	std::cout << polynomial << " at {x=5, y=4, z=0} = " << result << std::endl;
}

Matrix<FloatSemiring> get_first_matrix()
{
	FloatSemiring elems[] = {
		FloatSemiring(0.5), FloatSemiring(0.9),
		FloatSemiring(0.2), FloatSemiring(0.1)};
	return Matrix<FloatSemiring>(2,2,std::vector<FloatSemiring>(elems, elems+4));
}

Matrix<FloatSemiring> get_second_matrix()
{
	FloatSemiring elems[] = {
		FloatSemiring(0.3), FloatSemiring(0.1),
		FloatSemiring(0.7), FloatSemiring(0.2)};
	return Matrix<FloatSemiring>(2,2,std::vector<FloatSemiring>(elems, elems+4));
}

Matrix<FloatSemiring> get_star_test_matrix()
{
	FloatSemiring elems[] = {
		FloatSemiring(0.6), FloatSemiring(0.1), FloatSemiring(0.1),
		FloatSemiring(0.1), FloatSemiring(0.4), FloatSemiring(0.1),
		FloatSemiring(0.1), FloatSemiring(0.1), FloatSemiring(0.3)};
	return Matrix<FloatSemiring>(3,3,std::vector<FloatSemiring>(elems, elems+9));
}
void test_matrix_addition()
{
	Matrix<FloatSemiring> first = get_first_matrix();
	Matrix<FloatSemiring> second = get_second_matrix();
	Matrix<FloatSemiring> result = first + second;
	std::cout << "- matrix addition:" << std::endl;
	std::cout << first << " + " << std::endl << second << " = " << std::endl << result;
}

void test_matrix_multiplication()
{
	Matrix<FloatSemiring> first = get_first_matrix();
	Matrix<FloatSemiring> second = get_second_matrix();
	Matrix<FloatSemiring> result = first * second;
	std::cout << "- matrix multiplication:" << std::endl;
	std::cout << first << " * " << std::endl << second << " = " << std::endl << result;
}

void test_matrix_transpose()
{
	Matrix<FloatSemiring> matrix = get_first_matrix();
	Matrix<FloatSemiring> result = matrix.transpose();
	std::cout << "- matrix transpose:" << std::endl;
	std::cout << matrix << " transposed " << std::endl << result;
}

void test_matrix_star()
{
	Matrix<FloatSemiring> matrix = get_star_test_matrix();
	Matrix<FloatSemiring> result = matrix.star();
	std::cout << "- matrix star:" << std::endl;
	std::cout << matrix << " the matrix star of this is: " << std::endl << result;
}

void test_jacobian()
{
	std::vector<Var> variables;
	variables.push_back(Var("x"));
	variables.push_back(Var("y"));
	variables.push_back(Var("z"));
	std::vector<Polynomial<FloatSemiring> > polynomials;
	polynomials.push_back(get_first_polynomial());
	polynomials.push_back(get_second_polynomial());
	Matrix<Polynomial<FloatSemiring> > jacobian = Polynomial<FloatSemiring>::jacobian(polynomials, variables);
	std::cout << "- jacobian (x,y,z) of " << std::endl;
	std::cout << polynomials.at(0) << std::endl << polynomials.at(1) << std::endl << " = " << std::endl << jacobian << std::endl;
}

void test_polynomial_matrix_evaluation()
{
	std::vector<Var> variables;
	variables.push_back(Var("x"));
	variables.push_back(Var("y"));
	variables.push_back(Var("z"));
	std::vector<Polynomial<FloatSemiring> > polynomials;
	polynomials.push_back(get_first_polynomial());
	polynomials.push_back(get_second_polynomial());
	Matrix<Polynomial<FloatSemiring> > polynomial_matrix = Polynomial<FloatSemiring>::jacobian(polynomials, variables);
	std::map<Var,FloatSemiring> values = get_variable_values();
	Matrix<FloatSemiring> result = Polynomial<FloatSemiring>::eval(polynomial_matrix, values);
	std::cout << "- polynomial matrix evaluation" << std::endl;
	std::cout << polynomial_matrix << " at {x=5, y=4, z=0} = " << std::endl << result;
}

void test_freesemiring()
{
	std::cout << "- free semiring" << std::endl;
	FreeSemiring term1 = FreeSemiring(Var("a"));
	FreeSemiring term2 = FreeSemiring(Var("b"));
	std::cout << "some terms: " <<  term1 << " and " << term2 << std::endl;
	FreeSemiring termAdd = term1 + term2;
	FreeSemiring termMul = term1 * term2;
	std::cout << "addition&multiplication&star: " << (termAdd*termMul.star()) << std::endl;
}

void test_polynomial_to_freesemiring()
{
	std::cout << "- polynomial to free semiring conversion" <<std::endl;
	Polynomial<FloatSemiring> poly = get_first_polynomial();
	auto valuation = new std::unordered_map<FloatSemiring, FreeSemiring, FloatSemiring>();
	FreeSemiring elem = poly.make_free(valuation);
	for(auto v_it = valuation->begin(); v_it != valuation->end(); ++v_it)
	{
		std::cout << "valuation: " << v_it->first << " → " << v_it->second << std::endl;
	}
	std::cout << "polynomial: " << poly << std::endl;
	std::cout << "free semiring: " << elem << std::endl;
}


std::vector<Polynomial<FloatSemiring> > get_newton_test_polynomials()
{
	std::vector<Polynomial<FloatSemiring> > polynomials;

	// define new polynomial 0.4xy+0.6
	Polynomial<FloatSemiring> f1 = Polynomial<FloatSemiring>({
		Monomial<FloatSemiring>(FloatSemiring(0.4), {Var("x"),Var("y")}),
		Monomial<FloatSemiring>(FloatSemiring(0.6), {}) });

	// define new polynomial 0.3yz+0.4yx+0.3
	Polynomial<FloatSemiring> f2 = Polynomial<FloatSemiring>({
		Monomial<FloatSemiring>(FloatSemiring(0.3), {Var("y"),Var("z")}),
		Monomial<FloatSemiring>(FloatSemiring(0.4), {Var("y"),Var("x")}),
		Monomial<FloatSemiring>(FloatSemiring(0.3), {}) });

	// define new polynomial 0.3x+0.7
	Polynomial<FloatSemiring> f3 =  Polynomial<FloatSemiring>({
		Monomial<FloatSemiring>(FloatSemiring(0.3), {Var("x")}),
		Monomial<FloatSemiring>(FloatSemiring(0.7), {}) });
	polynomials.push_back(f1);
	polynomials.push_back(f2);
	polynomials.push_back(f3);

	return polynomials;
}

void test_newton()
{
	Newton<FloatSemiring> newton;
	std::vector<Var> variables;
	variables.push_back(Var("x"));
	variables.push_back(Var("y"));
	variables.push_back(Var("z"));
	std::cout << "- newton:" << std::endl;
	std::vector<Polynomial<FloatSemiring> > polynomials = get_newton_test_polynomials();
	Matrix<FloatSemiring> result = newton.solve_fixpoint(polynomials, variables, 10);
	std::cout << result << std::endl;
}


int main(int argc, char* argv[])
{
	std::cout << "testing..." << std::endl;
	//test_addition();
	//test_multiplication();
	//test_variables();
	//test_monomials();
	//test_polynomial_addition();
	//test_polynomial_multiplication();
	//test_polynomial_evaluation();
	//test_matrix_addition();
	//test_matrix_multiplication();
	//test_matrix_transpose();
	//test_matrix_star();
	//test_jacobian();
	//test_polynomial_matrix_evaluation();
	//test_freesemiring();
	test_polynomial_to_freesemiring();
	test_newton();

	return 0;
}
