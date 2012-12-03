#include <iostream>
#include <string>
#include <cstdlib>
#include "float-semiring.h"
#include "free-semiring.h"
#include "semilinSetExp.h"
#include "matrix.h"
#include "polynomial.h"
#include "newton.h"
#include "commutativeRExp.h"
#include "parser.h"

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
	VarPtr var1 = Var::getVar();
	VarPtr var2 = Var::getVar();
	std::cout << "- variables:" << std::endl;
	std::cout << var1 << "," << var2 << std::endl;
}

void test_monomials()
{
	std::cout << "- monomials:" << std::endl;
	Monomial<FloatSemiring> mon1(FloatSemiring(5),{Var::getVar("x"),Var::getVar("x"),Var::getVar("y")});
	Monomial<FloatSemiring> mon2(FloatSemiring(2),{Var::getVar("y"), Var::getVar("x")});
	Monomial<FloatSemiring> res1 = mon1 + mon1;
	Monomial<FloatSemiring> res2 = mon1 * mon2;
	std::cout << mon1 << " + " << mon1 << " = " << res1 << std::endl;
	std::cout << mon1 << " * " << mon2 << " = " << res2 << std::endl;
	std::cout << "derivative of " << mon1 << " = " << mon1.derivative(Var::getVar("x")) << std::endl;
}

Polynomial<FloatSemiring> get_first_polynomial()
{
	// define new polynomial 2xx+7y
	return Polynomial<FloatSemiring>(
			{	Monomial<FloatSemiring>(FloatSemiring(2),{Var::getVar("x"),Var::getVar("x")}),
				Monomial<FloatSemiring>(FloatSemiring(7),{Var::getVar("z")}) });
}

Polynomial<FloatSemiring> get_second_polynomial()
{
	// define new polynomial 5xx+3xy+6yy
	return Polynomial<FloatSemiring>(
			{	Monomial<FloatSemiring>(FloatSemiring(5),{Var::getVar("x"),Var::getVar("x")}),
				Monomial<FloatSemiring>(FloatSemiring(3),{Var::getVar("x"),Var::getVar("y")}),
				Monomial<FloatSemiring>(FloatSemiring(6),{Var::getVar("y"),Var::getVar("y")})});
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

std::map<VarPtr,FloatSemiring> get_variable_values()
{
	std::map<VarPtr,FloatSemiring> values;
	values.insert(values.begin(), std::pair<VarPtr,FloatSemiring>(Var::getVar("x"),FloatSemiring(5)));
	values.insert(values.begin(), std::pair<VarPtr,FloatSemiring>(Var::getVar("y"),FloatSemiring(4)));
	values.insert(values.begin(), std::pair<VarPtr,FloatSemiring>(Var::getVar("z"),FloatSemiring(0)));
	return values;
}

void test_polynomial_evaluation()
{
	Polynomial<FloatSemiring> polynomial = get_second_polynomial();
	std::map<VarPtr,FloatSemiring> values = get_variable_values();
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

Matrix<FreeSemiring> get_star_test_matrix()
{
/*	FloatSemiring elems[] = {
		FloatSemiring(0.6), FloatSemiring(0.1), FloatSemiring(0.1),
		FloatSemiring(0.1), FloatSemiring(0.4), FloatSemiring(0.1),
		FloatSemiring(0.1), FloatSemiring(0.1), FloatSemiring(0.3)};
	return Matrix<FloatSemiring>(3,3,std::vector<FloatSemiring>(elems, elems+9));
*/
	FreeSemiring elems[]  = {
		FreeSemiring(Var::getVar("a")),
		FreeSemiring(Var::getVar("b")),
		FreeSemiring(Var::getVar("c")),
		FreeSemiring(Var::getVar("d")) };
	return Matrix<FreeSemiring>(2,2,std::vector<FreeSemiring>(elems, elems+4));
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
	auto matrix = get_star_test_matrix();
	auto result = matrix.star();
	std::cout << "- matrix star:" << std::endl;
	std::cout << matrix << " the matrix star of this is: " << std::endl << result;
}

void test_jacobian()
{
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	variables.push_back(Var::getVar("y"));
	variables.push_back(Var::getVar("z"));
	std::vector<Polynomial<FloatSemiring> > polynomials;
	polynomials.push_back(get_first_polynomial());
	polynomials.push_back(get_second_polynomial());
	Matrix<Polynomial<FloatSemiring> > jacobian = Polynomial<FloatSemiring>::jacobian(polynomials, variables);
	std::cout << "- jacobian (x,y,z) of " << std::endl;
	std::cout << polynomials.at(0) << std::endl << polynomials.at(1) << std::endl << " = " << std::endl << jacobian << std::endl;
}

void test_polynomial_matrix_evaluation()
{
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	variables.push_back(Var::getVar("y"));
	variables.push_back(Var::getVar("z"));
	std::vector<Polynomial<FloatSemiring> > polynomials;
	polynomials.push_back(get_first_polynomial());
	polynomials.push_back(get_second_polynomial());
	Matrix<Polynomial<FloatSemiring> > polynomial_matrix = Polynomial<FloatSemiring>::jacobian(polynomials, variables);
	std::map<VarPtr,FloatSemiring> values = get_variable_values();
	Matrix<FloatSemiring> result = Polynomial<FloatSemiring>::eval(polynomial_matrix, values);
	std::cout << "- polynomial matrix evaluation" << std::endl;
	std::cout << polynomial_matrix << " at {x=5, y=4, z=0} = " << std::endl << result;
}

void test_freesemiring()
{
	std::cout << "- free semiring" << std::endl;
	FreeSemiring term1 = FreeSemiring(Var::getVar("a"));
	FreeSemiring term2 = FreeSemiring(Var::getVar("b"));
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
		Monomial<FloatSemiring>(FloatSemiring(0.4), {Var::getVar("x"),Var::getVar("y")}),
		Monomial<FloatSemiring>(FloatSemiring(0.6), {}) });

	// define new polynomial 0.3yz+0.4yx+0.3
	Polynomial<FloatSemiring> f2 = Polynomial<FloatSemiring>({
		Monomial<FloatSemiring>(FloatSemiring(0.3), {Var::getVar("y"),Var::getVar("z")}),
		Monomial<FloatSemiring>(FloatSemiring(0.4), {Var::getVar("y"),Var::getVar("x")}),
		Monomial<FloatSemiring>(FloatSemiring(0.3), {}) });

	// define new polynomial 0.3x+0.7
	Polynomial<FloatSemiring> f3 =  Polynomial<FloatSemiring>({
		Monomial<FloatSemiring>(FloatSemiring(0.3), {Var::getVar("x")}),
		Monomial<FloatSemiring>(FloatSemiring(0.7), {}) });
	polynomials.push_back(f1);
	polynomials.push_back(f2);
	polynomials.push_back(f3);

	return polynomials;
}

void test_newton()
{
	/*
	Newton<FloatSemiring> newton;
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	variables.push_back(Var::getVar("y"));
	variables.push_back(Var::getVar("z"));
	std::cout << "- newton (float):" << std::endl;
	std::vector<Polynomial<FloatSemiring> > polynomials = get_newton_test_polynomials();
	Matrix<FloatSemiring> result = newton.solve_fixpoint(polynomials, variables, 2);
	std::cout << result << std::endl;

*/

/*
	Newton<SemilinSetExp> newton;
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	std::cout << "- newton (cnt-SR):" << std::endl;

	std::vector<Polynomial<SemilinSetExp> > polynomials;
	Polynomial<SemilinSetExp> f1 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("a")), {Var::getVar("x"),Var::getVar("x")}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("c")), {}) });

	polynomials.push_back(f1);

	Matrix<SemilinSetExp> result = newton.solve_fixpoint(polynomials, variables, 6);
	std::cout << result << std::endl;
*/

	Newton<CommutativeRExp> newton;
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	variables.push_back(Var::getVar("y"));
	variables.push_back(Var::getVar("z"));
	std::cout << "- newton (counting-SR):" << std::endl;

	std::vector<Polynomial<CommutativeRExp> > polynomials;
	// define new polynomial axy+b
	Polynomial<CommutativeRExp> f1 = Polynomial<CommutativeRExp>({
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("a")), {Var::getVar("x"),Var::getVar("y")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("b")), {}) });

	// define new polynomial cyz+dyx+e
	Polynomial<CommutativeRExp> f2 = Polynomial<CommutativeRExp>({
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("c")), {Var::getVar("y"),Var::getVar("z")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("d")), {Var::getVar("y"),Var::getVar("x")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("e")), {}) });

	// define new polynomial fx+g
	Polynomial<CommutativeRExp> f3 =  Polynomial<CommutativeRExp>({
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("f")), {Var::getVar("x")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("g")), {}) });
	polynomials.push_back(f1);
	polynomials.push_back(f2);
	polynomials.push_back(f3);

	Matrix<CommutativeRExp> result = newton.solve_fixpoint(polynomials, variables, 2);
	std::cout << result << std::endl;


/*	Newton<SemilinSetExp> newton;
	std::vector<VarPtr> v;
	v.push_back(Var::getVar("x"));
	std::cout << "- newton (cnt-SR):" << std::endl;

	std::vector<Polynomial<SemilinSetExp> > polys;
	Polynomial<SemilinSetExp> f0 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("s")), {Var::getVar("x"),Var::getVar("x")}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("s")), {}) });

	polys.push_back(f0);
	Matrix<SemilinSetExp> r = newton.solve_fixpoint(polys, v, 2);
	std::vector<SemilinSetExp> sol1 = r.getElements();

	std::cout << sol1[0] << std::endl;

//	std::cout << (SemilinSetExp(Var::getVar("r"))*SemilinSetExp(Var::getVar("s"))*sol1[0])<< std::endl;

	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("e"));
	std::vector<Polynomial<SemilinSetExp> > polynomials;
	Polynomial<SemilinSetExp> f1 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("l"))*SemilinSetExp(Var::getVar("s")), {Var::getVar("e")}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("r"))*SemilinSetExp(Var::getVar("s"))*sol1[0], {Var::getVar("e")}) });

	polynomials.push_back(f1);

	Matrix<SemilinSetExp> result = newton.solve_fixpoint(polynomials, variables, 1);
	std::cout << result << std::endl;
	std::cout << (SemilinSetExp(Var::getVar("l"))*SemilinSetExp(Var::getVar("s")) + SemilinSetExp(Var::getVar("r"))*SemilinSetExp(Var::getVar("s"))*sol1[0]).star() << std::endl;
*/

/*
	SemilinSetExp x = SemilinSetExp(Var::getVar("a"));
	SemilinSetExp y = SemilinSetExp(Var::getVar("b"));
	SemilinSetExp z = (x*x*x+y*y*y).star();
	std::cout << z << std::endl;
*/

}

int main(int argc, char* argv[])
{
	//std::cout << "testing..." << std::endl;
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
	//test_polynomial_to_freesemiring();
	if(argc == 1) // no command-line arguments
		test_newton();
	else // there is an argument
	{
		Parser p;
		std::string input;
		if(std::string("-float").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_float(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-free").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_free(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-rexp").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_rexp(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-polyrexp").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_polyrexp(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-grammar").compare(argv[1]) == 0)
		{
			Newton<CommutativeRExp> newton;

			// accumulate the input in rules
			std::vector<Polynomial<CommutativeRExp> > rules;
			// and all the variables in vars
			std::vector<VarPtr> vars;

			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_grammar(input); // (string, polynomial)
				rules.push_back(result.second);
				vars.push_back(Var::getVar(result.first));

				std::cout << result.first << " → " << result.second << std::endl;
			}
			int i = 2; // default value
			if(argc >= 2) // set the given iteration count
				i = std::atoi(argv[2]);

			Matrix<CommutativeRExp> result = newton.solve_fixpoint(rules, vars, i);
			std::cout << i << "-th newton iteration: " << std::endl;
			std::cout << result << std::endl;
		}
		else
		{
			std::cout << "unknown argument" << std::endl;
		}
	}

	return 0;
}
