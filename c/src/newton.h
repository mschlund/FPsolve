#ifndef NEWTON_H
#define NEWTON_H

#include <cstdio>
#include "matrix.h"
#include "polynomial.h"

template <typename SR>
class Newton
{
private:
	Matrix<Polynomial<SR> > compute_symbolic_delta(Matrix<Polynomial<SR> > u, std::vector<Polynomial<SR> > F)
	{
		std::vector<Polynomial<SR> > ret;
		for(int i = 0; i<F.size(); i++)
		{
			Matrix<Polynomial<SR> > H = F.at(i).hessian();
			// FIXME: 0.5 cannot be represented universally
			ret.push_back(0.5*u.transpose()*H*u);
		}

		return Matrix<Polynomial<SR> >(1, F.size(), ret);
	}
public:
	//Matrix<Polynomial<SR> > step(Matrix<Polynomial<SR> > F, std::set<char> poly_vars, Matrix<Polynomial<SR> > J_s,
	//		Matrix<Polynomial<SR> > v, Matrix<Polynomial<SR> > delta);
	Matrix<Polynomial<SR> > solve_fixpoint(std::vector<Polynomial<SR> > F, std::vector<char> poly_vars, int max_iter)
	{
		Matrix<Polynomial<SR> > J = Polynomial<SR>::jacobian(F, poly_vars);
		Matrix<Polynomial<SR> > J_s = J.star();
		// define new symbolic vector [u1,u2,u3] TODO: this is ugly...
		std::vector<Polynomial<SR> > ret;
		std::map<std::string, FloatSemiring> coeff;
		std::set<char> variables;
		char var[] = "0"; // initialize char*
		for(int i=0; i<F.size(); i++)
		{
			coeff.clear();
			sprintf(var,"%d",i);
			coeff[var] = SR::one(); // variable "1", "2",...
			variables.insert(var[0]);
			Polynomial<SR> f = Polynomial<SR>(variables, coeff);
			ret.push_back(f);
		}
		Matrix<Polynomial<SR> > u = Matrix<Polynomial<SR> >(1, F.size(), ret);
		Matrix<Polynomial<SR> > delta = compute_symbolic_delta(u, F);

		return delta;
	}
};

#endif
