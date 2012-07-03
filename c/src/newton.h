#ifndef NEWTON_H
#define NEWTON_H

#include <cstdio>
#include "matrix.h"
#include "polynomial.h"

template <typename SR>
class Newton
{
private:
	Matrix<Polynomial<SR> > compute_symbolic_delta(const Matrix<Polynomial<SR> >& u, const std::vector<Polynomial<SR> >& F)
	{
		std::vector<Polynomial<SR> > ret;
		for(int i = 0; i<F.size(); i++)
		{
			Matrix<Polynomial<SR> > H = F.at(i).hessian();
			// FIXME: 0.5 cannot be represented universally
			Matrix<Polynomial<SR> > d = u.transpose() * H * u;
			ret.push_back(d.getElements().at(0)); // TODO: beautify this
		}

		return Matrix<Polynomial<SR> >(1, F.size(), ret);
	}
public:
	// calculate the next newton iterand
	Matrix<Polynomial<SR> > step(const std::vector<char>& poly_vars, const Matrix<Polynomial<SR> >& J_s,
			const Matrix<Polynomial<SR> >& v, const Matrix<Polynomial<SR> >& delta)
	{
		assert(poly_vars.size() == v.getRows());
		std::map<char,Polynomial<SR> > values;
		for(int i = 0; i<poly_vars.size(); i++)
		{
			values.insert(values.begin(), std::pair<char,Polynomial<SR> >(poly_vars.at(i),v.getElements().at(i)));
		}
		Matrix<Polynomial<SR> > J_s_new = Polynomial<SR>::eval(J_s, values);
		Matrix<Polynomial<SR> > result = J_s_new * delta;
		return result;
	}

	// iterate until convergence
	Matrix<Polynomial<SR> > solve_fixpoint(const std::vector<Polynomial<SR> >& F, const std::vector<char>& poly_vars, int max_iter)
	{
		Matrix<Polynomial<SR> > J = Polynomial<SR>::jacobian(F, poly_vars);
		Matrix<Polynomial<SR> > J_s = J.star();
		// define new symbolic vector [u1,u2,...,un] TODO: this is ugly...
		std::vector<Polynomial<SR> > ret;
		std::map<std::string, FloatSemiring> coeff;
		std::set<char> variables;
		char var[] = "0"; // initialize char*
		for(int i=0; i<poly_vars.size(); i++)
		{
			coeff.clear();
			variables.clear();
			sprintf(var,"%d",i);
			coeff[var] = SR::one(); // variable "1", "2",...
			variables.insert(var[0]);
			Polynomial<SR> f = Polynomial<SR>(variables, coeff);
			ret.push_back(f);
		}
		Matrix<Polynomial<SR> > u = Matrix<Polynomial<SR> >(1, F.size(), ret);
		Matrix<Polynomial<SR> > delta = compute_symbolic_delta(u, F);

		Matrix<Polynomial<SR> > v = Matrix<Polynomial<SR> >(1,(int)F.size()); // v^0 = 0

		Matrix<Polynomial<SR> > v_upd = step(poly_vars, J_s, v, delta);

		std::cout << v_upd;
		return v_upd;
	}
};

#endif
