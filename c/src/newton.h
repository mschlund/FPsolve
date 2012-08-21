#ifndef NEWTON_H
#define NEWTON_H

#include <cstdio>
#include "matrix.h"
#include "polynomial.h"

template <typename SR>
class Newton
{
private:
	int sum(int* array, int n)
	{
		int s = 0;
		for(int i=0; i<=n; i++)
		{
			s += array[i];
		}
		return s;
	};

	// generate vectors of vectors [[x0,x1,...xn,]] such that each possible permutation
	// of integers with 0 <= xi <max is in the result. also each permutation has a sum with
	// min_sum <= sum <= max_sum
	std::vector<std::vector<int> > genIdx(int max, int n, int min_sum, int max_sum)
		{
		std::vector<std::vector<int> > result;

		// initialize base array
		int* base = new int[n+1];

		while(base[n] < max)
		{
			int k = 0;
			base[k]++;
			if(base[k]>=max)
			{
				while(k<n)
				{
					base[k] = 0;
					k++;
					base[k]++;
				}
			}
			int s = sum(base, n);
			if(min_sum <= s && s <= max_sum)
			{
				std::vector<int> tmp(base,base+n);
				result.push_back(tmp);
			}
		}

		delete[] base;
		return result;
	}

	Matrix<Polynomial<SR> > compute_symbolic_delta(const std::vector<Var>& v,
							const std::vector<Var>& v_upd,
							const std::vector<Polynomial<SR> >& F,
							const std::vector<Var>& poly_vars)
	{
		std::vector<Polynomial<SR> > delta;
		int n = v.size();

		for(int i=0; i<n; ++i)
		{
			Polynomial<SR> delta_i = Polynomial<SR>::null();
			Polynomial<SR> f = F.at(i);
			int deg = f.get_degree();

			// create [[0,...,0],[0,...,deg],[deg,...,deg]]
			std::vector<std::vector<int> > p = genIdx(deg,n,2,deg);

			//iterate over (x,y,...,z) with x,y,...,z \in [0,deg+1]
			for(std::vector<std::vector<int> >::const_iterator it = p.begin(); it != p.end(); ++it)
			{
				std::vector<Var> dx;
				Polynomial<SR> prod = Polynomial<SR>(SR::one());

				std::vector<Var>::const_iterator var = poly_vars.begin();
				std::vector<Var>::const_iterator elem = v_upd.begin();
				for(std::vector<int>::const_iterator z = it->begin(); z != it->end(); ++z)
				{
					for(int j=0; j<(*z); ++j)
					{
						dx.push_back(*var); // generate a multiset of variables like "xxxyyzzzz";
						prod = prod * (*elem);	// prod = prod * elem^z ... (elem is a variable)
					}
					++var; // next variable
					++elem; // next element of vector v_upd
				}

				// eval f.derivative(dx) at v
				std::map<Var,Var> values;
				int i = 0;
				for(std::vector<Var>::const_iterator poly_var = poly_vars.begin(); poly_var != poly_vars.end(); ++poly_var)
				{
					values.insert(values.begin(), std::pair<Var,Var>((*poly_var),v.at(i++)));
				}
				Polynomial<SR> f_eval = f.derivative(dx).subst(values);

				delta_i = delta_i + f_eval * prod;
			}

			delta.push_back(delta_i);
		}

		return Matrix<Polynomial<SR> >(1, delta.size(), delta);
	}

	std::vector<Var> get_symbolic_vector(int size, std::string prefix)
	{
		// define new symbolic vector [u1,u2,...,un] TODO: this is ugly...
		std::vector<Var> ret;
		for(int i=0; i<size; i++)
		{
			std::stringstream ss;
			ss << prefix << "_" << i;
			ret.push_back(Var(ss.str()));
		}
		return ret;
	}

public:
	// calculate the next newton iterand
	Matrix<SR> step(const std::vector<Var>& poly_vars, const Matrix<Polynomial<SR> >& J_s,
			const Matrix<SR>& v, const Matrix<SR>& delta)
	{
		assert(poly_vars.size() == v.getRows());
		std::map<Var,SR> values;
		int i=0;
		for(std::vector<Var>::const_iterator poly_var = poly_vars.begin(); poly_var != poly_vars.end(); ++poly_var)
		{
			values.insert(values.begin(), std::pair<Var,SR>((*poly_var),v.getElements().at(i++)));
		}
		Matrix<SR> J_s_new = Polynomial<SR>::eval(J_s, values);
		Matrix<SR> result = J_s_new * delta;
		return result;
	}

	// iterate until convergence
	Matrix<SR> solve_fixpoint(const std::vector<Polynomial<SR> >& F, const std::vector<Var>& poly_vars, int max_iter)
	{
		Matrix<Polynomial<SR> > F_mat = Matrix<Polynomial<SR> >(1,F.size(),F);
		Matrix<Polynomial<SR> > J = Polynomial<SR>::jacobian(F, poly_vars);
		Matrix<Polynomial<SR> > J_s = J.star();

		// define new symbolic vectors [u1,u2,...,un] TODO: this is ugly...
		std::vector<Var> u = this->get_symbolic_vector(poly_vars.size(), "u");
		std::vector<Var> u_upd = this->get_symbolic_vector(poly_vars.size(), "u_upd");

		Matrix<Polynomial<SR> > delta = compute_symbolic_delta(u, u_upd, F, poly_vars);
		Matrix<SR> v = Matrix<SR>(1,(int)F.size()); // v^0 = 0

		// d^0 = F(0)
		std::map<Var,SR> values;
		for(std::vector<Var>::const_iterator poly_var = poly_vars.begin(); poly_var != poly_vars.end(); ++poly_var)
		{
			values.insert(values.begin(), std::pair<Var,SR>(*poly_var, SR::null()));
		}
		Matrix<SR> delta_new = Polynomial<SR>::eval(F_mat, values);

		Matrix<SR> v_upd = step(poly_vars, J_s, v, delta_new);

		for(int i=2; i<max_iter; ++i)
		{
			values.clear();
			for(int i = 0; i<u.size(); i++)
			{
				values.insert(values.begin(), std::pair<Var,SR>(u_upd.at(i), v_upd.getElements().at(i)));
				values.insert(values.begin(), std::pair<Var,SR>(u.at(i), v.getElements().at(i)));
			}
			delta_new = Polynomial<SR>::eval(delta,values);

			v = v + v_upd;
			v_upd = step(poly_vars, J_s, v, delta_new);
		}

		return v;
	}
};

#endif
