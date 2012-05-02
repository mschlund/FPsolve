from sage.symbolic.ring import SR

# n = newton_fixpoint_solve(F, [x,y,z], 10)
# n.subs( a=0.5, b=0.3 )
# TODO (plan, long run):
# substitute either functions (p_(i,j)(n) = transition probabilities) or constants for some parameters
# and substitute the function f(x) = 1/(1-x) for s() to evaluate :) --> we get a rational expression in the free parameters
# Ultimately: Use it to approximate e.g. the expected value of the "runtime" (or some other "cost") as well as its variance (in terms of the input size n)
# Also:
# - Derive the equations from Flowgraphs, Remopla/PDS (which are derived from e.g. Java Bytecode)
# - More semirings (other analyses), integrate semirings, polynomials over them, and Newton's method nicely into Sage


# symbolic functions for the Kleene star of elements
from sage.symbolic.function_factory import function
s = function('s', nargs=1, eval_func=lambda self, x: 1 if x == 0 else s(x,hold=True))
#s = 1/(1-x)

from sage.matrix.constructor import block_matrix
def compute_mat_star(M) :
    if M.nrows() == 1 and M.ncols() == 1: # just a scalar in a matrix
        return s(M[0,0]) # apply the Kleene star to the only element and return
    else: # there is a matrix to deconstruct
        M = copy(M)
        if M.nrows() % 2 == 0 and M.ncols() % 2: # even number of rows/columns, split in half
            M.subdivide(M.nrows()/2,M.ncols()/2)
        else:   # odd number of rows/columns, peeling mode
            M.subdivide(M.nrows()-1,M.ncols()-1)
        a_11 = M.subdivision(0,0)
        a_12 = M.subdivision(0,1)
        a_21 = M.subdivision(1,0)
        a_22 = M.subdivision(1,1)
        as_11 = compute_mat_star(a_11)
        as_22 = compute_mat_star(a_22)
        # matrix divided and precalculated, now do the rest
        A_11 = compute_mat_star(a_11 + a_12 * as_22 * a_21)
        A_22 = compute_mat_star(a_22 + a_21 * as_11 * a_12)
        A_12 = as_11 * a_12 * A_22
        A_21 = as_22 * a_21 * A_11
        return block_matrix(SR,[[A_11,A_12],[A_21,A_22]],  subdivide=False)

# taken from sage/symbolic/expression.pyx
# modified it to be able to specify the variables
# should be fixed in lib
def nhessian(self, poly_vars=None):
    from sage.matrix.constructor import matrix
    if poly_vars is None:
        poly_vars = self.arguments()
    return matrix([[g.derivative(x) for x in poly_vars] for g in self.gradient(poly_vars)])


# compute symbolic delta with a parameter vector u
# NOTE: this is for quadratic polynomials only!!
def compute_symbolic_delta(u,F,v) :
    d = []
    for i in F: # iterate over equations
        H = nhessian(i[0], v)
        d.append((1/2*u.transpose()*H*u)[0,0]) #TODO: make this nice :)
    return vector(SR,d).column()

# compute the degree of f when viewed as a polynomial in the variables v
# (necessary since some "symbolic variables" may actually be symbolic constants for us
def compute_degree(f,poly_vars) :
    variables = poly_vars[:]
    maxdeg = sum([f.degree(t) for t in variables])
    deg = 0
    while deg < maxdeg and len(variables) > 0:
        x = variables.pop()
        degx = f.degree(x)
        deg += degx
        f = f.diff([x]*degx)
        if (f==0) :
            break
    return deg

import itertools as it
# symbolic delta computation also for non-quadratic polynomials F
# parameter vectors are v and v_upd (which may also be concrete values if wanted :))
# v should be the (d-1)-st newton iterand and v_upd should be the (d)-th newton-update
def compute_symbolic_delta_general(v, v_upd, F, poly_vars) :
    n = len(v)
    assert(len(poly_vars) == n)
    delta = vector(SR,n)
    for i in range(n) :
        f = F[i][0]
        deg = compute_degree(f, poly_vars)

        for idx in it.product(range(0,deg+1), repeat=n) :
            if sum(idx) <=deg and sum(idx) >= 2 :
#                print str(idx)
                dx = reduce(lambda x,y : x + ([y[0]]*y[1]), zip(poly_vars,idx), [])
                prod = reduce(lambda p,x : p*(x[0]**x[1]), zip(v_upd,idx), 1)

                sub_dict = dict(zip(poly_vars,v))
                delta[i] = delta[i] + diff(f,dx).subs(sub_dict) * prod
    return delta.column()


# given a vector of polynomials F in variables poly_vars, its Jacobian,
# a starting value v, and the "update" delta, compute the update
# v_update = J^*|v * delta
# such that v_new = v + v_update

def newton_step(F, poly_vars, J_s, v, delta) :
    assert(len(poly_vars) == v.nrows())

    sub_dict = dict(zip(poly_vars,v.list()) )
    J_s = J_s.subs(sub_dict)

#    v_new = v + J_s*delta
    v_upd = J_s*delta
    return v_upd


# TODO: iterate until convergence, but for at most max_iter iterations
# 1) compute the concrete delta (from the precomputed symbolic expression)
# 2) execute a newton_step

from sage.calculus.functions import jacobian 
def newton_fixpoint_solve(F, poly_vars, max_iter=10) :
    J = jacobian(F, poly_vars)
    J_s = compute_mat_star(J) #only compute matrix star once

    u = var(join(['u%d' %i for i in range(J.ncols())]))
    u_upd = var(join(['u_upd%d' %i for i in range(J.ncols())]))

#    delta = compute_symbolic_delta(vector(u).column(),F,poly_vars)
    delta = compute_symbolic_delta_general(u, u_upd, F, poly_vars)

    v = matrix(SR,F.nrows(),1) # v^0 = 0
    delta_new = F.subs( dict( (v,0) for v in poly_vars ))

    v_upd = newton_step(F,poly_vars,J_s,v,delta_new)
#    v = v + v_upd
    
    # define symbolic variables for v^[i] and v^(i)
    v_s = matrix(var(join(['vs%d_%d' %(1,j) for j in range(v.nrows())]))).transpose()
    vu_s = matrix(var(join(['vus%d_%d' %(1,j) for j in range(v_upd.nrows())]))).transpose()
    # save the current values
    v_list = zip(v_s.list(), v.list())
    vu_list = zip(vu_s.list(), v_upd.list())

    # newton-iteration..
    for i in range(2,max_iter+1) :
 #       delta_new = delta.subs( dict( zip(u,v_upd.list()) ) )
        delta_new = delta.subs(dict( zip(u_upd,vu_s.list()) + zip(u, v_s.list())) )

        v = v_s + vu_s
        v_upd = newton_step(F,poly_vars,J_s,v_s,delta_new)

        v_s = matrix(var(join(['vs%d_%d' %(i,j) for j in range(v.nrows())]))).transpose()
        vu_s = matrix(var(join(['vus%d_%d' %(i,j) for j in range(v_upd.nrows())]))).transpose()
        v_list += zip(v_s.list(), v.list())
        vu_list += zip(vu_s.list(), v_upd.list())

    v_dict = dict(v_list)
    vu_dict = dict(vu_list)
    # try intelligent backsubstitution
    # but this is an ugly implementation
    # this evaluates the symbolic variables bottom up
    for i in range(1,max_iter):
        v_var = var(join(['vs%d_%d' %(i,j) for j in range(v.nrows())]))
        vu_var = var(join(['vus%d_%d' %(i,j) for j in range(v.nrows())]))
        for j in range(F.nrows()):
            v_new_var = var('vs%d_%d' %(i+1,j))
            vu_new_var = var('vus%d_%d' %(i+1,j))
            # substitute all occurrences 
            for variable in v_var:
                # replace each dictionary key with a new one where the symbolic variables are
                # substituted with the previous iteration
                v_dict[v_new_var] = v_dict[v_new_var].subs(dict([(variable,v_dict[variable])]))
                vu_dict[vu_new_var] = vu_dict[vu_new_var].subs(dict([(variable,v_dict[variable])]))
            for variable in vu_var:
                v_dict[v_new_var] = v_dict[v_new_var].subs(dict([(variable,vu_dict[variable])]))
                vu_dict[vu_new_var] = vu_dict[vu_new_var].subs(dict([(variable,vu_dict[variable])]))

    # last step
    v = v + v_upd
    v = v.subs(v_dict)
    v = v.subs(vu_dict)

    return v

# Newton's method as a textbook-implementation
# find a zero of the multivariate polynomial (vector) F starting with v_0
# TODO: nicefy... or rather throw away and write new :)
from numpy import linalg
import numpy
from sage.symbolic.ring import NumpyToSRMorphism

# x = newton_numerical(F_diff, [x,y,z])
# vgl. mit x_sym = newton_fixpoint_solve(F_c, [x,y,z])

def newton_numerical(F, poly_vars, max_iter=10) :
    np_to_SR = NumpyToSRMorphism(numpy.float64, SR)
    J = jacobian(F,poly_vars)
    
    v_0 = numpy.array([0,0,0])
    v = vector(map(np_to_SR,v_0)).column()
    for i in range(1,max_iter+1) :
        v = vector(map(np_to_SR,v)).column()
        sub = dict(zip(poly_vars,v.list()))
        A = J.subs(sub)
        b = F.subs(sub)
        x = linalg.solve(A,b) #TODO: solve system in software rather than with numpy (finite precision!)??
        v = v - x

    return v









