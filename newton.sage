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
# TODO: perhaps add an evaluation function for "s" (so that we have s(0) = 1 etc.)?
from sage.symbolic.function_factory import function
#s = function('s', nargs=1)
s = 1/(1-x)

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
def nhessian(self, variables=None):
    from sage.matrix.constructor import matrix
    if variables is None:
        variables = self.arguments()
    return matrix([[g.derivative(x) for x in variables] for g in self.gradient(variables)])


# compute symbolic delta with a parameter vector u
def compute_symbolic_delta(u,F,var) :
    d = []
    for i in F: # iterate over equations
        H = nhessian(i[0], var)
        d.append((1/2*u.transpose()*H*u)[0,0]) #TODO: make this nice :)
    return vector(SR,d).column()

# symbolic delta computation also for non-quadratic polynomials F
#def compute_symbolic_delta_general (u,F,var) :
    


# given a vector of polynomials F in variables poly_vars, its Jacobian,
# a starting value v, and the "update" delta, compute the update
# v_update = J^*|v * delta
# such that v_new = v + v_update

def newton_step(F, poly_vars, J_s, v, delta) :
    assert(len(poly_vars) == len(v.list()))

    sub_dict = dict(zip(poly_vars,v.list()))
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
    
    delta = compute_symbolic_delta(vector(u).column(),F,poly_vars)
    
    # v^0 = F(0)
#    v = F.subs( dict( (v,0) for v in poly_vars )) 
    v0 = matrix(SR,F.nrows(),1)
    sub = dict(zip(poly_vars,v0.list()))
    v_upd1 = J_s.subs(sub)* F.subs( dict( (v,0) for v in poly_vars ))
    v1 = v0 + v_upd1
    
    v_upd = v_upd1
    v=v1

    # TODO: iteration..
    for i in range(2,max_iter+1) :
        delta_new = delta.subs( dict( zip(u,v_upd.list()) ) )
        v_upd = newton_step(F,poly_vars,J_s,v,delta_new)
        v = v + v_upd

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
        x = linalg.solve(A,b) #TODO: andere Funktion... Gleichungssystem in Software lösen?
        v = v - x

    return v









