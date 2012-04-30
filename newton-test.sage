# Test inputs :)
var('x,y,z')
var('a,b,c,d,e,f,g,h,i')

f1 = a*x*y + b
f2 = c*y*z + d*y*x + e
f3 = g*x*h + i

# to compute the termination probability of the above system:
probabilistic_subs = dict( [(a,0.4),(b,0.6),(c,0.3),(d,0.4),(e,0.3),(g,0.3),(h,1),(i,0.7) ] )
#probabilistic_subs = dict( [(a,2/5),(b,3/5),(c,3/10),(d,2/5),(e,3/10),(g,3/10),(h,1),(i,7/10) ] )

F = vector(SR,[f1,f2,f3]).column()

F_c = F.subs(probabilistic_subs)

variables = [x,y,z]
F_diff = F_c - vector(SR,variables).column()


# test functions

# test the fixpoint method with i steps
def test_newton_fixpoint(i) :
    return newton_fixpoint_solve(F_c, variables, i)

# test the numerical approach with i steps
def test_newton_numeric(i) :
    return newton_numerical(F_diff, variables, i)

# return the difference between both method after i steps
def newton_error_at(i) :
    fp = test_newton_fixpoint(i)
    num = test_newton_numeric(i)
    return num - fp

# wrapper to see the difference between both methods after each step
def test_newton(max_iter=20) :
    for i in range(1,max_iter+1):
        error = newton_error_at(i)
        print "Step ", i, ": error (min,max) = ", error.min(), error.max()


# generate a (sparse) system of n quadratic equations in n variables with real coefficients in (0,1)
# return the system (as a vector of symbolic expressions) as well as its variables
# eps is the "density", i.e. eps*n of the coefficients will be non-zero
import numpy.random
import itertools
import random

def gen_random_quadratic_eqns(n, eps) :
    poly_vars = list( var(join(['x%d' %i for i in range(n)])) )

    # n variables ==> n choose 2 monomials, select eps*n of those which will be non-zero
    monomials = [m for m in itertools.combinations(poly_vars,2)]

    k = int(eps*n)
    
    #TODO: make coefficients symbolic constants
    F = []

    #for the non-zero monomials choose random coefficients from (0,1) that add up to 1
    # this can be done by sampling from a (/the) k-dimensional (symmetric) Dirichlet-distribution
    
    for i in range(n) :
        non_zero = random.sample(monomials, k)
        coeff = numpy.random.dirichlet([1]*k)
        f = sum(map(lambda x : x[0]*x[1][0]*x[1][1], zip(coeff,non_zero) ) )
        F.append(f)

    return F,poly_vars







