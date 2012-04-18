# Test inputs :)
var('x,y,z')
var('a,b,c,d,e,f,g,h,i')

f1 = a*x*y + b
f2 = c*y*z + d*y*x + e
f3 = g*x*h + i

# to compute the termination probability of the above system:
probabilistic_subs = dict( [(a,0.4),(b,0.6),(c,0.3),(d,0.4),(e,0.3),(g,0.3),(h,1),(i,0.7) ] )
#probabilistic_subs = dict( [(a,2/5),(b,3/5),(c,3/10),(d,2/5),(e,3/10),(g,3/10),(h,1),(i,7/10) ] )

F = vector(SR,[f1,f2,f3]).transpose()

F_c = F.subs(probabilistic_subs)

variables = [x,y,z]
F_diff = F_c - vector(SR,variables).transpose()


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

