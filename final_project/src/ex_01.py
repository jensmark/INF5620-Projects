from diffu_nonlin_module import *

domain = UnitDomain([10])
alpha = lambda u: 1.0
f = Constant(0.0)
I = Expression("cos(pi*x[0])")
p = 1
rho = 1.0
dt = 0.0005
T = 25

solver(dt, T, domain, p, rho, alpha, f, I, 0.0, True)
