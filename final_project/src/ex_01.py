from diffu_nonlin_module import *

domain = UnitDomain([10,10,10])
alpha = Constant(1.0)
f = Constant(0.0)
I = Expression("cos(pi*x[0])")
p = 1
rho = 1.0
dt = 0.05
T = 25

solver(dt, T, domain, p, rho, alpha, f, I, True)
