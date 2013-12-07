from diffu_nonlin_module import *

def action(t, u):
	plot(u)

domain = UnitDomain([10,10])
beta = 0.5
alpha = lambda u: 1.0+beta*pow(u,2.0)
f = Constant(0.0)
sigma = 1.0
I = Expression("exp(-(1.0/2*pow(sigma,2.0))*(pow(x[0],2.0)+pow(x[1],2.0)))", sigma=sigma)
p = 2
rho = 1.0
dt = 0.0005
T = 25

solver(dt, T, domain, p, rho, alpha, f, I, action)
