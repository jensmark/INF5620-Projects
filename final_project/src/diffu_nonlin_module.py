from dolfin import *
from numpy import *
import time
from nose.tools import *

class UnitDomain:
	def __init__(self, d = [1]):
		self.d = d
	
	def dim(self):
		return len(self.d)
	def x(self):
		return self.d[0]
	def y(self):
		return self.d[1]
	def z(self):
		return self.d[2]
	

def solver(dt, T, domain = UnitDomain(), degree = 1, rho = 1.0, alpha = None, f = None, I = None, action = None):
	N = int(T/dt)
	
	if f == None:
		f = Constant(1.0)
	if I == None:
		I = Constant(1.0)
	if alpha == None:
		alpha = lambda u: 1.0
	
	if domain.dim() == 1:
		mesh = UnitIntervalMesh(domain.x())	
	elif domain.dim() == 2:
		mesh = UnitSquareMesh(domain.x(),domain.y())	
	elif domain.dim() == 3:
		mesh = UnitCubeMesh(domain.x(),domain.y(),domain.z())
	else:
		print 'Unsupported domain dimentions'
		return

	V = FunctionSpace(mesh, 'Lagrange', degree)
	f = interpolate(f, V)
	I = interpolate(I, V)
	u_ = interpolate(I, V)
	u1 = interpolate(I, V)
	
	u = TrialFunction(V)
	v = TestFunction(V)
	g = Constant(0.0)           # Neuman boundary value.
	
	f.t = 0
	F = -inner(rho*((u-u1)/dt),v)*dx - inner(alpha(u_)*nabla_grad(u), nabla_grad(v))*dx \
		+ inner(f,v)*dx + alpha(u_)*g*v*ds
	
	u = Function(V)   # the unknown at a new time level
	t = dt
	
	eps = 1.0
	tol = 1.0E-5        # tolerance
	c = 0           	# iteration counter
	maxiter = 25        # max no of iterations allowed	
	
	A = assemble(lhs(F)) #Do we need to re-assemble matrix?
	
	while t < T + DOLFIN_EPS:
		# Picard iterations
		while eps > tol and c < maxiter:
			b = assemble(rhs(F))
			c += 1
			solve(A, u.vector(), b)
			diff = u.vector().array() - u_.vector().array()
			eps = linalg.norm(diff, ord=Inf)
			u_.assign(u)   # update for next iteration
		
		b = assemble(rhs(F))
		f.t = t
		solve(A, u.vector(), b)

		t += dt
		u1.assign(u)
		if action is not None:
			action(t, u)
		time.sleep(dt)
				
	return u, V
	
def test_convergance_rate():
	domain = UnitDomain([10,10])
	alpha = lambda u: 1.0
	f = Constant(0.0)
	I = Expression("cos(pi*x[0])")
	p = 1
	rho = 1.0
	dt = 0.05
	T = 0.05
	
	u, V = solver(dt, T, domain, p, rho, alpha, f, I, None)
	
	u_e = interpolate(Expression("pow(E,pow(-pi,2.0)*0.05)*cos(pi*x[0])",E=2.71828), V)
	
	e = u_e.vector().array() - u.vector().array()
	E = sqrt(sum(e**2)/u.vector().array().size)
	assert_almost_equal(E, pow(dt,p))
	
def test_mms():
	domain = UnitDomain([10])
	alpha = lambda u: 1.0 + pow(u,2.0)
	#f = Constant(0.0)
	I = Expression("cos(pi*x[0])")
	p = 1
	rho = 1.0
	dt = 0.05
	T = 0.5
	f = Expression("-rho*pow(x[0],3)/3 + rho*pow(x[0],2)/2 + 8*pow(t,3)*pow(x[0],7)/9 - 28*pow(t,3)*pow(x[0],6)/9+7*pow(t,3)*pow(x[0],5)/2 - 5*pow(t,3)*pow(x[0],4)/4 + 2*t*x[0] - t",t=0.0,rho=rho)
	
	def u_e(x,t):
		return (x**2.0)*(0.5 - x/3.0)*t
	
	def action(t, u):
		for x in range(10):
			assert_almost_equal(u_e(x,t),u.vector().array()[x])
	
	solver(dt, T, domain, p, rho, alpha, f, I, action)	

if __name__ == "__main__":
    import sys
    from scitools.misc import function_UI
    cmd = function_UI([solver,
                       test_convergance_rate, 
                       test_mms], sys.argv)
    eval(cmd)
