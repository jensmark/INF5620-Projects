from dolfin import *
from numpy import *
import time

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
	

def solver(dt, T, domain = UnitDomain(), degree = 1, rho = 1.0, alpha = None, f = None, I = None, b = 0.0, show = True):
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
	
	u_ = interpolate(I, V)
	u1 = interpolate(I, V)
	
	u = TrialFunction(V)
	v = TestFunction(V)
	
	f.t = 0
	plot(u1)
	F = -rho*((u-u1)/dt)*v*dx - inner(alpha(u_)*nabla_grad(u), nabla_grad(v))*dx \
		+ inner(f,v)*dx + alpha(u_)*I*v*ds
	
	u = Function(V)   # the unknown at a new time level
	t = dt
	
	eps = 1.0
	tol = 1.0E-5        # tolerance
	c = 0           	# iteration counter
	maxiter = 25        # max no of iterations allowed	
	
	while t < T + DOLFIN_EPS:
		# Picard iterations
		while eps > tol and c < maxiter:
			A = assemble(lhs(F)) #Do we need to re-assemble matrix?
			b = assemble(rhs(F))
			c += 1
			solve(A, u.vector(), b)
			diff = u.vector().array() - u_.vector().array()
			eps = linalg.norm(diff, ord=Inf)
			u_.assign(u)   # update for next iteration
		
		#u_ has changed, do we need to re-assemble matrix?
		A = assemble(lhs(F)) 
		b = assemble(rhs(F))
		f.t = t
		solve(A, u.vector(), b)

		t += dt
		u1.assign(u)
		if show:
			plot(u)
		time.sleep(dt*10)
				
	return u

if __name__ == "__main__":
	print 'TODO: Implement module execution'
