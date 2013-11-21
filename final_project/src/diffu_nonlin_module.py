from dolfin import *
import time

class Domain:
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
	

def solver(dt, T, domain = Domain(), degree = 1, rho = 1.0, alpha = None, f = None, I = None, show = True):
	N = int(T/dt)
	
	if domain.dim() == 1:
		mesh = UnitIntervalMesh(domain.x())
		if f == None:
			f = lambda x,t: 1.0
		if I == None:
			I = lambda x: 1.0

		
	if domain.dim() == 2:
		mesh = UnitSquareMesh(domain.x(),domain.y())
		if f == None:
			f = lambda x,y,t: 1.0
		if I == None:
			I = lambda x,y: 1.0
			
	if domain.dim() == 3:
		mesh = UnitCubeMesh(domain.x(),domain.y(),domain.z())
		if f == None:
			f = lambda x,y,z,t: 1.0
		if I == None:
			I = lambda x,y,z: 1.0
	
	if alpha == None:
		alpha = lambda u: 1.0
	
	V = FunctionSpace(mesh, 'Lagrange', degree)
	
	#TODO compute real u0
	u0 = Expression('1')
	u0.t = 0
	ma
	u1 = interpolate(u0, V)
	
	u = TrialFunction(V)
	v = TestFunction(V)
	
	#TODO create f from python func
	f = Constant(1.0)
	
	#TODO insert a and L (variational form) from notes
	a = u*v*dx + dt*inner(nabla_grad(u), nabla_grad(v))*dx
	L = (u1 + dt*f)*v*dx + u0*v*ds

	A = assemble(a)   # assemble only once, before the time stepping
	
	u = Function(V)   # the unknown at a new time level
	t = dt

	while t <= T:
		b = assemble(L)
		u0.t = t
		solve(A, u.vector(), b)

		t += dt
		u1.assign(u)
		if show:
			plot(u)
		time.sleep(dt) # 3 seconds
				
	return u

if __name__ == "__main__":
	u = solver(0.1, 50);
