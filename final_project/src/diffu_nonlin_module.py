from dolfin import *
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
	

def solver(dt, T, domain = UnitDomain(), degree = 1, rho = 1.0, alpha = None, f = None, I = None, show = True):
	N = int(T/dt)
	
	if f == None:
		f = Constant(1.0)
	if I == None:
		I = Constant(1.0)
	if alpha == None:
		alpha = Constant(1.0)
	
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
	
	#TODO compute real u0
	u1 = interpolate(I, V)
	
	u = TrialFunction(V)
	v = TestFunction(V)
	
	#TODO insert variational formulation from notes
	F = u*v*dx + dt*inner(nabla_grad(u), nabla_grad(v))*dx - (u1 + dt*f)*v*dx + I*v*ds
	
	u = Function(V)   # the unknown at a new time level
	t = dt

	while t < T + DOLFIN_EPS:
		#TODO: Solve non-linearity
		
		A = assemble(lhs(F)) #Do we need to re-assemble matrix?
		b = assemble(rhs(F))
		
		f.t = t
		solve(A, u.vector(), b)

		t += dt
		u1.assign(u)
		if show:
			plot(u)
		time.sleep(dt) # 3 seconds
				
	return u

if __name__ == "__main__":
	print 'TODO: Implement module execution'
