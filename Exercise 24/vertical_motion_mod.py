"""
	## Variables ##

	dt	= timestep
	m 	= mass of body (0.43 kg)
	g 	= gravitional acceleration
	rho = fluid density (1000kg/m^3) (8.9e-4 Pas => 1.0e3 kg/m^3)
	rhob= body density
	V	= volume of body (V = (4/3) * PI * r^3)
	r	= radius (11cm)
	d 	= diameter of body perpendicular to flow of fluid (r*2)
	v	= velocity of fluid
	mu	= dynamic fluid viscosity (3.2039999999999997 kg/mh) (0.00089 Pas)
	CD	= drag coefficient (0.45 for sphere)
	A	= cross-selection produced by cut plane perpendicular to the motion (PI * (r*2))
	
	## Functions ##
	
	Fg 	= - m * g
	Fb	= rho * g * V
	
	Re = (rho * d * abs(v)) / mu
	
	# Drag based on Reynolds number
	# when (Re < 1) 
	Fd(s) = -3 * PI * d * rho * v
	# when (Re > 10**3) 
	Fd(q) = -(1/2) * CD * mu * A * abs(v) * v
	
	
	# Newtons law of motion
	ma = Fg + Fd(s) + Fb
	ma = (-mg) - (3 * PI * d * rho * v) + (rho * g * V)
	
	# Differensial equation stokes drag:
	v' 	= -a * v + b
	a 	= (3 * PI * d * mu) / (rhob * V)
	b 	= g * ((rho / rhob) - 1)
	
	# apply CN scheme
	v^(n+1) = (v^n - (1/2) * dt * a * v^n + dt * b) / (1 + (1/2) * dt * a) 
	
	# Differensial equation quadratic drag:
	v' 	= -a * abs(v) * v + b
	a 	= (1/2) * CD * ((rho * A) / (rhob * V))
	b	= g * ((rho / rhob) - 1) 
	
	# apply CN scheme
	v^(n+1) = (v^n + dt * b^(n+(1/2))) / (1 + dt * a^(n+(1/2)) * abs(v^n))
"""

import numpy as np
from math import pi, sqrt
from numpy import exp

def solver(dt, T, m, g, rho, mu, r, CD, drag = True):
	T = float(T)
	dt = float(dt)
	Nt = int(round(T/dt))
	T =	Nt*dt
	v = np.zeros(Nt+1)
	t = np.linspace(0, T, Nt+1)
	
	V = (4/3) * pi * ((r * 0.01)**3)
	rhob = m / V
	d = (r * 0.01) * 2
	A = pi * (r * 0.01)**2
	
	# Start in idle position
	v[0] = 1
	for n in range(0, Nt):
		Re = (rho * d * abs(v[n])) / mu
		
		def b():
			return g * ((rho / rhob) - 1)
		
		# Use Stokes drag model
		if Re < 1 or not drag:
			if(drag):
				def a():
					return (3.0 * pi * d * mu) / (rhob * V)
			else:
				def a():
					return (rhob * V)
				
			v[n+1] = ((v[n] - 0.5 * dt * a() * v[n] + dt * b()) / (1.0 + 0.5 * dt * a()))
			
		# Use Quadratic drag model
		else:
			def a():
				return (0.5) * CD * ((rhob * V) / (rho * A))
			
			v[n+1] = ((v[n] + dt * b()) / (1.0 + dt * a() * abs(v[n])))
			
	return v, t

def exact_solution_neg_drag(t, m, g, rho, mu, r):
	V = (4/3) * pi * ((r * 0.01)**3)
	rhob = m / V
	d = (r * 0.01) * 2
	A = pi * (r * 0.01)**2
	
	def a():
		return (rhob * V)
	def b():
		return g * ((rho / rhob) - 1)
	
	return (b() * exp(-a() * t) * (exp(a()*t) - 1.0)) / a()

def verify_terminal_velocity(m, g, rho, mu, r, CD):
	V = (4/3) * pi * ((r * 0.01)**3)
	rhob = m / V
	d = (r * 0.01) * 2
	A = pi * (r * 0.01)**2

	def b():
		return g * ((rho / rhob) - 1)
	def a():
		return (3.0 * pi * d * mu) / (rhob * V)
			
	vT = b() / a()
	v1 = -a() * vT + b() 
	
	def a():
		return (0.5) * CD * ((rho * A) / (rhob * V))
		
	vT = sqrt(b() / a())
	v2 = -a() * abs(vT) * vT + b()
	
	return v1, v2
	
	
"""
m = len(dt_values)
r = [log(E_values[i-1]/E_values[i])/
                log(dt_values[i-1]/dt_values[i])
                for i in range(1, m, 1)]
"""	
def verify_convergence_rate(r):
	tol = 0.1
	expected_rate = 2
	r_final = r[-1]
	diff = abs(expected_rate - r_final)
	if diff > tol:
		return False
	return True 
