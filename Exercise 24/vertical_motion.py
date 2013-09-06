import numpy as np
from math import pi
from numpy import exp

class Problem:
	def __init__(self, m = 0.43, g = 9.81, rho = 1000, mu = 0.00089, r = 11, CD = 0.45):
		self.m, self.g, self.rho, self.mu, self.r, self.CD = m, g, rho, mu, r, CD

	def define_command_line_options(self, parser = None):
		if parser is None:
			import argparse
			parser = argparse.ArgumentParser()

		parser.add_argument("--m", "--body_mass", 
				type=float, default=0.43, help="mass of body in kg")
		parser.add_argument("--g", "--gravity", 
				type=float, default=9.81, help="gravity in m/s^2")
		parser.add_argument("--rho", "--water_density", 
				type=float, default=1000, help="density of water in kg/m^3")
		parser.add_argument("--mu", "--dynamic_viscosity", 
				type=float, default=0.00089, help="dynamic viscosity of fluid in Pa-s")
		parser.add_argument("--r", "--radius", 
				type=float, default=11, help="radius of ball in cm")
		parser.add_argument("--CD", "--drag_coeficient", 
				type=float, default=0.45, help="Drag coeficient")		
						
		return parser

	def init_from_command_line(self, args):
		self.m, self.g, self.rho, self.mu, self.r, self.CD = args.m, args.g, args.rho, args.mu, args.r, args.CD

	def exact_solution(self, t):
		# compute exact solution at time		
		return t

class Solver:
	def __init__(self, problem, dt=0.5):
		self.problem = problem
		self.dt = dt

	def define_command_line_options(self, parser):
		if parser is None:
			import argparse
			parser = argparse.ArgumentParser()

		parser.add_argument("--dt", "--time_step_value", 
				type=float, default=0.5, help="time step value")

		return parser

	def init_from_command_line(self, args):
		self.dt = args.dt
	
	def solve(self):
		from vertical_motion_mod import solver
		self.v, self.t = solver(self.dt, self.problem.m, self.problem.g, self.problem.rho, self.problem.mu, self.problem.r, self.problem.CD)

	def error(self):
		# return numerical error
		return None

class Visualizer:
	def __init__(self, problem, solver):
		self.problem, self.solver = problem, solver

	def plot(self, plt = None):
		if plt is None:
			import scitools.std as plt

		plt.plot(self.solver.t, self.solver.v, 'b--o')
		"""
		t_e = np.linspace(0, 10, 1001)
		from vertical_motion_mod import exact_solution
		u_e = exact_solution(t_e, self.problem.m, self.problem.g, self.problem.rho, self.problem.mu, self.problem.r)
		plt.plot(t_e, u_e, 'b-')
		"""
		
		plt.hold("on")
		plt.xlabel("t")
		plt.ylabel("v")
		plt.savefig("%g.png" % (self.solver.dt))

		return plt

def main():
	problem = Problem()
	solver = Solver(problem)
	viz = Visualizer(problem, solver)

	parser = problem.define_command_line_options()
	parser = solver.define_command_line_options(parser)
	args = parser.parse_args()
	problem.init_from_command_line(args)
	solver.init_from_command_line(args)
	
	solver.solve()
	import matplotlib.pyplot as plt
	plt = viz.plot(plt = plt)
	E = solver.error()
	if E is not None:
		print "Error: %.4E" % E
	plt.show()	

main()
