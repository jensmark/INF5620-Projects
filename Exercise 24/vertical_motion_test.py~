from nose import *
import numpy as np
import vertical_motion_mod as mod

def test_solver():
	dt = 0.1; T = 100; m = 0.43; g = 9.81; 
	rho = 1000; mu = 0.00089; r = 11; CD = 0.45
	plot = False
	
	v, t = mod.solver(dt, T, m, g, rho, mu, r, CD, False)
	t_e = np.linspace(0, 100, 1001)
	v_e = mod.exact_solution_neg_drag(t_e, m, g, rho, mu, r)
	
	if (plot):
		import matplotlib.pyplot as plt

		plt.plot(t, v, 'b--o')
		plt.plot(t_e, v_e, 'b-')
	
		plt.hold("on")
		plt.xlabel("t")
		plt.ylabel("v")
	
		plt.show()
	
	difference = abs(v_e - v).max() 
	tol = 1.0
	success = difference <= tol
	assert success == True 

def test_verify_terminal_velocity():
	m = 0.43; g = 9.81; rho = 1000; mu = 0.00089; r = 11; CD = 0.45
	v1, v2 = mod.verify_terminal_velocity(m, g, rho, mu, r, CD)
	assert v1 == 0
	assert v2 == 0
