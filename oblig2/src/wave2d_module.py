#!/usr/bin/env python

import time
from scitools.std import *
from scitools import easyviz
from scitools.all import *

def solver(I, V, f, b, q, Lx, Ly, Nx, Ny, dt, T,
           user_action=None, version='scalar'):
	if version == 'compiled':
		print 'Pre-compiled loops not yet implemented'
		sys.exit(1)
	elif version == 'vectorized':
		advance = advance_vectorized
	elif version == 'scalar':
		advance = advance_scalar

	x = linspace(0, Lx, Nx+1)  # mesh points in x dir
	y = linspace(0, Ly, Ny+1)  # mesh points in y dir
	dx = x[1] - x[0]
	dy = y[1] - y[0]
	dx2 = dx**2
	dy2 = dy**2

	xv = x[:,newaxis]          # for vectorized function evaluations
	yv = y[newaxis,:]

	
	# Allow f and V to be None or 0
	if f is None or f == 0:
		f = (lambda x, y, t: 0) if version == 'scalar' else \
		    lambda x, y, t: zeros((x.shape[0], y.shape[1]))
		# or simpler: x*y*0
	if V is None or V == 0:
		V = (lambda x, y: 0) if version == 'scalar' else \
		    lambda x, y: zeros((x.shape[0], y.shape[1]))
	if q is None or q == 0:
		q = (lambda x, y: 1.0) 

	u   = zeros((Nx+1,Ny+1), order='C')   # solution array
	u_1 = zeros((Nx+1,Ny+1), order='C')   # solution at t-dt
	u_2 = zeros((Nx+1,Ny+1), order='C')   # solution at t-2*dt
	V_a = zeros((Nx+1,Ny+1), order='C')	  # initial condition V in array form
	f_a = zeros((Nx+1,Ny+1), order='C')	  # source term in array form
	q_a = zeros((Nx+1,Ny+1), order='C')
	q_a[:,:] = q(xv,yv);
		
	stability_limit = (1/float(q_a.max()))*(1/sqrt(1/dx2 + 1/dy2))
	if dt <= 0:
		safety_factor = -dt
		dt = safety_factor*stability_limit
	elif dt < stability_limit:
		print "Error, invalid timestep for simulating"
	
	print "T", T
	print "dt", dt
	Nt = int(round(T/float(dt)))
	t = linspace(0, Nt*dt, Nt+1)    # mesh points in time

	Ix = range(0, u.shape[0])
	Iy = range(0, u.shape[1])
	It = range(0, t.shape[0])

	import time; t0 = time.clock()          # for measuring CPU time
	# Load initial condition into u_1
	if version == 'scalar':
		for i in Ix:
		    for j in Iy:
		        u_1[i,j] = I(x[i], y[j])
	else: # use vectorized version
		u_1[:,:] = I(xv, yv)

	if user_action is not None:
		user_action(u_1, x, xv, y, yv, t, 0)

	# Special formula for first time step
	n = 0
	# Can use advance function with adjusted parameters (note: u_2=0)
	if version == 'scalar':
		u = advance(u, u_1, u_2, f, x, y, t, n,
		            b, q, dt, V, step1=True)

	else:  # use vectorized version
		f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
		V_a = V(xv, yv)
		u = advance(u, u_1, u_2, f_a, b, q_a, dt, V_a, step1=True)

	if user_action is not None:
		user_action(u, x, xv, y, yv, t, 1)

	u_2[:,:] = u_1; u_1[:,:] = u

	for n in It[1:-1]:
		if version == 'scalar':
		    # use f(x,y,t) function
		    u = advance(u, u_1, u_2, f, x, y, t, n, b, q, dt)
		else:
		    f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
		    u = advance(u, u_1, u_2, f_a, b, q_a, dt)

		if user_action is not None:
		    if user_action(u, x, xv, y, yv, t, n+1):
		        break

		u_2[:,:], u_1[:,:] = u_1, u

	t1 = time.clock()
	print "t0", t0
	print "t1", t1
	# dt might be computed in this function so return the value
	return dt, t1 - t0

def advance_scalar(u, u_1, u_2, f, x, y, t, n, b, q, dt,
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    dt2 = dt**2 
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dx2 = dx**2
    dy2 = dy**2
    B = b*(dt/2.0)
    
    if step1:
    	for i in Ix[1:-1]:
    		for j in Iy[1:-1]:
    			u[i,j] = u_1[i,j]+(dt2/4.0*dx2)*\
    				((q(x[i+1], y[j])+q(x[i], y[j]))*(u_1[i+1,j]-u_1[i,j])-\
    					(q(x[i], y[j])+q(x[i-1], y[j]))*(u_1[i,j]-u_1[i-1,j]))+\
    				(dt2/4.0*dy2)*\
    				((q(x[i], y[j+1])+q(x[i], y[j]))*(u_1[i,j+1]-u_1[i,j])-\
    					(q(x[i], y[j])+q(x[i], y[j-1]))*(u_1[i,j]-u_1[i,j-1]))
    			u[i,j] += dt*V(x[i], y[j])
    else:
    	for i in Ix[1:-1]:
    		for j in Iy[1:-1]:
    			u[i,j] = (1.0+B)**(-1)*(2.0*u_1[i,j]+u_2[i,j]*(B-1.0)+\
    				(dt2/2.0*dx2)*\
    					((q(x[i+1], y[j])+q(x[i], y[j]))*(u_1[i+1,j]-u_1[i,j])-\
    						(q(x[i], y[j])+q(x[i-1], y[j]))*(u_1[i,j]-u_1[i-1,j]))+\
					(dt2/2.0*dy2)*\
						((q(x[i], y[j+1])+q(x[i], y[j]))*(u_1[i,j+1]-u_1[i,j])-\
    						(q(x[i], y[j])+q(x[i], y[j-1]))*(u_1[i,j]-u_1[i,j-1])))+\
					dt*f(x[i], y[j], t[n])
    						
           
    # Boundary condition du/dn=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0

    return u

def advance_vectorized(u, u_1, u_2, f_a, b, q_a, dt,
                       V=None, step1=False):
    dt2 = dt**2
    if step1:
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    u_xx = u_1[:-2,1:-1] - 2*u_1[1:-1,1:-1] + u_1[2:,1:-1]
    u_yy = u_1[1:-1,:-2] - 2*u_1[1:-1,1:-1] + u_1[1:-1,2:]
    u[1:-1,1:-1] = D1*u_1[1:-1,1:-1] - D2*u_2[1:-1,1:-1] + \
                   Cx2*u_xx + Cy2*u_yy + dt2*f_a[1:-1,1:-1]
    if step1:
        u[1:-1,1:-1] += dt*V[1:-1, 1:-1]
    # Boundary condition u=0
    j = 0
    u[:,j] = 0
    j = u.shape[1]-1
    u[:,j] = 0
    i = 0
    u[i,:] = 0
    i = u.shape[0]-1
    u[i,:] = 0
    return u

import nose.tools as nt

def test_manufactored_solution():
	print "Not implemented"

def test_standing_damped_waves():
	print "Not implemented"

def test_standing_undamped_waves():
	print "Not implemented"

def test_constant_solution(verbose=False,version='scalar'):
	Lx = 2
	Ly = 2
	Nx = 4
	Ny = 4
	T = 5
	b = 0.0
    
	def q(x, y):
		return 1.0
		
	# Constant initial condition gives constant solution
	def I(x, y):
		return 0.2
	
	def action_u(u, x, xv, y, yv, t, n):
		if verbose:
			print "u", u
			print "t", t
	
	solver(I, None, None, b, q, Lx, Ly, Nx, Ny, -1, T, 
		user_action=action_u, version=version)
	
def test_cubic_solution():
	print "Not implemented"

def test_1d_plug_solution(verbose=False,version='scalar'):
	Lx = 2
	Ly = 2
	Nx = 4
	Ny = 4
	T = 5
	b = 0.0
    
	def q(x, y):
		return 1.0
		
	def I(x, y):
		if abs(x-Lx/2.0) > 0.1:
			return 0
		else:
			return 1
	
	def action_u(u, x, xv, y, yv, t, n):
		if verbose:
			print "u", u
			print "t", t
	
	solver(I, None, None, b, q, Lx, Ly, Nx, Ny, -1, T, 
		user_action=action_u, version=version)

def run_efficiency_tests(nrefinements=4):
    def I(x, y):
        return sin(pi*x/Lx)*sin(pi*y/Ly)
        
    def q(x, y):
    	return 1.5
	
    Lx = 10;  Ly = 10
    b = 0.5
    T = 100
    ##versions = ['scalar', 'vectorized', 'compiled']
    versions = ['scalar']
    print ' '*15, ''.join(['%-13s' % v for v in versions])
    for Nx in 15, 30, 60, 120:
        cpu = {}
        for version in versions:
            dt, cpu_ = solver(I, None, None, b, q, Lx, Ly, Nx, Nx,
                              -1, T, user_action=None,
                              version=version)
            cpu[version] = cpu_
        cpu_min = min(list(cpu.values()))
        if cpu_min < 1E-6:
            print 'Ignored %dx%d grid (too small execution time)' \
                  % (Nx, Nx)
        else:
            cpu = {version: cpu[version]/cpu_min for version in cpu}
            print '%-15s' % '%dx%d' % (Nx, Nx),
            print ''.join(['%13.1f' % cpu[version] for version in versions])

def run_Gaussian(plot_method=2, version='scalar', save_plot=False):
    """
    Initial Gaussian bell in the middle of the domain.
    plot_method=1 applies mesh function, =2 means surf, =0 means no plot.
    """
    # Clean up plot files
    for name in glob('tmp_*.png'):
        os.remove(name)

    Lx = 20
    Ly = 20
    b = 0.2
    
    def q(x, y):
		return 1.0
	
    def I(x, y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)

    if plot_method == 3:
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt
        from matplotlib import cm
        plt.ion()
        fig = plt.figure()
        u_surf = None

    def plot_u(u, x, xv, y, yv, t, n):
        if t[n] == 0:
            time.sleep(2)
        if plot_method == 1:
            mesh(x, y, u, title='t=%g' % t[n], zlim=[-1,1],
                 caxis=[-1,1])
        elif plot_method == 2:
			easyviz.surfc(xv, yv, u, title='t=%g' % t[n], zlim=[-1, 1],
				colorbar=True, colormap=hot(), caxis=[-1,1],
				shading='flat')
        elif plot_method == 3:
            print 'Experimental 3D matplotlib...under development...'
            #plt.clf()
            ax = fig.add_subplot(111, projection='3d')
            u_surf = ax.plot_surface(xv, yv, u, alpha=0.3)
            #ax.contourf(xv, yv, u, zdir='z', offset=-100, cmap=cm.coolwarm)
            #ax.set_zlim(-1, 1)
            # Remove old surface before drawing
            if u_surf is not None:
                ax.collections.remove(u_surf)
            plt.draw()
            time.sleep(1)
        if plot_method > 0:
            time.sleep(0) # pause between frames
            if save_plot:
                filename = 'tmp_%04d.png' % n
                savefig(filename)  # time consuming!

    Nx = 40; Ny = 40; T = 20
    dt, cpu = solver(I, None, None, b, q, Lx, Ly, Nx, Ny, -1, T,
                     user_action=plot_u, version=version)

def run_physical_problem(plot_method=2, version='scalar', save_plot=False):
    # Clean up plot files
	for name in glob('tmp_*.png'):
		os.remove(name)
	
	Lx = 20
	Ly = 20
	b = 0.1
	g = 9.81
	
	def H(x,y):
		return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)
	
	def q(x, y):
		return g*H(x,y)
	
	def I(x, y):
		return exp(-0.5*x**2)
	
	if plot_method == 3:
		from mpl_toolkits.mplot3d import axes3d
		import matplotlib.pyplot as plt
		from matplotlib import cm
		plt.ion()
		fig = plt.figure()
		u_surf = None
	
	def plot_u(u, x, xv, y, yv, t, n):
		if t[n] == 0:
			time.sleep(2)
		if plot_method == 1:
			mesh(x, y, u, title='t=%g' % t[n], zlim=[-1,1],
				 caxis=[-1,1])
		elif plot_method == 2:
			easyviz.surfc(xv, yv, u, title='t=%g' % t[n], zlim=[-1, 1],
				colorbar=True, colormap=hot(), caxis=[-1,1],
				shading='flat')
		elif plot_method == 3:
			print 'Experimental 3D matplotlib...under development...'
			#plt.clf()
			ax = fig.add_subplot(111, projection='3d')
			u_surf = ax.plot_surface(xv, yv, u, alpha=0.3)
			#ax.contourf(xv, yv, u, zdir='z', offset=-100, cmap=cm.coolwarm)
			#ax.set_zlim(-1, 1)
			# Remove old surface before drawing
			if u_surf is not None:
				ax.collections.remove(u_surf)
			plt.draw()
			time.sleep(1)
		if plot_method > 0:
			time.sleep(0) # pause between frames
			if save_plot:
				filename = 'tmp_%04d.png' % n
				savefig(filename)  # time consuming!

	Nx = 40; Ny = 40; T = 20
	dt, cpu = solver(I, None, None, b, q, Lx, Ly, Nx, Ny, -1, T,
		             user_action=plot_u, version=version)


if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI([run_efficiency_tests,
                       run_Gaussian, 
                       run_physical_problem, ], sys.argv)
    eval(cmd)
