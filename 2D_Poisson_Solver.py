def 2DPoisson_Solver(x_s,x_e,y_s,y_e,N_x,N_y):
	"""Solves the 2D Poisson Equation using an implicit method 
	on a uniform staggered grid

	Params:
	-------
	x_s, x_e: x direction starting and ending points
	y_s, y_e: y direction starting and ending points
	N_x, N_y: Number of points in x and y direction

	Returns:
	-------
	x:        x solution to the 2D Poisson Equation
	"""

#Grid Building
x_s, y_s = 0, 0 #starting points of x and y
x_e, y_e = 1.0, 1.0 #ending points of x and y
L_x, L_y = x_e - x_s, y_e - y_s #length of domain
N_x, N_y = 5, 5 #number of points
N_GC = 1 #thickness of ghost cell boundary

N_fx = N_x + 1 #number of face cells in x
N_cx = N_x     #number of center cells in x
N_fy = N_y + 1 #number of face cells in y
N_cy = N_y     #number of center cells in y

xf = np.zeros((N_fx,1),dtype=float) #initialize x face points
dxc = np.zeros((N_cx,1),dtype=float) #initializing x center delta points
xc = np.zeros_like(dxc,dtype=float) #intializing x center points

yf = np.zeros((N_fy,1),dtype=float) #ditto
dyc = np.zeros((N_cy,1),dtype=float)
yc = np.zeros_like(dyc,dtype=float) 

#x grid
for i in range(0,N_fx):
    xf[i] = x_s + (i*(L_x/N_x))

for i in range(0,N_cx):
    dxc[i] = xf[i+1] - xf[i] #delta x for face centers
    xc[i] = xf[i] + 0.5*dxc[i]
    
dx1 = dxc.copy() #setting face and center dx values equal
dx = np.insert(dx1,0,dx1[0]) #or append

#y grid
for j in range(0,N_fy):
    yf[j] = y_s + (j*(L_y/N_y))
    
for j in range(0,N_cy):
    dyc[j] = yf[j+1] - yf[j] #delta y for centers
    yc[j] = yf[j] + 0.5*dyc[j]
    
dy1 = dyc.copy() #setting face and center dy values equal
dy = np.insert(dy1,0,dy1[0])

#Ghost cells
N_cxGC = N_cx + 2.0*N_GC
N_cyGC = N_cy + 2.0*N_GC

top = N_GC + 1
bottom = N_GC + N_cx
left = N_GC + 1
right = N_GC + N_cy

[X,Y] = np.meshgrid(xc,yc)

#X grid building values for A
a_xd = np.zeros((N_cy,N_cx),dtype=float)
b_xd = np.zeros_like(a_xd,dtype=float)
c_xd = np.zeros_like(a_xd,dtype=float)
for j in range(0,N_cy):
    for i in range(0,N_cx):
        step = 1.0/(dx[i])
        a_xd[j,i] = step/dx[i]
        b_xd[j,i] = -step * (2.0*dx[i]) / (dx[i+1]*dx[i])
        c_xd[j,i] = step/dx[i]

#Y grid ditto
a_yd = np.zeros_like(a_xd,dtype=float)
b_yd = np.zeros_like(a_xd,dtype=float)
c_yd = np.zeros_like(a_xd,dtype=float)
for j in range(0,N_cy):
    for i in range(0,N_cx):
        step = 1.0/(dy[i])
        a_yd[j,i] = step/dy[i]
        b_yd[j,i] = -step * (2.0*dy[i]) / (dy[i+1]*dy[i])
        c_yd[j,i] = step/dy[i+1]

