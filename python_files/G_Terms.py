import numpy as np

def Build_G1(u, v, dx, dy, nu):
    """
    Calculates the G1 term for time marching NS Solver
    
    Params:
    -------
    u, v     2D array of float, u and v velocities
    dx, dy   float, step size in x and y
    nu       float, viscosity of the fluid
    
    Returns:
    -------
    Gu       2D array of float, G in x direction
    
    """
    
    Gu = np.zeros_like(u, dtype=float)  
        
    Gu[1:-1, 1:-1] = \
    (u[1:-1,2:] - 2 * u[1:-1, 1:-1] + u[1:-1,:-2]) / dx**2 \
  + (u[2:,1:-1] - 2 * u[1:-1, 1:-1] + u[:-2,1:-1]) / dy**2 \
  - 0.25 * ((u[1:-1,2:] + u[1:-1, 1:-1])**2 -
            (u[1:-1, 1:-1] + u[1:-1,:-2])**2) / dx \
  - 0.25 * ((u[2:,1:-1] + u[1:-1, 1:-1])*(v[2:-1,:-1]+v[2:-1,1:]) -
            (u[:-2,1:-1] + u[1:-1, 1:-1])*(v[1:-2,:-1]+v[1:-2,1:])) / dy
      
    return Gu
    
def Build_G2(u, v, dx, dy, nu):
    """
    Calculates the G2 term for time marching NS Solver
    
    Params:
    -------
    u, v     2D array of float, u and v velocities
    dx, dy   float, step size in x and y
    nu       float, viscosity of the fluid
    
    Returns:
    -------
    Gv       2D array of float, G in y direction
    
    """
    Gv = np.zeros_like(v, dtype=float)

    Gv[1:-1, 1:-1] = \
    (v[1:-1,2:] - 2 * v[1:-1, 1:-1] + v[1:-1,:-2]) / dx**2 \
  + (v[2:,1:-1] - 2 * v[1:-1, 1:-1] + v[:-2,1:-1]) / dy**2 \
  - 0.25 * ((u[1:,2:-1] + u[:-1,2:-1])*(v[1:-1,2:]+v[1:-1,1:-1]) -
            (u[1:,1:-2] + u[:-1,1:-2])*(v[1:-1, 1:-1]+v[1:-1,:-2])) / dx \
  - 0.25 * ((v[2:,1:-1] + v[1:-1, 1:-1])*(v[2:,1:-1] + v[1:-1, 1:-1]) -
            (v[:-2,1:-1]+ v[1:-1, 1:-1])*(v[:-2,1:-1]+v[1:-1, 1:-1])) / dy
    
    return Gv
