def RK3(u, v, Nx, Ny, dx, dy, nu, nt):
    for t in range(0,nt):
        un = np.empty_like(u)
        vn = np.empty_like(v)
        #copy info from last loop, needed for current loop into variable_n
        un = u.copy()
        vn = v.copy()
        #Pn = p.copy()
        T = 0 + t*dt

        ####STEP 1
        ### t -----> t + (1/3)*dt
        un, vn = BCs(un, vn)
        G1 = Build_G1(un, vn, dx, dy, nu)
        G2 = Build_G2(un, vn, dx, dy, nu)
        u1 = un + (1.0/3.0*dt*G1)
        v1 = vn + (1.0/3.0*dt*G2)
        u1, v1 = BCs(u1, v1)
        #Solve for Pressure Field
        b1 = np.zeros((Ny+2,Nx+2), dtype=float)
        b1[1:-1,1:-1] = ((u1[1:-1,2:-1] - u1[1:-1,1:-2])/dx +\
            (v1[2:-1,1:-1] - v1[1:-2,1:-1])/dy)*(dx**2)
        b1 = np.reshape(b1[1:-1,1:-1], Ny*Nx)
        #A1 = Poisson_Solver_Neumann(Nx, Ny)
        A1 = Build_Sparse_A_Neumann(Nx, Ny)
        b1 = b1[:]*(3.0/dt)
        #temp1 = np.linalg.solve(A1,b1)
        temp1 = scipy.sparse.linalg.spsolve(A1,b1)
        
        p_star = np.reshape(temp1, (Ny,Nx))
        p_star = p_star[:,:] - p_star[Ny/2,Nx/2]

        F1_Pstar, F2_Pstar = Pressure_Terms(p_star, dx, dy)

        #calc predictor velocities
        u_star = u1.copy()
        v_star = v1.copy()

        u_star[1:-1,2:-2] = u1[1:-1,2:-2] - (1.0/3.0*dt*F1_Pstar)
        v_star[2:-2,1:-1] = v1[2:-2,1:-1] - (1.0/3.0*dt*F2_Pstar)

        ####STEP 2
        ### t + (1/3)*dt -----> t + (3/4)*dt
        u_star, v_star = BCs(u_star, v_star)
        G1_star = ((-5.0/9.0)*G1) + Build_G1(u_star, v_star, dx, dy, nu)
        G2_star = ((-5.0/9.0)*G2) + Build_G2(u_star, v_star, dx, dy, nu)
        u2 = u_star + ((15.0/16.0)*dt*G1_star)
        v2 = v_star + ((15.0/16.0)*dt*G2_star)
        u2, v2 = BCs(u2, v2)
        #Solve for Pressure Field
        b2 = np.zeros((Ny+2,Nx+2), dtype=float)
        b2[1:-1,1:-1] = ((u2[1:-1,2:-1] - u2[1:-1,1:-2])/dx +\
            (v2[2:-1,1:-1] - v2[1:-2,1:-1])/dy)*(dx**2)
        b2 = np.reshape(b2[1:-1,1:-1], Ny*Nx)
        #A2 = Poisson_Solver_Neumann(Nx, Ny)
        A2 = Build_Sparse_A_Neumann(Nx, Ny)
        b2 = b2[:]*(12.0/(5.0*dt))
        #temp2 = np.linalg.solve(A2,b2)
        temp2 = scipy.sparse.linalg.spsolve(A2,b2)
        p_dstar = np.reshape(temp2, (Ny,Nx))
        p_dstar = p_dstar[:,:] - p_dstar[Ny/2,Nx/2]

        F1_Pdstar, F2_Pdstar = Pressure_Terms(p_dstar, dx, dy)

        #calc predictor velocities
        u_dstar = u2.copy()
        v_dstar = v2.copy()

        u_dstar[1:-1,2:-2] = u2[1:-1,2:-2] - ((5.0/12.0)*dt*F1_Pdstar)
        v_dstar[2:-2,1:-1] = v2[2:-2,1:-1] - ((5.0/12.0)*dt*F2_Pdstar)

        ####STEP 3
        ### t + (3/4)*dt ------> t + dt
        u_dstar, v_dstar = BCs(u_dstar, v_dstar)
        G1_dstar = ((-153.0/128.0)*G1_star) + Build_G1(u_dstar, v_dstar, dx, dy, nu)
        G2_dstar = ((-153.0/128.0)*G2_star) + Build_G2(u_dstar, v_dstar, dx, dy, nu)
        u3 = u_dstar + ((8.0/15.0)*dt*G1_dstar)
        v3 = v_dstar + ((8.0/15.0)*dt*G2_dstar)
        u3, v3 = BCs(u3, v3)
        #Solve for the Pressure Field
        b3 = np.zeros((Ny+2,Nx+2), dtype=float)
        b3[1:-1,1:-1] = ((u3[1:-1,2:-1] - u3[1:-1,1:-2])/dx +\
            (v3[2:-1,1:-1] - v3[1:-2,1:-1])/dy)*(dx**2)
        b3 = np.reshape(b3[1:-1,1:-1], Ny*Nx)
        #A3 = Poisson_Solver_Neumann(Nx, Ny)
        A3 = Build_Sparse_A_Neumann(Nx, Ny)
        b3 = b3[:]*(4.0/dt)
        #temp3 = np.linalg.solve(A3,b3)
        temp3 = scipy.sparse.linalg.spsolve(A3,b3)
        p = np.reshape(temp3, (Ny,Nx))
        p = p[:,:] - p[Ny/2,Nx/2]

        #Calculate the Pressure Terms
        F1_P, F2_P = Pressure_Terms(p, dx, dy)

        #calc predictor velocities
        u = u3.copy()
        v = v3.copy()

        u[1:-1,2:-2] = u3[1:-1,2:-2] - ((1.0/4.0)*dt*F1_P)
        v[2:-2,1:-1] = v3[2:-2,1:-1] - ((1.0/4.0)*dt*F2_P)
        #u, v = BCs(u, v)

        if t % 100 == 0:
            print 'time step is %1.0f' % t
            
    return u, v, p