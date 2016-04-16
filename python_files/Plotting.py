def Pressure_Plotting(p, xp, yp, Nx, nt):
    p_min = np.min(p)
    p_max = np.max(p)

    %matplotlib inline
    size=8
    plt.figure(figsize=(size,size))
    plt.xlim(0,2.0*math.pi)
    plt.ylim(0,2.0*math.pi)
    plt.xlabel('x')
    plt.ylabel('y')
    contf = plt.contourf(xp, yp, p, levels = np.linspace(p_min\
        ,p_max,40))
    cbar = plt.colorbar(contf)
    cbar.set_label('Pressure', fontsize=16)
    cbar.set_ticks(np.linspace(p_min,p_max,7))
    loc = './Plots/Taylor_Green/Pressure_Plot_%1.0f_Sparse.png' % Nx
    #loc = './Plots/Taylor_Green/Pressure_Plot_%1.0f.png' % Nx
    #loc = './Plots/Taylor_Green/Pressure_Plot_%1.0f_%1.0f.png' % (Nx, nt)
    plt.savefig(loc);
    
def u_Plotting(u, xu, yu, Nx, nt):
    u_min = np.min(u)
    u_max = np.max(u)

    %matplotlib inline
    size=8
    plt.figure(figsize=(size,size))
    plt.xlim(0,2.0*math.pi)
    plt.ylim(0,2.0*math.pi)
    plt.xlabel('x')
    plt.ylabel('y')
    contf = plt.contourf(xu, yu, u, levels = np.linspace(u_min\
        ,u_max,40))
    cbar = plt.colorbar(contf)
    cbar.set_label('u-Velocity', fontsize=16)
    cbar.set_ticks(np.linspace(u_min,u_max,7))
    loc = './Plots/Taylor_Green/u_Plot_%1.0f_Sparse.png' % Nx
    #loc = './Plots/Taylor_Green/u_Plot_%1.0f.png' % Nx
    #loc = './Plots/Taylor_Green/u_Plot_%1.0f_%1.0f.png' % (Nx, nt)
    plt.savefig(loc);
    
def v_Plotting(v, xv, yv, Nx, nt):
    v_min = np.min(v)
    v_max = np.max(v)

    %matplotlib inline
    size=8
    plt.figure(figsize=(size,size))
    plt.xlim(0,2.0*math.pi)
    plt.ylim(0,2.0*math.pi)
    plt.xlabel('x')
    plt.ylabel('y')
    contf = plt.contourf(xv, yv, v, levels = np.linspace(v_min\
        ,v_max,40))
    cbar = plt.colorbar(contf)
    cbar.set_label('v-Velocity', fontsize=16)
    cbar.set_ticks(np.linspace(v_min,v_max,7))
    loc = './Plots/Taylor_Green/v_Plot_%1.0f_Sparse.png' % Nx
    #loc = './Plots/Taylor_Green/v_Plot_%1.0f.png' % Nx
    #loc = './Plots/Taylor_Green/v_Plot_%1.0f_%1.0f.png' % (Nx, nt)
    plt.savefig(loc);
    
def Vorticity_Plotting(u, v, vorticity, xv, yv, Nx, nt):
    vort_min = np.min(vorticity)
    vort_max = np.max(vorticity)

    %matplotlib inline
    size=8
    plt.figure(figsize=(size,size))
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.xlim(0.0,2.0*math.pi)
    plt.ylim(0.0,2.0*math.pi)
    plt.title('Taylor-Green Vortices \n Vorticity and Velocity Vectors', fontsize=20)
    plt.quiver(xv, yv[1:] ,u[:,1:], v[1:,:], units='width', zorder=2)
    contf = plt.contourf(xv, yv[1:], vorticity, levels = np.linspace(vort_min\
        ,vort_max,40))
    cbar = plt.colorbar(contf)
    cbar.set_label('Vorticity', fontsize=16)
    cbar.set_ticks(np.linspace(vort_min,vort_max,7))
    loc = './Plots/Taylor_Green/Vorticity_Plot_%1.0f_Sparse.png' % Nx
    #loc = './Plots/Taylor_Green/Vorticity_Plot_%1.0f.png' % Nx
    #loc = './Plots/Taylor_Green/Vorticity_Plot_%1.0f_%1.0f.png' % (Nx, nt)
    plt.savefig(loc);