import numpy as np
# Note, this advection scheme is positive definite given stable courant conditions
#  However, it may be worth it to look into Adams Bashforth
#  also Marker Based: http://www.msi.umn.edu/~lilli/taras-pepi.pdf


def F(l,r,U):
    """Calculate the donor cell flux function
    l = left gridcell scalar 
    r = right gridcell scalar
    U = Courant number (u*dt/dx)
    
    If U is positive, return l*U if U is negative return r*U
    By using the mathematical form instead of the logical form, 
    we can run on the entire grid simultaneously. 
    """
    # Logical form:
    # if U>0:
    #     return U*l
    # else:
    #     return r*U
    # Mathematical form:
    #       U (or 0) *l  + U (or 0)*r
    # Uabs=np.abs(U)
    return (U+np.abs(U))/2 * l + (U-np.abs(U))/2 * r

def advect(q,u,v,dx,dt,courant=False):
    """advect a scalar q by wind field (u,v)
    q = input scalar field (e.g. cloud water [kg])
    u,v = horizontal and vertical wind speeds [m/s] on a staggered grid. 
    dx = grid size [m]
    dt = time step [s] 
            note if the Courant number exceeds 0.5 this will be reduced
            automagically and run twice, the output after the complete 
            time step will still be returned (i.e. dt/2 will be run twice)
    S = spatially varying source/sink term (same units as input scalar term /s)
    
    Algorithm from MPDATA: 
        Smolarkiewicz, PK and Margolin, LG (1998)
            MPDATA: A Finite-Differnce Solver for Geophysical Flows
            Journal of Computational Physics 140, p459-480. CP985901
    """
    
    U=u*dt/dx
    V=v*dt/dx
    # ideally the Courant condition should be maximized outside of the advection scheme, 
    # but this is here just in case
    if not courant:
        if np.array(U).size>1:
            maxcourant=np.max((np.abs(U)+np.abs(V)))
        else:
            maxcourant=np.abs(U[1:,:])+np.abs(V[:,1:])
        if maxcourant>0.495:
            # this is unstable, so we need to reduce the timestep 
            # (recursively so we will continue to reduce it if necessary)
            print(str(dt)+' --> '+str(dt/2))
            advect(q,u,v,dx,dt/2)
            advect(q,u,v,dx,dt/2)
    
    # First Pass: standard upwind scheme
    q1=q.copy()
    if np.array(U).size>1:
        Ux0=U[1:-1,:-1]
        Ux1=U[1:-1,1:]
        Vy0=V[:-1,1:-1]
        Vy1=V[1:,1:-1]
        q[1:-1,1:-1]=(q1[1:-1,1:-1]
            -(F(q1[1:-1,1:-1],q1[1:-1,2:],Ux1)
            -F(q1[1:-1,:-2],q1[1:-1,1:-1],Ux0))
            -(F(q1[1:-1,1:-1],q1[2:,1:-1],Vy1)
            -F(q1[:-2,1:-1],q1[1:-1,1:-1],Vy0)))
    else:
        q[1:-1,1:-1]=(q[1:-1,1:-1]
            -(F(q[1:-1,1:-1],q[1:-1,2:],U)
            -F(q[1:-1,:-2],q[1:-1,1:-1],U))
            -(F(q[1:-1,1:-1],q[2:,1:-1],V)
            -F(q[:-2,1:-1],q[1:-1,1:-1],V)))
     
    # previously, I was using q1 as an intermediate, this isn't necessary if it is just a single step, revert to q1 for below
    # in fortran q1 will be required? 
    #bad=np.where((q1<0)|(~np.isfinite(q1)))
    #if len(bad[0])>0:
    #    q1[bad]=0
    
    # The above calculation is just the standard upwind formulations
    # below, calculate the diffusivity correction term of MPDATA
    # f=np.double(0.5)
    # 
    # # define A,B [l,r,u,b] (see MPDATA review reference)
    # # l,r,u,b = left, right, upper, bottom edges of the grid cells
    # Al=(q1[1:-1,:-2]-q1[1:-1,1:-1])/(q1[1:-1,:-2]+q1[1:-1,1:-1])
    # Ar=(q1[1:-1,2:]-q1[1:-1,1:-1])/(q1[1:-1,2:]+q1[1:-1,1:-1])
    # Au=(q1[:-2,1:-1]-q1[1:-1,1:-1])/(q1[:-2,1:-1]+q1[1:-1,1:-1])
    # Ab=(q1[2:,1:-1]-q1[1:-1,1:-1])/(q1[2:,1:-1]+q1[1:-1,1:-1])
    # 
    # q11=q1[2:,2:]+q1[2:,1:-1]
    # Br=0.5*((q11-q1[:-2,2:]-q1[:-2,1:-1])
    #     /(q11+q1[:-2,2:]+q1[:-2,1:-1]))
    # q11=q1[2:,:-2]+q1[2:,1:-1]
    # Bl=0.5*((q11-q1[:-2,:-2]-q1[:-2,1:-1])
    #     /(q11+q1[:-2,:-2]+q1[:-2,1:-1]))
    # q11=q1[:-2,2:]+q1[1:-1,2:]
    # Bu=0.5*((q11-q1[:-2,:-2]-q1[1:-1,:-2])
    #     /(q11+q1[:-2,:-2]+q1[1:-1,:-2]))
    # q11=q1[2:,2:]+q1[1:-1,2:]
    # Bb=0.5*((q11-q1[2:,:-2]-q1[1:-1,:-2])
    #     /(q11+q1[2:,:-2]+q1[1:-1,:-2]))
    # 
    # # compute diffusion correction U/V terms (see MPDATA review reference)
    # Uabs=np.abs(U)
    # # first find U/V terms on the grid cell borders
    # curUabs=(Uabs[1:-1,:-2]+Uabs[1:-1,1:-1])/2
    # curU=Ux0
    # curV=(V[1:-1,:-2]+V[1:-1,1:-1])/2
    # # then compute Ul
    # Ul=curUabs*(1-curUabs)*Al - 2*f*curU*curV*Bl
    # # compute Ur using the same two steps
    # curUabs=(Uabs[1:-1,2:]+Uabs[1:-1,1:-1])/2
    # curU=Ux1
    # curV=(V[1:-1,2:]+V[1:-1,1:-1])/2
    # Ur=curUabs*(1-curUabs)*Ar - 2*f*curU*curV*Br
    # # compute Vu
    # Vabs=np.abs(V)
    # curVabs=(Vabs[:-2,1:-1]+Vabs[1:-1,1:-1])/2
    # curV=Vy0
    # curU=(U[:-2,1:-1]+U[1:-1,1:-1])/2
    # Vu=curVabs*(1-curVabs)*Bu - 2*f*curU*curV*Bu
    # # compute Vb
    # curVabs=(Vabs[2:,1:-1]+Vabs[1:-1,1:-1])/2
    # curV=Vy1
    # curU=(U[2:,1:-1]+U[1:-1,1:-1])/2
    # Vb=curVabs*(1-curVabs)*Bb - 2*f*curU*curV*Bb
    # 
    # q[1:-1,1:-1]=(q1[1:-1,1:-1]
    #     -(F(q1[1:-1,1:-1],q1[1:-1,2:],Ul)
    #     -F(q1[1:-1,:-2],q1[1:-1,1:-1],Ur))
    #     -(F(q1[1:-1,1:-1],q1[2:,1:-1],Vu)
    #     -F(q1[:-2,1:-1],q1[1:-1,1:-1],Vb)))
    
    # return q
    
def advectvertical(q,w,dx,dt,courant=False):
    W=w*dt/dx
    # if the courant stability condition is not met, decrease dt and try again
    if not courant:
        if np.max(np.abs(W))>0.95:
            q1=advectvertical(q,w,dx,dt/2.0,S)
            q2=advectvertical(q1,w,dx,dt/2.0,S)
            return q2

    q1=q.copy()
    q[:,1:-1,:]=(q1[:,1:-1,:]
            - (F(q1[:,1:-1,:],q1[:,2:,:],W[:,1:-1,:])
            - F(q1[:,:-2,:],q1[:,1:-1,:],W[:,:-2,:])))
    # assume a no flow bottom boundary condition
    q[:,0,:]=q1[:,0,:]- F(q1[:,0,:],q1[:,1,:],W[:,0,:])
    # and assume the top layer flows in and out of an identical layer above it. 
    q[:,-1,:]=q1[:,-1,:]- ((q1[:,-1,:]*W[:,-1,:])-F(q1[:,-2,:],q1[:,-1,:],W[:,-2,:]))
    
    # this is the MPDATA diffusion correction term for 1D flow
    # U2=np.abs(U)-U**2
    # Vl=(U2[:-2]+U2[1:-1])/2*(q1[:-2]-q1[1:-1])/(q1[:-2]+q1[1:-1])
    # Vr=(U2[2:]+U2[1:-1])/2*(q1[2:]-q1[1:-1])/(q1[2:]+q1[1:-1])
    # q[1:-1]=(q1[1:-1]+S[1:-1]*dt/2
    #         - (F(q1[1:-1],q1[2:],Vr) - F(q1[1:-1],q1[:-2],Vl)))
    # return q

def upwind1d(q,u,dx,dt,S):
    """docstring for upwind"""
    U=u*dt/dx
    # if the courant stability condition is not met, decrease dt and try again
    if np.max(np.abs(U))>0.95:
        q1=upwind1d(q,u,dx,dt/2.0,S)
        q2=upwind1d(q1,u,dx,dt/2.0,S)
        return q2
    
    q[1:-1]=(q[1:-1]+S[1:-1]*dt/2 
            - (F(q[1:-1],q[2:],(U[1:-1]+U[2:])/2.0) 
            - F(q[:-2],q[1:-1],(U[1:-1]+U[:-2])/2.0)))
    return q

# def init():
#   q=np.abs(5-np.arange(15,dtype='f'))
#   q2=np.abs(5-np.arange(15,dtype='f'))
#   q[0]=q[-2]
#   q[-1]=q[1]
#   q2[0]=q2[-2]
#   q2[-1]=q2[1]
#   return q,q2


def advect1d(q,u,dx,dt,S,verbose=False):
    U=u*dt/dx
    # if the courant stability condition is not met, decrease dt and try again
    if np.max(np.abs(U))>0.95:
        q1=advect1d(q,u,dx,dt/2.0,S)
        q2=advect1d(q1,u,dx,dt/2.0,S)
        return q2
    
    # format string for printing verbose data at runtime
    fmt='{0[0]:.2f} {0[1]:.2f} {0[2]:.2f} {0[3]:.2f} {0[4]:.2f}'
    q1=q.copy()
    q1[1:-1]=(q[1:-1] #+S[1:-1]*dt/2 
            - (F(q[1:-1],q[2:],(U[1:-1]+U[2:])/2.0) 
            - F(q[:-2],q[1:-1],(U[1:-1]+U[:-2])/2.0)))
    U2=np.abs(U)-U**2
    
    # useful for debugging
    if verbose:
        print("qin = "+fmt.format(q[5:10]))
        print("qup = "+fmt.format(q1[5:10]))
        print("upf = "+fmt.format(-(F(q[1:-1],q[2:],(U[1:-1]+U[2:])/2.0)
                            - F(q[:-2],q[1:-1],(U[1:-1]+U[:-2])/2.0))[4:9]))
    
    Vl=(U2[:-2]+U2[1:-1])/2 * (q1[:-2]- q1[1:-1]) / (q1[:-2]+ q1[1:-1])
    Vr=(U2[2:]+U2[1:-1])/2  * (q1[2:] - q1[1:-1]) / (q1[2:] + q1[1:-1])
    
    q[1:-1]=(q1[1:-1] #+S[1:-1]*dt/2
            - (F(q1[1:-1],q1[2:],Vr) - F(q1[:-2],q1[1:-1],Vl)))

    if verbose:
        print("mpf = "+fmt.format(-(F(q1[1:-1],q1[2:],Vr) - F(q1[:-2],q1[1:-1],Vl))[4:9]))
        print("out = "+fmt.format(q[5:10]))
    
    return q