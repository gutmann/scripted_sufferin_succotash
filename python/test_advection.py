#!/usr/bin/env python
from bunch import Bunch
import numpy as np
import matplotlib.pyplot as plt
import advect

def init_domain(nx=20,cfl=0.2, init="sine"):
    """docstring for init_domain"""
    # step function
    if init=="step":
        q=np.ones(nx,dtype='d')
        q[nx/2:]=0
    
    # sine function
    if init=="sine":
        q=np.sin((np.arange(nx)-1)/float(nx-2) * 2* np.pi)+1

    # V function
    if init.lower()=="v":
        q=np.abs(np.arange(nx)-(nx/2))/(nx/2.0) *2
    
    u=np.zeros(nx-1,dtype='d')+cfl
    
    return (u,q,nx)

def flux(q0,q1,u):
    """docstring for upwind"""
    f=q0*u
    backwards=np.where(u<0)
    
    f[backwards]=q1[backwards]*u[backwards]
    
    return f

def ftcs_flux(q,u):
    return u[1:] * (q[:-2] - q[2:]) /2

def ftcs_flux_fourth_order(q,u):
    return u[2:-1] * (4.0/3.0 * (q[1:-3] - q[3:-1]) /2 - 
                    1.0/3.0 * (q[:-4]  - q[4:])/4 )

def mpdata(nx=20, cfl=0.2, nrotations=1, 
            init="sine", delay=False,limit=None):
    (u,q,nx)=init_domain(nx=nx+1,cfl=cfl,init=init)
    nx2=nx-1
    q_up=q.copy()
    q0=q.copy()
    
    n_time_steps=int((nx-1.6)/u[0]) * nrotations
    print("ntimesteps="+str(n_time_steps))
    print("nx="+str(nx))
    
    
    for i in range(n_time_steps):
        
        f=flux(q_up[:nx-1],q_up[1:],u)
        q_up[1:-1]+=(f[:nx-2]-f[1:])
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        
        if i==0:
            q=advect.advect1d(q[1:],u,1,1,0)
        else:
            q=advect.advect1d(q,u,1,1,0)
        
        # wrap around boundary conditions
        q[0]=q[nx-2]
        q[nx-1]=q[1]
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        
        if limit!=None:
            q[q<limit]=limit
        
        if (i%int(nrotations/u[0]))==0:
            plt.clf()
            plt.plot(q,label="MPDATA", color="green")
            plt.plot(q_up,label="upwind",color="red")
            plt.title(i)
            plt.legend()
            plt.ylim(-0.1,2.1)
            plt.draw()
            if delay:
                _=raw_input()
        
    plt.clf()
    plt.plot(q,label="MPDATA", color="green")
    plt.plot(q_up,label="upwind",color="red")
    plt.title(i)
    plt.legend()
    plt.ylim(-0.1,2.1)
    plt.plot(q0,color="black",linewidth=2)
    plt.draw()
    
    print("Initial mean = {:.3}".format(q0.mean()))
    print(" Final mean  = {:.3}".format(q.mean()))

    

def rungakutta(order=4, space_order=2, nx=20, cfl=0.2, nrotations=1, 
                    init="sine", delay=False,limit=None):
    """docstring for main"""
    (u,q,nx)=init_domain(nx=nx,cfl=cfl,init=init)
    q0=q.copy()
    
    n_time_steps=int((nx-1.6)/u[0]) * nrotations
    print("ntimesteps="+str(n_time_steps))
    print("nx="+str(nx))
    print("order="+str(order))
    
    q_up=q.copy()
    for i in range(n_time_steps):
        f=flux(q_up[:nx-1],q_up[1:],u)
        q_up[1:-1]+=(f[:nx-2]-f[1:])
        
        f0=ftcs_flux(q,u)
        if space_order==4:
           f0[1:-1]=ftcs_flux_fourth_order(q,u)
        # f0=flux(q[:nx-1],q[1:],u)
        # f0=f0[:nx-2]-f0[1:]

        # first order
        if order==1:
            q[1:nx-1]+=f0
        # second order
        elif order==2:
            q1=q.copy()
            q1[1:nx-1]+=f0
            f2=ftcs_flux(q1,u)
            if space_order==4:
               f2[1:-1]=ftcs_flux_fourth_order(q1,u)
            f2-=f0
            # f2=flux(q1[:nx-1],q1[1:],u)
            # f2=(f2[:nx-2]-f2[1:]) - f0
            
            q[1:nx-1]=q1[1:nx-1] + f2/2.0
            
        # third order
        elif order==3:
            q1=q.copy()
            q1[1:nx-1]+=f0/3.0
            
            f2=ftcs_flux(q1,u)
            if space_order==4:
               f2[1:-1]=ftcs_flux_fourth_order(q1,u)
            f2-= 5.0*f0/9.0
            # f2=flux(q1[:nx-1],q1[1:],u)
            # f2=(f2[:nx-2]-f2[1:]) - 5.0*f0/9.0
            
            q2=q1.copy()
            q2[1:nx-1]+=15.0/16.0 * f2
            
            f3=ftcs_flux(q2,u)
            if space_order==4:
               f3[1:-1]=ftcs_flux_fourth_order(q2,u)
            f3-= 153.0/128.0*f2
            # f3=flux(q2[:nx-1],q2[1:],u)
            # f3=(f3[:nx-2]-f3[1:]) - 153.0/128.0*f2
            
            q[1:nx-1]=q2[1:nx-1] + 8.0*f3/15.0
        
        # fourth order
        elif order==4:
            # f0 = k0
            q1=q.copy()
            q1[1:nx-1]+=f0/2.0
            q1[0]=q1[nx-2]
            q1[nx-1]=q1[1]
            
            # f2 = k1
            f2=ftcs_flux(q1,u)
            if space_order==4:
               f2[1:-1]=ftcs_flux_fourth_order(q1,u)
            
            q1=q.copy()
            q1[1:nx-1]+=f2/2.0
            q1[0]=q1[nx-2]
            q1[nx-1]=q1[1]
            f0+=2.0*f2
            
            # f2 = k2
            f2=ftcs_flux(q1,u)
            if space_order==4:
               f2[1:-1]=ftcs_flux_fourth_order(q1,u)
            
            q1=q.copy()
            q1[1:nx-1]+=f2
            q1[0]=q1[nx-2]
            q1[nx-1]=q1[1]
            f0+=2.0*f2

            # f2 = k3
            f2=ftcs_flux(q1,u)
            if space_order==4:
               f2[1:-1]=ftcs_flux_fourth_order(q1,u)
            
            q1=q.copy()
            q1[1:nx-1]+=f2
            f0+=f2
            
            q[1:nx-1]=q[1:nx-1] + f0/6.0
        
        
        # wrap around boundary conditions
        q[0]=q[nx-2]
        q[nx-1]=q[1]
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        
        if limit!=None:
            q[q<limit]=limit
        
        if (i%np.ceil(nrotations/u[0]))==0:
            plt.clf()
            # plt.plot(f,label="f")
            # plt.plot(f1,label="f1")
            # plt.plot(f0,label="f0")
            plt.plot(q,label="RK-"+str(order), color="green")
            plt.plot(q_up,label="upwind",color="red")
            plt.title(i)
            plt.legend()
            plt.ylim(-0.1,2.1)
            plt.draw()
            if delay:
                _=raw_input()
        
        
    plt.clf()
    plt.plot(q,label="RK-"+str(order), color="green")
    plt.plot(q_up,label="upwind",color="red")
    plt.title(i)
    plt.legend()
    plt.ylim(-0.1,2.1)
    plt.plot(q0,color="black",linewidth=2)
    plt.draw()
    
    print("Initial mean = {:.3}".format(q0.mean()))
    print(" Final mean  = {:.3}".format(q.mean()))

def leapfrog(order=1, space_order=4, nx=20, cfl=0.2, nrotations=1, 
                    init="sine", delay=False,limit=None):
    """docstring for main"""
    (u,q,nx)=init_domain(nx=nx,cfl=cfl,init=init)
    q0=q.copy()
    
    n_time_steps=int((nx-1.6)/u[0]) * nrotations
    print("ntimesteps="+str(n_time_steps))
    print("nx="+str(nx))
    print("order="+str(order))
    
    
    q_up=q.copy()
    q_last=q.copy()
    q_new=q.copy()
    for i in range(n_time_steps):
        
        f=flux(q_up[:nx-1],q_up[1:],u)
        q_up[1:-1]+=(f[:nx-2]-f[1:])
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        if i==0:
            q_last=q_up.copy()
        
        f=ftcs_flux(q,u)
        if space_order==4:
           f[1:-1]=ftcs_flux_fourth_order(q,u) 
        
        q_new[1:-1]=q_last[1:-1] + 2*f
        q_last=q
        q=q_new
        
        # wrap around boundary conditions
        q[0]=q[nx-2]
        q[nx-1]=q[1]
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        
        if limit!=None:
            q[q<limit]=limit
        
        if (i%int(nrotations/u[0]))==0:
            plt.clf()
            plt.plot(q,label="Leap-frog", color="green")
            plt.plot(q_up,label="upwind",color="red")
            plt.title(i)
            plt.legend()
            plt.ylim(-0.1,2.1)
            plt.draw()
            if delay:
                _=raw_input()
        
    plt.clf()
    plt.plot(q,label="Leap-frog", color="green")
    plt.plot(q_up,label="upwind",color="red")
    plt.title(i)
    plt.legend()
    plt.ylim(-0.1,2.1)
    plt.plot(q0,color="black",linewidth=2)
    plt.draw()
    
    print("Initial mean = {:.3}".format(q0.mean()))
    print(" Final mean  = {:.3}".format(q.mean()))

def leapfrog_asselin_filter(gamma=0.2, space_order=4, nx=20, cfl=0.2, nrotations=1, 
                    init="sine", delay=False,limit=None):
    """docstring for main"""
    (u,q,nx)=init_domain(nx=nx,cfl=cfl,init=init)
    q0=q.copy()
    
    n_time_steps=int((nx-1.6)/u[0]) * nrotations
    print("ntimesteps="+str(n_time_steps))
    print("nx="+str(nx))
    print("gamma="+str(gamma))
    
    
    q_up=q.copy()
    q_last=q.copy()
    q_new=q.copy()
    for i in range(n_time_steps-1):
        
        # the upwind scheme
        f=flux(q_up[:nx-1],q_up[1:],u)
        q_up[1:-1]+=(f[:nx-2]-f[1:])
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        if i==0:
            q[:]=q_up[:]
        
        f=ftcs_flux(q,u)
        if space_order==4:
           f[1:-1]=ftcs_flux_fourth_order(q,u) 
        
        q_new[1:-1]=q_last[1:-1] + 2*f
        
        # asselin filter
        q_last = q + gamma*(q_last - 2*q + q_new)
        q=q_new
        
        # wrap around boundary conditions
        q[0]=q[nx-2]
        q[nx-1]=q[1]
        
        if limit!=None:
            q[q<limit]=limit
        
        if (i%int(nrotations/u[0]))==0:
            plt.clf()
            plt.plot(q,label="Leap-frog", color="green")
            plt.plot(q_up,label="upwind",color="red")
            plt.title(i)
            plt.legend()
            plt.ylim(-0.1,2.1)
            plt.draw()
            if delay:
                _=raw_input()
        
    plt.clf()
    plt.plot(q,label="Leap-frog", color="green")
    plt.plot(q_up,label="upwind",color="red")
    plt.title(i)
    plt.legend()
    plt.ylim(-0.1,2.1)
    plt.plot(q0,color="black",linewidth=2)
    plt.draw()
    
    print("Initial mean = {:.3}".format(q0.mean()))
    print(" Final mean  = {:.3}".format(q.mean()))



def adamsbashforth(order=1, space_order=2, nx=20, cfl=0.2, nrotations=1,
                     init="sine", delay=False,limit=None):
    """docstring for main"""
    (u,q,nx)=init_domain(nx=nx,cfl=cfl,init=init)
    q0=q.copy()
    
    n_time_steps=int((nx-1.6)/u[0]) * nrotations
    print("ntimesteps="+str(n_time_steps))
    print("nx="+str(nx))
    print("order="+str(order))
    
    f0=ftcs_flux(q,u)
    if space_order==4:
       f0[1:-1]=ftcs_flux_fourth_order(q,u) 
    # f0=flux(q[:nx-1],q[1:],u)
    # f0=f0[:nx-2]-f0[1:]
    f1=f0.copy()
    
    q_up=q.copy()
    for i in range(n_time_steps):
        f=flux(q_up[:nx-1],q_up[1:],u)
        q_up[1:-1]+=(f[:nx-2]-f[1:])
        
        f=ftcs_flux(q,u)
        if space_order==4:
           f[1:-1]=ftcs_flux_fourth_order(q,u) 
        # f=flux(q[:nx-1],q[1:],u)
        # f=f[:nx-2]-f[1:]

        # first order
        if order==1:
            q[1:nx-1]+=f
        # second order
        elif order==2:
            q[1:nx-1]+=(3.0*f - f1)/2.0
        # third order
        elif order==3:
            q[1:nx-1]+=(23.0*f - 16.0*f1 + 5.0*f0)/12.0
        
        # wrap around boundary conditions
        q[0]=q[nx-2]
        q[nx-1]=q[1]
        q_up[0]=q_up[nx-2]
        q_up[nx-1]=q_up[1]
        
        if limit!=None:
            q[q<limit]=limit
        
        if (i%int(nrotations/u[0]))==0:
            plt.clf()
            plt.plot(q,label="Adams", color="green")
            plt.plot(q_up,label="upwind",color="red")
            plt.title(i)
            plt.legend()
            plt.ylim(-0.1,2.1)
            plt.draw()
            if delay:
                _=raw_input()
        
        if i<n_time_steps-1:
            f0[:]=f1[:]
            f1[:]=f[:]
        
    plt.clf()
    plt.plot(q,label="Adams", color="green")
    plt.plot(q_up,label="upwind",color="red")
    plt.title(i)
    plt.legend()
    plt.ylim(-0.1,2.1)
    plt.plot(q0,color="black",linewidth=2)
    plt.draw()
    
    print("Initial mean = {:.3}".format(q0.mean()))
    print(" Final mean  = {:.3}".format(q.mean()))
    print("Upwind mean  = {:.3}".format(q_up.mean()))


main=adamsbashforth

if __name__ == '__main__':
    main()