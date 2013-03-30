from __future__ import print_function
import numpy as np
# import matplotlib.pyplot as plt

# from nc import NC_writer
import units
from bunch import Bunch
import advect
# import mp #mp has to be passed in from the driver so that init only has to be called once

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

# note these are both simple solutions that mostly work but could be done better(?)
def get_newt(atm,dz):
    return atm.t-9.8*dz/1000.0
def get_newp(atm,dz):
    rho = atm.p/(atm.ta*287.058)
    return atm.p-rho*g*dz

def get_real_t(atm,potentialT):
    return potentialT/((100000.0/atm.p)**(R/cp))

def get_potential_t(atm):
    return atm.ta*(100000.0/atm.p)**(R/cp)


dx=2000.0

# def test(realP=None):
#     import mygis
#     topo=mygis.read_nc('/Volumes/G-SAFE/headwaters/wrf_output/terrain/2km_wrf_input_d01',var='HGT')
#     z=topo.data[0,:,:]
#     atm=Bunch(sh=0.004,u=15.0,v=-2.0,rh=0.95,dp=269.0,ta=270.0,p=850.0)
#     return swim2d(z,atm,realP=realP)
    
# this is just a dummy wrapper around advect to allow parallelization
# pool.map only allows a single argument to the procedure it calls.  
def padvect(args):
    # inputq=args[0][100,100]
    # if inputq>200:
    #     print("advecting? ",args[0][100,100])
    if len(args)==5:
        advect.advect(args[0],args[1],args[2],args[3],args[4],courant=True)
    elif len(args)==4:
        advect.advectvertical(args[0],args[1],args[2],args[3],courant=True)
    return args[0]
    # if inputq>200:
    #     print("advected:  ",args[0][100,100])

def swim2d(z,atm,oldatm,newp,newu,newv,neww,mp,dTdt,processPool, timestep=None):# ,fname=None):
    # dTdt=np.array(dTdt.transpose([1,0,2]),order="F")
    N=z.shape
    pos=0
    # if fname!=None:
    #     ncout=NC_writer(fname,N[1],N[0],atm.qv.shape[0])
    #     ncout.addVar('qv')
    #     # ncout.addVar('qc')
    #     # ncout.addVar('qr')
    #     # ncout.addVar('nr')
    #     # ncout.addVar('qi')
    #     # ncout.addVar('ni')
    #     # ncout.addVar('qg')
    #     # ncout.addVar('qs')
    #     ncout.addVar('th')
    if timestep==None:
        timestep=3*60.0*60.0 # seconds
    th=np.array(oldatm.th.transpose([1,0,2]), order='F') # theta = potential temperature
    # all of these should be mixing ratios kg/kg-air (not kg/kg-dry-air)
    qv=np.array(oldatm.qv.transpose([1,0,2]),order='F') # vapor pressure
    qc=np.array(oldatm.qc.transpose([1,0,2]),order='F') # cloud water
    qr=np.array(oldatm.qr.transpose([1,0,2]),order='F') # rain in air
    nr=np.array(oldatm.nr.transpose([1,0,2]),order='F') # number concentration of rain drops
    qs=np.array(oldatm.qs.transpose([1,0,2]),order='F') # snow in air
    qi=np.array(oldatm.qi.transpose([1,0,2]),order='F') # cloud ice
    ni=np.array(oldatm.ni.transpose([1,0,2]),order='F') # number concentration of ice crystals
    qg=np.array(oldatm.qg.transpose([1,0,2]),order='F') # graupel in air
    pptrain=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptsnow=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptgraul=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptrainV=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptsnowV=np.zeros(N[0:2],order='F',dtype=np.float32)
    pptgraulV=np.zeros(N[0:2],order='F',dtype=np.float32)
    SR=np.zeros(N[0:2],order='F',dtype=np.float32)
    # zilch=np.zeros(N[0:2],order='F',dtype=np.float32)

    dz=np.array(np.gradient(atm.hgt)[0].transpose([1,0,2]),order="F",dtype=np.float32)
    # atm.p=get_newp(atm,z-zinit) #this needs to be done at the update step instead of in here. 
    # dzdx=np.gradient(z)
    # w=(U*dzdx[1]/dx + V*dzdx[0]/dx)/5
    # Fzs=fft.fftshift(fft.fft2(z))/(N[0]*N[1])
    # (u,v,w)=lop.winds(Fzs,atm)
    # U=oldatm.u.transpose([1,0,2])
    # V=oldatm.v.transpose([1,0,2])
    # W=oldatm.w.transpose([1,0,2])
    U=np.array(oldatm.u.transpose([1,0,2]),order='F') # vapor pressure
    V=np.array(oldatm.v.transpose([1,0,2]),order='F') # cloud water
    W=np.array(oldatm.w.transpose([1,0,2]),order='F') # rain in air
    newU=np.array(newu.transpose([1,0,2]),order='F')
    newV=np.array(newv.transpose([1,0,2]),order='F')
    newW=np.array(neww.transpose([1,0,2]),order='F')
    # thU=U.mean(axis=2).mean(axis=0)
    # thV=V.mean(axis=2).mean(axis=0)
    # thU2=newU.mean(axis=2).mean(axis=0)
    # thV2=newV.mean(axis=2).mean(axis=0)
    # thU=U.mean(axis=2).mean(axis=0)[np.newaxis,:].repeat(U.shape[0],axis=0)[:,:,np.newaxis].repeat(U.shape[2],axis=2)
    # thV=V.mean(axis=2).mean(axis=0)[np.newaxis,:].repeat(V.shape[0],axis=0)[:,:,np.newaxis].repeat(V.shape[2],axis=2)
    # thU2=newU.mean(axis=2).mean(axis=0)[np.newaxis,:].repeat(U.shape[0],axis=0)[:,:,np.newaxis].repeat(U.shape[2],axis=2)
    # thV2=newV.mean(axis=2).mean(axis=0)[np.newaxis,:].repeat(V.shape[0],axis=0)[:,:,np.newaxis].repeat(V.shape[2],axis=2)
    p=np.array(oldatm.p.transpose([1,0,2]), order='F')
    p2=np.array(atm.p.transpose([1,0,2]), order='F')
    
    t=0.0
    itimestep=0
    
    # calculate the maximum advection dt which still meets the Courant condition (here 0.5)
    dt=0.49/np.max(np.abs(U[1:,:,:])/dx+np.abs(V[:,:,1:])/dx)
    if dt>30:dt=30.0 # maximum dt = 30 seconds even if winds happen to be 0 so that the microphysics will be stable? 
    n=U.shape
    
    dth=np.array(((atm.th-oldatm.th)*dt/timestep).transpose([1,0,2]),order="F")
    dqv=np.array(((atm.qv-oldatm.qv)*dt/timestep).transpose([1,0,2]),order="F")
    du=(newU-U)*dt/timestep
    dv=(newV-V)*dt/timestep
    dw=(newW-W)*dt/timestep
    dp=np.array((p2-p)*dt/timestep,order="F")
    # set up domain parameters and assume a two layer model 
    ids=1;ide=N[0];jds=1;jde=N[1];kds=1;kde=n[1];
    ims=ids;ime=ide;jms=jds;jme=jde;kms=kds;kme=kde;
    # for now set it up to run the entire domain (except the outer border that is used for the boundary conditions)
    its=ids+1;ite=ide-1;jts=jds+1;jte=jde-1;kts=kds;kte=kde;
    # its=50;ite=350;jts=200;jte=400;kts=kds;kte=kde;
    # nextval=10
    verticalFirst=False
    # th=th*0+300
    # dth*=0
    # if fname:
    #     ncout.appendToVar(qv.transpose([1,0,2]),varname='qv',date=int(2*t/dt),pos=pos)
    #     # ncout.appendToVar(qc.transpose([1,0,2]),varname='qc',pos=pos)
    #     # ncout.appendToVar(qr.transpose([1,0,2]),varname='qr',pos=pos)
    #     # ncout.appendToVar(nr.transpose([1,0,2]),varname='nr',pos=pos)
    #     # ncout.appendToVar(qi.transpose([1,0,2]),varname='qi',pos=pos)
    #     # ncout.appendToVar(ni.transpose([1,0,2]),varname='ni',pos=pos)
    #     # ncout.appendToVar(qg.transpose([1,0,2]),varname='qg',pos=pos)
    #     # ncout.appendToVar(qs.transpose([1,0,2]),varname='qs',pos=pos)
    #     ncout.appendToVar(th.transpose([1,0,2]),varname='th',pos=pos)
    #     pos+=1
    
    while (t < timestep) :# & (dqc > 0.00001):
        # if ((t/timestep*100)>=nextval):
        #     print(np.round(t/timestep*100))
        #     nextval+=10
        # update the atmospheric grid p,u,v and boundary condition th,qv
        p=np.array(p+dp,order="F",dtype=np.float32)
        U+=du
        V+=dv
        W+=dw
        th=np.array(th+dth,order="F",dtype=np.float32)
        qv=np.array(qv+dqv,order="F",dtype=np.float32)
        # alpha=t/timestep
        # thU=thU1*(1-alpha)+thU2*alpha
        # thV=thV1*(1-alpha)+thV2*alpha
        # p=p1*(1-alpha)+p2*alpha
        pii=np.array(1.0/((100000.0/p)**(R/cp)),order="F",dtype=np.float32) # conversion from potential temperature to temperature
            
        # if fname:
        #     ncout.appendToVar(qv.transpose([1,0,2]),varname='qv',date=int(2*t/dt+1),pos=pos)
        #     ncout.appendToVar(qc.transpose([1,0,2]),varname='qc',pos=pos)
        #     ncout.appendToVar(qr.transpose([1,0,2]),varname='qr',pos=pos)
        #     ncout.appendToVar(nr.transpose([1,0,2]),varname='nr',pos=pos)
        #     ncout.appendToVar(qi.transpose([1,0,2]),varname='qi',pos=pos)
        #     ncout.appendToVar(ni.transpose([1,0,2]),varname='ni',pos=pos)
        #     ncout.appendToVar(qg.transpose([1,0,2]),varname='qg',pos=pos)
        #     ncout.appendToVar(qs.transpose([1,0,2]),varname='qs',pos=pos)
        #     ncout.appendToVar(th.transpose([1,0,2]),varname='th',pos=pos)
        #     pos+=1
        # print(dt, timestep, 100*t/timestep)
        # print(pptrain.max(),pptsnow.max(),qc.max(),qi.max(),qr.max(),qs.max(),U.max(),V.max())
        # print(pptrain.min(),pptsnow.min(),qc.min(),qi.min(),qr.min(),qs.min(),U.min(),V.min())
        if (t+dt)>timestep:
            dt=timestep-t
        # Run microphysics driver converts 3D to 1D and runs each column independantly
        # print("Microphysics...")
        # print(qv.dtype,qv.flags)
        # print(qc.dtype,qc.flags)
        # print(qr.dtype,qr.flags)
        # print(qi.dtype,qi.flags)
        # print(qs.dtype,qs.flags)
        # print(qg.dtype,qg.flags)
        # print(pii.dtype,pii.flags)
        # print(p.dtype,p.flags)
        # print(th.dtype,th.flags)
        # print(ni.dtype,ni.flags)
        # print(nr.dtype,nr.flags)
        # print(SR.dtype,SR.flags)
        # print(pptrain.dtype,pptrain.flags)
        # print(pptsnow.dtype,pptsnow.flags)
        # print(pptgraul.dtype,pptgraul.flags)
        
        if mp!=None:
            output=mp.module_mp_thompson.mp_gt_driver(qv, qc, qr, qi, qs, qg, ni, nr, 
                         th, pii, p, dz, dt, itimestep, 
                         pptrain, pptrainV, pptsnow, pptsnowV, pptgraul, pptgraulV, SR, 
                         int(ids),int(ide), int(jds),int(jde), int(kds),int(kde),              # domain dims
                         int(ims),int(ime), int(jms),int(jme), int(kms),int(kme),              # memory dims
                         int(its),int(ite), int(jts),int(jte), int(kts),int(kte)) 
        
        # perform advection for each species, eventually I should be able to remove graupel...
        # currently assumes U,V are terrain following (and they aren't too far off because NARR terrain is so smooth)
        # this may need to be cast as a single fortran call that loops over arrays and levels before returning, 
        # and eventually the whole time stepping scheme (advection+microphysics) should be a single fortran call. 
        # print("Advection...")
        # if fname:
        #     ncout.appendToVar(qv.transpose([1,0,2]),varname='qv',date=int(2*t/dt+1),pos=pos)
        #     # ncout.appendToVar(qc.transpose([1,0,2]),varname='qc',pos=pos)
        #     # ncout.appendToVar(qr.transpose([1,0,2]),varname='qr',pos=pos)
        #     # ncout.appendToVar(nr.transpose([1,0,2]),varname='nr',pos=pos)
        #     # ncout.appendToVar(qi.transpose([1,0,2]),varname='qi',pos=pos)
        #     # ncout.appendToVar(ni.transpose([1,0,2]),varname='ni',pos=pos)
        #     # ncout.appendToVar(qg.transpose([1,0,2]),varname='qg',pos=pos)
        #     # ncout.appendToVar(qs.transpose([1,0,2]),varname='qs',pos=pos)
        #     ncout.appendToVar(th.transpose([1,0,2]),varname='th',pos=pos)
        #     pos+=1
        qv[qv<1E-15]=1E-15
        qc[qc<1E-15]=1E-15
        qr[qr<1E-15]=1E-15
        nr[nr<1E-15]=1E-15
        qs[qs<1E-15]=1E-15
        qi[qi<1E-15]=1E-15
        ni[ni<1E-15]=1E-15
        qg[qg<1E-15]=1E-15
        
        # alternate vertical and horizontal advection order for stability? (there's a term for this)
        if verticalFirst:
            args=[(qv,W,dx,dt),(qc,W,dx,dt),(qr,W,dx,dt),(nr,W,dx,dt),(qs,W,dx,dt),
                  (qi,W,dx,dt),(ni,W,dx,dt),(qg,W,dx,dt),(th,W,dx,dt)]
            output=processPool.map(padvect,args)
            (qv,qc,qr,nr,qs,qi,ni,qg,th)=output
            args=[(qv[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])]
            args.extend([(qc[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(qr[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(nr[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(qs[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(qi[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(ni[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(qg[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            args.extend([(th[:,j,:],U[:,j,:],V[:,j,:],dx,dt) for j in range(U.shape[1])])
            # args=[(qv[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qv[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qv[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qv[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (qc[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qc[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qc[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qc[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (qr[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qr[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qr[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qr[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (nr[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(nr[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(nr[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(nr[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (qs[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qs[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qs[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qs[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (qi[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qi[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qi[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qi[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (ni[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(ni[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(ni[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(ni[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (qg[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qg[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qg[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qg[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
            #       (th[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(th[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(th[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(th[:,3,:],U[:,3,:],V[:,3,:],dx,dt)]
            output=processPool.map(padvect,args)
            k=0
            for j in range(U.shape[1]):
                qv[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                qc[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                qr[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                nr[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                qs[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                qi[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                ni[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                qg[:,j,:]=output[k]
                k+=1
            for j in range(U.shape[1]):
                th[:,j,:]=output[k]
                k+=1
            # (qv[:,0,:],qv[:,1,:],qv[:,2,:],qv[:,3,:],qc[:,0,:],qc[:,1,:],qc[:,2,:],qc[:,3,:],qr[:,0,:],qr[:,1,:],qr[:,2,:],qr[:,3,:],
            #  nr[:,0,:],nr[:,1,:],nr[:,2,:],nr[:,3,:],qs[:,0,:],qs[:,1,:],qs[:,2,:],qs[:,3,:],qi[:,0,:],qi[:,1,:],qi[:,2,:],qi[:,3,:],
            #  ni[:,0,:],ni[:,1,:],ni[:,2,:],ni[:,3,:],qg[:,0,:],qg[:,1,:],qg[:,2,:],qg[:,3,:],th[:,0,:],th[:,1,:],th[:,2,:],th[:,3,:])=output
            verticalFirst=False
        else:
            args=[(qv[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qv[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qv[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qv[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (qc[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qc[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qc[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qc[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (qr[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qr[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qr[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qr[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (nr[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(nr[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(nr[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(nr[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (qs[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qs[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qs[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qs[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (qi[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qi[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qi[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qi[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (ni[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(ni[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(ni[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(ni[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (qg[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(qg[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(qg[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(qg[:,3,:],U[:,3,:],V[:,3,:],dx,dt),
                  (th[:,0,:],U[:,0,:],V[:,0,:],dx,dt),(th[:,1,:],U[:,1,:],V[:,1,:],dx,dt),(th[:,2,:],U[:,2,:],V[:,2,:],dx,dt),(th[:,3,:],U[:,3,:],V[:,3,:],dx,dt)]
            output=processPool.map(padvect,args)
            (qv[:,0,:],qv[:,1,:],qv[:,2,:],qv[:,3,:],qc[:,0,:],qc[:,1,:],qc[:,2,:],qc[:,3,:],qr[:,0,:],qr[:,1,:],qr[:,2,:],qr[:,3,:],
             nr[:,0,:],nr[:,1,:],nr[:,2,:],nr[:,3,:],qs[:,0,:],qs[:,1,:],qs[:,2,:],qs[:,3,:],qi[:,0,:],qi[:,1,:],qi[:,2,:],qi[:,3,:],
             ni[:,0,:],ni[:,1,:],ni[:,2,:],ni[:,3,:],qg[:,0,:],qg[:,1,:],qg[:,2,:],qg[:,3,:],th[:,0,:],th[:,1,:],th[:,2,:],th[:,3,:])=output
            args=[(qv,W,dx,dt),(qc,W,dx,dt),(qr,W,dx,dt),(nr,W,dx,dt),(qs,W,dx,dt),
                  (qi,W,dx,dt),(ni,W,dx,dt),(qg,W,dx,dt),(th,W,dx,dt)]
            (qv,qc,qr,nr,qs,qi,ni,qg,th)=processPool.map(padvect,args)
            verticalFirst=True
        # ideally this should just be ensured by working entirely with F ordered arrays internally so that we don't have this overhead
        # translate the entire loop to fortran and it won't matter? 
        qv=np.array(qv,order="F",dtype=np.float32)
        qc=np.array(qc,order="F",dtype=np.float32)
        qr=np.array(qr,order="F",dtype=np.float32)
        nr=np.array(nr,order="F",dtype=np.float32)
        qs=np.array(qs,order="F",dtype=np.float32)
        qi=np.array(qi,order="F",dtype=np.float32)
        ni=np.array(ni,order="F",dtype=np.float32)
        qg=np.array(qg,order="F",dtype=np.float32)
        th=np.array(th,order="F",dtype=np.float32)
        
        # for i in range(n[1]):
            # args=[(qv[:,i,:],thU[i],thV[i],dx,dt),
            #       (qc[:,i,:],thU[i],thV[i],dx,dt),
            #       (qr[:,i,:],thU[i],thV[i],dx,dt),
            #       (nr[:,i,:],thU[i],thV[i],dx,dt),
            #       (qs[:,i,:],thU[i],thV[i],dx,dt),
            #       (qi[:,i,:],thU[i],thV[i],dx,dt),
            #       (ni[:,i,:],thU[i],thV[i],dx,dt),
            #       (qg[:,i,:],thU[i],thV[i],dx,dt),
            #       (th[:,i,:],thU[i],thV[i],dx,dt)]
            # processorPool.map(padvect,args)
            # print("  advecting level "+str(i))
            # qv[qv<1E-15]=1E-15
            # qv[:,i,:]=advect.advect(qv[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # qc[qc<1E-15]=1E-15
            # qc[:,i,:]=advect.advect(qc[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # qr[qr<1E-15]=1E-15
            # qr[:,i,:]=advect.advect(qr[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # nr[nr<1E-15]=1E-15
            # nr[:,i,:]=advect.advect(nr[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # qs[qs<1E-15]=1E-15
            # qs[:,i,:]=advect.advect(qs[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # qi[qi<1E-15]=1E-15
            # qi[:,i,:]=advect.advect(qi[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # ni[ni<1E-15]=1E-15
            # ni[:,i,:]=advect.advect(ni[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # qg[qg<1E-15]=1E-15
            # qg[:,i,:]=advect.advect(qg[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
            # th[:,i,:]=advect.advect(th[:,i,:],thU[:,i,:],thV[:,i,:],dx,dt,courant=True)
        # plt.clf();plt.imshow(qv[0:450,0,100:500],vmin=0,vmax=0.003);plt.colorbar();plt.title(str(t))
        # plt.draw()
        # th+=dTdt*dt
        t+=dt
    # if fname:
    #     ncout.close()
    atm.th=th.transpose([1,0,2])
    atm.qv=qv.transpose([1,0,2])
    atm.qc=qc.transpose([1,0,2])
    atm.qr=qr.transpose([1,0,2])
    atm.nr=nr.transpose([1,0,2])
    atm.qi=qi.transpose([1,0,2])
    atm.ni=ni.transpose([1,0,2])
    atm.qs=qs.transpose([1,0,2])
    atm.qg=qg.transpose([1,0,2])
    atm.u=newu
    atm.v=newv
    atm.p=newp
    return pptrain+pptgraul+pptsnow
        