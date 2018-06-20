#!/usr/bin/env python
# encoding: utf-8
"""
units.py

converts between various formulations of atmospheric humidity
    specific humidity (sh)
    relative humidity (rh)
    dewpoint (dp)

    calculate saturated vapor pressure from a temperature

    xy2a = convert x,y point (or u/v vector) to the corresponding angle (0,0) to (x,y)

    Various utilities to calculate elevation as a function of pressure, temperature, humidity, and sea level pressure

    for atmospheric purposes can compute lifting condensation level (LCL)
    and the temperature at the LCL

    could also explicitly add
    vapor pressure (vp)

Created by Ethan Gutmann on 2011-08-15.
Copyright (c) 2011 National Center for Atmospheric Research. All rights reserved.
"""
import numpy as np

def rtod(angle):
    return angle*360.0/(2*np.pi)

def dtor(angle):
    return angle/360.0*(2*np.pi)

def moist_BV_frequency_squared(t,z,p,ql):
    # return None
    g = 9.81 #m/s^2
    R = 287.058 # J/kg/K
    # L = 2260.0*1000 # J/kg
    # L = 2.5e6 # J/kg
    L = (2500.8 - 2.36*(t-273.15) + 0.0016*(t-273.15)**2 - 0.00006*(t-273.15)**3)*1000 # J/kg

    qs = rh2mr(t,p,1.0)
    gamma_moist = sat_lapse_rate(t,mr=qs)
    print("Gmoist=",gamma_moist[0])
    print("qsat=  ",qs[0])

    dtdz = np.zeros(t.shape)
    dtdz[:-1] = (t[1:]-t[:-1]) / (z[1:]-z[:-1])
    dtdz[-1]=dtdz[-2]
    print("dtdz=  ",dtdz[0])

    qw=qs+ql
    dqwdz = np.zeros(qw.shape)
    dqwdz[:-1] = (qw[1:]-qw[:-1]) / (z[1:]-z[:-1])
    dqwdz[-1] = dqwdz[-2]
    print("dqwdz= ",dqwdz[0])
    print("t =    ",t[0])
    print(L/R)
    bv = g/t * (dtdz + gamma_moist) * (1 + (L*qs)/(R*t)) - (g/(1+qw))*dqwdz
    return bv

def dry_BV_frequency_squared(theta,z,bottom=0,top=None,t=None,p=None):
    g=9.81

    if (theta==None):
        if t is not None and p is not None:
            theta=t/exner(p)
        else:
            raise ValueError("Needs either potential temperature[K] or real T[K] and P[Pa]")

    lntheta=np.log(theta)
    if top==None:
        bv=g * (lntheta[bottom+1:top]-lntheta[bottom:-1]) / (z[bottom+1:top] - z[bottom:-1])
    else:
        bv=g * (lntheta[bottom+1:top]-lntheta[bottom:top-1]) / (z[bottom+1:top] - z[bottom:top-1])

    return bv

def exner(p):
    """use this to get the exner function
    exner * potential temperature = real temperature
    """
    Rd=287.058
    cp=1003.5
    p0=1000.0
    try:
        if p.max()>1200:
            p0*=100.0
    except:
        if p>1200:
            p0*=100.0

    return (p/p0)**(Rd/cp)

# def p2z(p,t,p0=101325.0):
#     '''Calculate elevation [m] for a given pressure [Pa] and column mean temperature [K],
#         and optionally sea level pressure (p0 [Pa])
#
#         Based off Babinet's Formula (doesn't seem to work.)
#         '''
#     z=np.zeros(p.shape)
#
#     z[0,...]=(16000.0+64.0*t[0,...])*(p0-p[0,...])/(p0+p[0,...]) #Babinet's Formula
#     for i in range(1,z.shape[0]):
#         p0=p[i-1,...]
#         z[i,...]=(16000.0+64.0*t[i,...])*(p0-p[i,...])/(p0+p[i,...])+z[i-1,...] #Babinet's Formula
#
#     return z

# return the moist / saturated adiabatic lapse rate for a given
# Temperature and mixing ratio (really MR could be calculated as f(T))
# from http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
def sat_lapse_rate(T,mr=None, p=None):
    L  = 2.5e6  # J /kg
    g  = 9.81   # m/s^2
    Rd = 287.0    # J /kg /K
    Rw = 461.5  # J /kg /K
    ratio=Rd/Rw # []
    cp = 1005.0 # J /kg /K

    # T [=] K
    # mr[=] kg/kg

    if mr==None:
        mr = rh2mr(T,p,100.0)

    return g* ((1 + ((L*mr) / (Rd*T)) )
            / (cp + (((L**2)*mr*ratio) / (Rd*(T**2))) ))


def calc_tv(t,mr=None,e=None,p=None):
    '''calculate the virtual temperature

    Uses a real temperature and mixing ratio or (vapor pressure and barometric pressure)

    t [=] K
    mr[=] kg/kg
    e [=] any (same as p) e.g. Pa, hPa
    p [=] any (same as e)
    '''
    # from http://en.wikipedia.org/wiki/Virtual_temperature
    RoR = 0.622 # Rdry / Rvapor

    # if mixing ratio was specified:
    if mr is not None:
        return t*(mr+RoR)/(RoR*(1+mr))
    # if vapor pressure and barometric pressure were given:
    if (e is not None) and (p is not None):
        T/(1-e/p*(1-RoR))

def calc_z(slp,p,t,mr, timeseries=False):
    """Calculate geopotential height [m] from a 3D atm field and sea level pressure
    the 3D fields must include pressure, temperature, and humidity.

    slp = sea level pressure    [Pa]
    p   = 3D pressure field     [Pa]
    ta  = 3D temperature field  [K]
    hus = 3D specific humidity  [kg/kg]
    mr  = 3D mixing ratio       [kg/kg]

    returns z = 3D geopotential height field [m]

    based off WMO CIMO guide: http://www.wmo.int/pages/prog/www/IMOP/CIMO-Guide.html
                    Part I Chapter 3
    """
    tv = calc_tv(t, mr=mr)
    Kp = 0.0148275 # K / gpm
    # if these are likely 2D data, not 3D just do a simple calculation
    if len(tv.shape) <= 2:
        return tv * np.log10(slp / p) / Kp

    # for 3D data integrate up from the surface
    z = np.zeros(tv.shape)
    # do the simple calculation at the bottom
    if timeseries:
        z[:,0] = tv[:,0] * np.log10(slp / p[:,0]) / Kp
        # now compute the elevation of each successive layer
        for i in range(1,z.shape[1]):
            z[:,i] = (tv[:,i] + tv[:,i-1]) / 2  * np.log10(p[:,i-1] / p[:,i]) / Kp + z[:,i-1]

    else:
        z[0] = tv[0] * np.log10(slp / p[0]) / Kp
        # now compute the elevation of each successive layer
        for i in range(1,z.shape[0]):
            z[i] = (tv[i] + tv[i-1]) / 2  * np.log10(p[i-1] / p[i]) / Kp + z[i-1]

    return z


def z2p(p,h):
    '''Convert p [Pa] at elevation h [m] by shifting its elevation by dz [m]'''
    # p    in pascals (or hPa or...)
    # h    in meters
    return p*(1 - 2.25577E-5*h)**5.25588

def zt2p(z,p0=101325.0,t0=288.15,dtdz= -0.0065,zaxis=1,use_z_axis=None):
    # p0=101325    #Pa
    # t0=288.15   #K
    # dtdz= -0.0065 #K/m
    g=9.807
    M=0.029
    R=8.314
    if (not isinstance(t0,np.ndarray)):
        p= p0*(t0/(t0+dtdz*z))**((g*M)/(R*dtdz))
        return p

    if use_z_axis==None:
        use_z_axis=len(t0.shape)>3

    if use_z_axis:
        p=np.zeros(t0.shape)
        for i in range(t0.shape[zaxis]):
            if zaxis==0:
                p[i]=p0*(t0.take(i,axis=zaxis)/(t0.take(i,axis=zaxis)+dtdz*z.take(i,axis=zaxis)))**((g*M)/(R*dtdz))
            elif zaxis==1:
                p[:,i]=p0*(t0.take(i,axis=zaxis)/(t0.take(i,axis=zaxis)+dtdz*z.take(i,axis=zaxis)))**((g*M)/(R*dtdz))
            elif zaxis==2:
                p[:,:,i]=p0*(t0.take(i,axis=zaxis)/(t0.take(i,axis=zaxis)+dtdz*z.take(i,axis=zaxis)))**((g*M)/(R*dtdz))
            elif zaxis==3:
                p[:,:,i]=p0*(t0.take(i,axis=zaxis)/(t0.take(i,axis=zaxis)+dtdz*z.take(i,axis=zaxis)))**((g*M)/(R*dtdz))
            else:
                raise ValueError("not setup to process higher axes")



            # p.take(i,axis=zaxis)=p0*(t0.take(i,axis=zaxis)/(t0.take(i,axis=zaxis)+dtdz*z.take(i,axis=zaxis)))**((g*M)/(R*dtdz))
    else:
        p= p0*(t0/(t0+dtdz*z))**((g*M)/(R*dtdz))
            # p[:,i,...]=p0*(t0[:,i,...]/(t0[:,i,...]+dtdz*z[:,i,...]))**((g*M)/(R*dtdz))

    return p

# def calc_z(slp,p,t,hus):
#     """Calculate geopotential height [m] from a 3D atm field and sea level pressure
#     the 3D fields must include pressure, temperature, and humidity.
#
#     slp = sea level pressure    [Pa]
#     p   = 3D pressure field     [Pa]
#     ta  = 3D temperature field  [K]
#     hus = 3D specific humidity  [kg/kg]
#
#     returns z = 3D geopotential height field [m]
#
#     Based off equation found on http://www.meteormetrics.com/correctiontosealevel.htm
#     Which states that it comes from:
#         WMO tech note 91, Methods in use for the reduction of atmospheric pressure
#
#     """
#
#     z0 = p2z(p,t,p0=slp) # rough first guess
#     z=np.zeros(p.shape)
#
#     dtdz=np.diff(t,axis=0)/np.diff(z0,axis=0)
#     K = 18400.0 # barometric / hypsometric constant
#     alpha = 0.0037 # thermal expansion coefficient for air
#     b = (p[0,...]+slp)/2 #mean barometric pressure of the air column
#     asphericity=1.00045 # = correction for asphericity of the earth at 40 degrees latitude
#     # asphericity = 1/(1-0.0026*cos(2*latitude))
#     mr=hus/(1-hus)
#     e=mr*p/(0.62197+mr)
#
#     z[0] = K * (1+alpha*t[0]+dtdz[0]*z0[0]/2) * 1.0/(1-0.378*e[0]/b) * asphericity * np.log10(slp/p[0])
#     for i in range(1,z.shape[0]):
#         b=(p[i]+p[i-1])/2
#         z[i] = K * (1+alpha*(t[i]+t[i-1])/2) * 1.0/(1-0.378*(e[i-1]+e[i])/2/b) * asphericity * np.log10(p[i-1]/p[i])+z[i-1]
#
#     return z,z0,e,dtdz




def calc_slp(ps,z,ts=288.15,dtdz= -0.0065,mr=None,e=None,sh=None,latitude=None,method=1):
    """ Calculate sea level pressure from as much information as available

    based off the WMO handbook:
    http://www.wmo.int/pages/prog/www/IMOP/meetings/SI/ET-Stand-1/Doc-10_Pressure-red.pdf
    excerpt from CIMO Guide, Part I, Chapter 3 (Edition 2008, Updated in 2010)
    equation 3.2
    see also: http://www.meteormetrics.com/correctiontosealevel.htm

    ps   = station pressure [Pa or hPa, output slp will have the same units]
    z    = station elevation [m]
    ts   = station temperature [K or C]   OPTIONAL
    dtdz = assumed fictitious temperature lapse rate down to sea level [K/m]   OPTIONAL
    mr   = station water vapor mixing ratio [kg/kg]   OPTIONAL
    e    = station water vapor pressure [Pa or hPa (consistent with ps)]   OPTIONAL
    sh   = station specific humidity [kg/kg]   OPTIONAL
    latitude=station latitude [deg] OPTIONAL
    """
    # p0=101325    #Pa
    # t0=288.15   #K
    # dtdz= -0.0065 #K/m
    g=9.80665
    M=0.029
    R=287.05

    #############################################
    #
    # This section is for converting units
    #
    p_inhPa=False
    try:
        if ps.max()<1200:
            ps*=100
            p_inhPa=True
    except:
        if ps<1200: # assume ps is in Pa, not hPa
            ps*=100
            p_inhPa=True

    if sh is not None:
        mr=sh/(1-sh)
    if mr is not None:
        #assumes mr is in kg/kg
        e=mr*ps/(0.62197+mr)

    if e is None:
        mr=0.01
        e=mr*ps/(0.62197+mr)

    T_inC=False
    try:
        if ts.min()<100:
            ts+=273.15
            T_inC=True
    except:
        if ts<100:
            ts+=273.15
            T_inC=True
    #
    # End of unit convertions
    #
    #############################################

    Ch = 0.0012 # K/Pa

    Hp=z # geopotential height
    a= dtdz # K/gpm

    # if latitude != None:
    #     # N.B. this may not be correct... removing for now
    #     k=0.0026 # earth shape factor
    #     lat_factor=1/(1-k*np.cos(2*dtor(latitude)))
    # else:
    #     lat_factor=1

    # convert p back to hPa
    if p_inhPa:
        ps/=100

    # these both come from CIMO Guide, and one would think they would be identical...
    if method==1:
        slp= ps*np.exp(((g/R)*Hp) / (ts - a*Hp/2.0 + e*Ch))
    elif method==2:
        Kp = 0.0148275 # K / gpm
        # based off the virtual temperature instead
        if mr is not None:
            tv=calc_tv(ts,mr=mr)
        elif e is not None:
            tv=calc_tv(ts,e=e,p=ps)
        slp= ps*10**(Kp*Hp/tv)

    # convert t back into degrees C
    if T_inC:
        t-=273.15

    return slp


def xy2a(x,y):
    '''
    convert an x,y coordinate to the angle to that coordinate (0-360)
    '''

    angle=rtod(np.arctan(np.cast['f'](x)/y))
    if np.array(x).size==1:
        if x>0:
            if y>0:
                pass
            else:
                angle=180+angle
        else:
            if y>0:
                angle=360+angle
            else:
                angle=180+angle

    else:
        for i in range(len(x)):
            if x[i]>0:
                if y[i]>0:
                    pass
                else:
                    angle[i]=180+angle[i]
            else:
                if y[i]>0:
                    angle[i]=360+angle[i]
                else:
                    angle[i]=180+angle[i]


    return angle


def t2vp(t):
    '''Saturated vapor pressure for a given temperature
    T in degrees C or K'''
    if np.array(t).max()>150:
        tIsInK=True
    else:
        t+=273.15
        tIsInK=False
    freezing=273.15

    if isinstance(t,np.ndarray):
        a=np.zeros(t.shape)+17.2693882
        b=np.zeros(t.shape)+35.86
        a[np.where(t<freezing)]=21.8745584
        b[np.where(t<freezing)]=7.66
    else:
        if (t<freezing):
            a=21.8745584
            b=7.66
        else:
            a=17.2693882
            b=35.86

    vp = 6.1078* np.exp(a*(t-273.16)/(t-b)) #(hPa)
    # vp=6.112*np.exp(17.67*t/(t+243.5)) #*10**((7.5*t)/(237.7+t))
    if not tIsInK:
        t-=273.15
    return vp

def sh2rh(t,p,sh):
    '''T(deg.C or K), p(mb), Specific Humidity(kg/kg)'''
    mr=sh/(1-sh)
    return mr2rh(t,p,mr)

def rh2sh(t,p,rh):
    '''T(deg.C or K), p(mb), Relative Humidity(0-1 or 0-100)'''
    mr=rh2mr(t,p,rh)
    return mr/(1+mr)

def sh2dp(t,p,sh):
    '''T(deg.C or K), p(mb), SH(kg/kg)'''
    return rh2dp(t,sh2rh(t,p,sh))

def mr2rh(t,p,mr):
    '''T(deg.C or K), p(mb), MR(kg/kg)'''
    # e=mr*p/(621.97+mr)
    e=mr*p/(0.62197+mr)
    rh=e/t2vp(t)
    return rh

def rh2mr(t,p,rh):
    '''T(deg.C or K), p(mb), rh(0-1 or 0-100 assumed if max is over 1.5)'''
    if np.array(rh).max()>1.5:
        multiplier=1/100.0
    else:
        multiplier=1.0

    e=t2vp(t)*rh*multiplier
    mr=0.62197*e/(p-e)

    return mr

def dp2mr(t,p,dp):
    '''T(deg.C or K), p(mb), dp(deg.C or K)'''
    return rh2mr(t,p,dp2rh(t,dp))

def mr2dp(t,p,mr):
    '''T(deg.C or K), p(mb), mr(kg/kg)'''
    if np.array(t).max()>150:
        tIsInK=True
        t-=273.15
    else:
        tIsInK=False
    dp=rh2dp(t,mr2rh(t,p,mr))
    if tIsInK:
        t+=273.15
        dp+=273.15
    return dp

def rh2dp(t,rh):
    '''T in deg.C or K RH:0-1 or 0-100 if RH>1'''
    if np.array(rh).max()>1:rh/=100.0
    if np.array(t).max()>150:
        tIsInK=True
        t-=273.15
    else:
        tIsInK=False
    e=t2vp(t) * rh
    dp=243.5*np.log(e/6.112)/(17.67-np.log(e/6.112))
    if tIsInK:
        t+=273.15
        dp+=273.15
    return dp

def dp2rh(t,dp):
    '''T and DP in deg.C or K (T and DP can be in different units)'''
    e=t2vp(dp)
    esat=t2vp(t)
    return np.max((np.min((0,e/esat)),1))

def find_lcl(p,t,dp):
    '''Find the lifting condensation level (m) p=[Pa,hPa],t[C,K],dp[C,K]

    Units are permitted to vary,
        p can be hPa, or Pa
        t and dp can be C or K
            but both t and dp must be the same (C or K)

    Main code from Greg Thompson's FORTRAN77 (meteo_funcs.f) based on Bolton (1980)?
    '''
    if np.array(t).max()<150:
        tIsInC=True
        t+=273.15
        dp+=273.15
    else:
        tIsInC=False
    if np.array(p).max()<1500:
        pIsInHPa=True
        p*=100.0
    else:
        pIsInHPa=False

    rcp=287.04/1004.0
    cpr=1.0/rcp

    tlcl=t_lcl(t,dp)
    theta=t*(100000.0/p)**rcp
    plcl = 100000.0 * (tlcl/theta)**cpr
    lcl=(1-((plcl/101325.0)**(1.0/5.25588)))/2.25577E-5

    if pIsInHPa:
        p/=100
    if tIsInC:
        t-=273.15
        dp-=273.15
    return lcl


def t_lcl(t,dp):
    '''Convert T and dew point to T at lifting condensation level

    Units can be in C or K, if T<150 assumes C and returns C
    if T>=150 assumes K and returns K
    Code from Greg Thompson FORTRAN77 (meteo_funcs.f) based on Bolton (1980) eqn #15
    claims accuracy of 0.1K for -35<T<35 C
    '''
    tIsInC=False
    if np.array(t).max()<150:
        tIsInC=True
        t+=273.15
        dp+=273.15
    t_lcl=(1.0/(1.0/((dp-56.0))+(np.log(t/dp)/800.))) + 56.0
    if tIsInC:
        t-=273.15
        dp-=273.15
        t_lcl-=273.15
    return t_lcl
