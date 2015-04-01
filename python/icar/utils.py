import numpy as np
import units
def exner(th,p):
    Rd=287.058
    cp=1004.0
    p0=100000
    pii=(p/p0)**(Rd/cp)
    return th * pii
    
def q2rh(qv,t,p):
    """docstring for q2rh"""
    return units.mr2rh(t,p,mr=qv)
    
