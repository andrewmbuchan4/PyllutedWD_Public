#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numba import jit

def S_disc(t_formation):
    #define the constants
    Mstar = 2.34
    s0 = (33*1.496*10**11) 
    sigma = (5.67*10**(-8)) 
    k = (1.38*10**(-23))
    mH = (1.67*10**(-27))
    alpha = 0.01
    gamma = 1.7
    mu = 2.4
    k0 = 0.3
    sun = (2.0e30)
    G = (6.674*10**(-11))
    #stellar evolution tracks from Siess 2000
    Tstar = 259.5822261222*(Mstar**3)-1424.0943845798*(Mstar**2)+2808.5873546861*(Mstar)+2644.9291158227
    Rstar = (6.955*10**8)*(0.1543971235*( Mstar**8) - 1.2349012901*(Mstar**7)+2.7854995*(Mstar**6)+3.0573592544*(Mstar**5)-25.2247969873*(Mstar**4)+47.9947080829*(Mstar**3)-42.5877850958*(Mstar**2)+18.4521562096*(Mstar)-0.3995121438)
    #disc model constants from Chambers 2009
    Tvis = (((27*k0)/(64*sigma))**(1/3))*(((alpha*gamma*k)/(mu*mH))**(1/3))*(((7*sun)/(100*3.14*(s0**2)))**(2/3))*(((G*sun)/(s0**3))**(1/6))*(Mstar**(5/6))
    Trad = ((4/7)**(1/4))*(((((Tstar)**8)*((Rstar)**(4))*k)/(G*sun*mu*mH*(s0**(3))))**(1/7))*(Mstar**(-1/7))
    Svis = (7*sun)/(100*3.14*(s0**2))*(Mstar) 
    Srad= Svis*(Tvis/Trad)
    tvis = (1/(16*3.14))*((mu*mH)/(alpha*gamma*k))*(sun/10)*(((G*sun)/((s0)**3))**(0.5))*(1/(Svis))*(1/(Tvis))*(Mstar**(3/2)) 
    invtvis = 1/tvis        
    t1 = tvis*(((Tvis/Trad)**(112/73))-1)
    r1 = s0*((Srad/Svis)**(70/33))*((1+t1*invtvis)**(-133/132))
    T1 = Trad*(r1/s0)**(-3/7)
    S1 = Srad*((r1/s0)**(-15/14))*((1+t1*invtvis)**(-19/16));
    Mdot1 = ((0.3*sun*Mstar*invtvis)/16)*(T1/Tvis)*(S1/Svis)*((r1/s0)**1.5)
    trad = (0.7*Mstar*sun*((Trad/Tvis)**(21/73)))/(13*Mdot1)
    invtrad = 1/trad;
    #Chambers Disc Model
    if 0 <= t_formation*(3.1536*10**13) <= t1:
        s = (s0*(1+(t_formation*(3.1536*10**13)*invtvis))**(6/16))            
    elif t_formation*(3.1536*10**13) > t1:
        s = (s0*((Tvis/Trad)**(42/73))*((1+((t_formation*(3.1536*10**13)-t1)*invtrad))**(14/13)))
    else:
        s = s0
    return s/(1.496*10**11)

@jit(nopython=True)    
def T_disc(d_formation, t_formation):
    #define the constants
    Mstar = 2.34
    s0 = (33*1.496*10**11) 
    sigma = (5.67*10**(-8)) 
    k = (1.38*10**(-23))
    mH = (1.67*10**(-27))
    alpha = 0.01
    gamma = 1.7
    mu = 2.4
    k0 = 0.3
    sun = (2.0e30)
    G = (6.674*10**(-11))
    #stellar evolution tracks from Siess 2000
    Tstar = 259.5822261222*(Mstar**3)-1424.0943845798*(Mstar**2)+2808.5873546861*(Mstar)+2644.9291158227
    Rstar = (6.955*10**8)*(0.1543971235*( Mstar**8) - 1.2349012901*(Mstar**7)+2.7854995*(Mstar**6)+3.0573592544*(Mstar**5)-25.2247969873*(Mstar**4)+47.9947080829*(Mstar**3)-42.5877850958*(Mstar**2)+18.4521562096*(Mstar)-0.3995121438)
    #disc model constants from Chambers 2009
    Te = 1380
    Tvis = (((27*k0)/(64*sigma))**(1/3))*(((alpha*gamma*k)/(mu*mH))**(1/3))*(((7*sun)/(100*3.14*(s0**2)))**(2/3))*(((G*sun)/(s0**3))**(1/6))*(Mstar**(5/6))
    Trad = ((4/7)**(1/4))*(((((Tstar)**8)*((Rstar)**(4))*k)/(G*sun*mu*mH*(s0**(3))))**(1/7))*(Mstar**(-1/7))
    Svis = (7*sun)/(100*3.14*(s0**2))*(Mstar) 
    Se = Svis*((Tvis/Te)**(14/19)) 
    Srad= Svis*(Tvis/Trad)
    tvis = (1/(16*3.14))*((mu*mH)/(alpha*gamma*k))*(sun/10)*(((G*sun)/((s0)**3))**(0.5))*(1/(Svis))*(1/(Tvis))*(Mstar**(3/2)) 
    invtvis = 1/tvis
    #Chambers Disc Model
    if 0 < d_formation < (((s0)*((Se/Svis)**(95/63))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-19/36)))/(1.496*10**11)):
        T = (Tvis**(5/19))*(Te**(14/19))*((((1.496*10**11)*d_formation)/s0)**(-9/38))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-1/8)) 
    elif (((s0)*((Se/Svis)**(95/63))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-19/36)))/(1.496*10**11)) < d_formation < (((s0)*((Srad/Svis)**(70/33))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-133/132)))/(1.496*10**11)):
        T = Tvis*((((1.496*10**11)*d_formation)/s0)**(-9/10))*((1+((t_formation*(3.1536*10**13))/(tvis)))**(-19/40))
    elif d_formation > (((s0)*((Srad/Svis)**(70/33))*((1+(t_formation*(3.1536*10**13)*(invtvis)))**(-133/132)))/(1.496*10**11)):
        T = Trad*(((((1.496*10**11)*d_formation)/s0)**(-3/7))) 
    else:
        T = 5000
    return T
