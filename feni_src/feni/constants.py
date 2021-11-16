# constants.py

import numpy as np
import types

from GenData import GenData as GenData

# non-bridging oxygens over tetrahedral cations
nbot = 2.7

# fixed activity coefficient of Fe in silicate (Rudge 2011)
gammaFe_sil = 3.0

# ---------- atomic masses -------------
M = GenData()
M.al = 26.9815385
M.c = 12.011
M.ca = 40.078
M.cr = 51.9961
M.fe = 55.845
M.k = 39.0983
M.mg = 24.305
M.mn = 54.938044
M.na = 22.98976928
M.n = 14.007
M.ni = 58.6934
M.o = 15.9994
M.p = 30.973761998
M.s = 32.06
M.si = 28.085
M.ti = 47.867

# ---------- molecular masses ----------
MO = GenData()
MO.si = M.si + M.o*2
MO.ti = M.ti + M.o*2
MO.al = M.al*2 + M.o*3
MO.fe = M.fe + M.o
MO.mn = M.mn + M.o
MO.mg = M.mg + M.o
MO.ca = M.ca + M.o
MO.na = M.na*2 + M.o
MO.k = M.k*2 + M.o
MO.p = M.p*2 + M.o*5
MO.c = M.c + M.o*2
MO.cr = M.cr*3 + M.o*2
MO.ni = M.ni + M.o

# ---------- chondrite abundances ----------
# lodders 2003
# atomic fraction
class CH(GenData):

    logal = 6.46-12   ; logalE = 0.02
    logc = 7.43-12    ; logcE = 0.06
    logca = 6.32-12   ; logcaE = 0.03
    logcr = 5.66-12   ; logcrE = 0.05
    logfe = 7.48-12   ; logfeE = 0.03
    logk = 5.09-12    ; logkE = 0.05
    logmg = 7.56-12   ; logmgE = 0.02
    logmn = 5.50-12   ; logmnE = 0.03
    logna = 6.30-12   ; lognaE = 0.03
    logn = 6.28-12    ; lognE = 0.07
    logni = 6.22-12   ; logniE = 0.03
    logo = 8.42-12    ; logoE = 0.02
    logp = 5.43-12    ; logpE = 0.04
    logs = 7.19-12    ; logsE = 0.04
    logsi = 7.54-12   ; logsiE = 0.02
    logti = 4.92-12   ; logtiE = 0.03

    chems = ['al', 'c', 'ca', 'cr', 'fe', 'k', 'mg', 'mn', 'na', 'n', 'ni', 'o', 'p', 's', 'si', 'ti']

    def log2x( self, eles=None ):
        if eles is None:
            eles = self.chems
            
        ss = 0.
        for x in eles:
            if not np.isnan(getattr(self, "log"+x)) :
                ss += np.power(10,getattr(self, "log"+x))

        for x in eles:
            setattr(self, x,     np.power(10,getattr(self, "log"+x))/ss )
            setattr(self, x+"E",     [] )

        # generate errors
        # run MC simulation propagating gaussian distributed errors on the log(X/H) terms
        #   through to variability in the atomic fraction terms: Xi's
        # number of MC iterations
        Ni = 1000
        c_tmp = np.zeros((Ni,len(eles)))
        c_tmp[:] = np.nan
        for k in range(0,Ni):
            l = 0
            for x in eles:
                c = getattr(self, "log"+x)
                ce = getattr(self, "log"+x+"E")
                if not np.isnan(c) :
                    if not np.isnan(ce):
                        c_tmp[k,l] = np.power(10,np.random.normal(loc=c,scale=ce))
                        l += 1
                    
        ss = np.nansum(c_tmp, axis=1)
        c_tmp = c_tmp/ss[:,None]
        ss_tmp = np.nanstd(c_tmp, axis=0)
    
        l=0
        for x in eles:
            getattr(self, x+"E").append(ss_tmp[l])
            l += 1


    def concs(self, MO=None, eles=None):

        if MO is not None :
            def x2c(self,MO,eles):
                sum = 0.
                for e in eles:
                    setattr(self, "C"+e,getattr(self, e)*getattr(M,e))
                    sum += getattr(self, e)*getattr(M,e)
                for e in eles:
                    setattr(self, "C"+e, getattr(self, "C"+e)/sum)
                
            if eles is None:
                # check if log2x conversion has been made
                if self.chems[0] not in locals():
                    self.log2x()
                    eles = self.chems
            else :
                # rerun conversion to x
                # as ele list may now be new
                self.log2x(eles)
                        
            x2c(self,MO,eles)
        
        else :
            if eles is None:
                eles = self.chems

            return [getattr(self, "C"+i) for i in eles]

    def xs(self, eles=None):
        if eles is None:
            eles = self.chems

        return [getattr(self, i) for i in eles]
    

# ---------- Peridotite liquidus ----------
def Tpdliq(P):
    # Schaefer 2016 ApJ construction
    # built on Hirschmann 2000 solidus (high P, >3.8785673242673817)
    # hiP: a = 26.53 K Gpa−1, b = 1825 K
    # loP: a = 104.42 K Gpa−1, b = 1420 K
    # liquidus assumed = solidus + 600
    # Check out Righter, K. Prediction of metal–silicate partition coefficients for siderophile elements: an update and assessment of PT conditions for metal–silicate equilibrium during accretion of the Earth. Earth
    a_hp = 26.53
    b_hp = 1825.
    a_lp = 104.42
    b_lp = 1420.

    pc = (b_hp-b_lp)/a_lp

    if P > pc:
        T = b_hp + a_hp*(P-pc)
    else:
        T = b_lp + a_lp*P

    return T #- 600  #Fudge to see what happens!
