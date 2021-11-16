import numpy as np
import csv

class D():

    # readin list of parameters (
    def __init__(self, file_params, file_interactions, file_composition, fischer_update=False):

        # load partition coefficient params
        with open(file_params) as csvfile:
            read = csv.reader(csvfile, delimiter=' ')
            
            i = 0
            for row in read:
                if i == 0:
                    self.params = []
                    for j in range(1,len(row)):
                        self.params.append(row[j].lower())
                        setattr(self, row[j].lower(),{})
                if i>0:
                    for j in range(1,len(row)):
                        getattr(self, self.params[j-1])[row[0].lower()] = float(row[j])
                i += 1
        
        # get activity coefficients for elements as a function of T (assuming comp)
        with open(file_interactions) as csvfile:
            read = csv.reader(csvfile, delimiter=' ')

            i = 0
            for row in read:
                if i == 0:
                    # record reference temperature
                    self.T0 = float(row[0])

                    eles = []
                    for j in range(1,len(row)):
                        eles.append(row[j].lower())
                    setattr(self, "eps", np.zeros((len(eles),len(eles))))

                if i>0 and row[0] != "g0":
                    self.eps[i-1] = row[1:]

                if row[0] == "g0":
                    self.logg0 = np.array([])
                    self.logg0name = []
                    for j in range(1,len(row)):
                        (self.logg0name).append(eles[j-1])
                        self.logg0 = np.append(self.logg0, float(row[j]))
                        
                i += 1

        # load composition of alloy
        with open(file_composition) as csvfile:
            read = csv.reader(csvfile, delimiter=' ')

            self.C = np.array([])
            self.Cname = []
            for row in read:
                self.Cname.append(row[0])
                self.C = np.append(self.C, np.array(float(row[1])))
        # generate log gamma's
        # !! Composition of alloy has to be same length, and ordered the same as the
        #    the epsilon interaction parameters
        # 
        # !! This part of the calculation is not self consistent:
        #    The gammas are generated assuming a metal composition
        #    then the metal composition is calculated from the D's
        #    -> a calculated metal composition that is different from the one assumed in
        #       calculating the D's
        #    This should be an iterative process, whereby an initial guess composition
        #    is used, the D's calculated, an equilibrium metal composition calculated
        #    and then a new set of gammas derived, a new set of D's, a new metal composition
        #    etc.
        #    until the metal composition stabilises !!
        self.fischer_update = fischer_update
        self.g = self.calculate_gammas()

    # If composition is None, assume that we should set gammas based on inital composition, self.C
    def calculate_gammas(self, composition=None):
        if composition is None:
            x =  self.C[1:]
        else:
            x = composition[1:]
        xk = np.tile(x,[len(x),1])
        xj = np.transpose(xk)

        with np.errstate(divide='ignore',invalid='ignore'):
            one = np.sum( np.diag(self.eps)*(x + np.log(1-x)) )
            
            two = - self.eps*xk*xj * (1+np.log(1-xj)/xj+np.log(1-xk)/xk) + \
                  0.5 * self.eps * xj**2 * xk**2 * (1/(1-xj) + 1/(1-xk) - 1)
            two[np.isnan(two)] = 0.
            two = np.triu(two,1)
            two = np.sum(two)

            three = self.eps*xj*xk*(1+np.log(1-xk)/xk-1/(1-xj)) - \
                    self.eps* xj**2 * xk**2 * (1/(1-xj) + 1/(1-xk) + xj/(2.*(1-xj)**2) - 1)
            three[np.isnan(three)] = 0.
            three = three-np.diag(np.diag(three))
            three = np.sum(three)

            loggFe = one + two + three
            
            one = - self.eps*xk*(1 + np.log(1-xk)/xk - 1/(1-xj)) +\
                  self.eps*xk**2*xj*(1/(1-xj) + 1/(1-xk) + xj/(2*(1-xj)**2)-1)
            one[np.isnan(one)] = 0.
            one = one - np.diag(np.diag(one))
            one = np.sum(one,axis=1)

            logg = loggFe + self.logg0 - np.transpose(np.diag(self.eps))*np.log(1-x) + np.transpose(one)
            gammas = dict(zip(self.logg0name,np.exp(logg)))
            # If we're using either fischer update, certain elements are
            # parametrised assuming loggFe and logg0 are absorbed into a, b, c.
            # logg0 is set to 0 in the input file. loggFe is effectively set to 0 here
            if self.fischer_update:
                for fu_element in ['ni', 'co', 'v', 'cr', 'si', 'o']:
                    gammas[fu_element] /= np.exp(loggFe)
            gammas['fe'] = np.exp(loggFe)
            return gammas

    def mkd(self, P, T, dIW, nbot, gammaFe_sil, eles, composition=None, params=None):
        # allow parameters to be updated to propagate error through model
        if params is None:
            p_a = self.a
            p_b = self.b
            p_c = self.c
            p_d = self.d
            p_v = self.v                                    
        else:
            p_a = params[0]
            p_b = params[1]
            p_c = params[2]
            p_d = params[3]
            p_v = params[4]
        if composition is None:
            # Assume we are just going to use gammas corresponding to initial composition
            gammas = self.g
        else:
            gammas = self.calculate_gammas(composition)
        print(gammas)
        # molar partition coefficient
        mkd_o = []
        mkdSd_o = []
        for e in eles:
            if e != "fe" and e in self.logg0name:
                lg10d = p_a[e] + p_b[e]/T + p_c[e]*P/T + p_d[e]*nbot - 0.25*p_v[e]*dIW - ((self.T0/T)*np.log10(gammas[e])) + 0.5*p_v[e]*np.log10(gammaFe_sil)
                mkd_o.append(10**lg10d)

                lg10dSd = np.sqrt(self.siga[e]**2 + nbot**2*self.sigd[e]**2 + self.sigb[e]**2*(1/T**2) + self.sigc[e]**2*(P/T)**2 + ((self.sigv[e]/4.)**2)*dIW**2)
                mkdSd_o.append([10**(lg10d+lg10dSd), 10**(lg10d-lg10dSd)])
            elif e == "fe":
                # this value assumed fixed
                lg10d = -0.5*dIW - ((self.T0/T)*np.log10(gammas['fe'])) + np.log10(gammaFe_sil)
                mkd_o.append(10**lg10d)

                lg10dSd = np.sqrt(self.siga[e]**2 + nbot**2*self.sigd[e]**2 + self.sigb[e]**2*(1/T**2) + self.sigc[e]**2*(P/T)**2 + ((self.sigv[e]/4.)**2)*dIW**2)
                mkdSd_o.append([10**(lg10d+lg10dSd), 10**(lg10d-lg10dSd)])
            else :
                mkd_o.append(0.)
                mkdSd_o.append([0., 0.])

        # errors on molar partition coefficients
        

        return dict(zip(eles,mkd_o)), dict(zip(eles,mkdSd_o))
