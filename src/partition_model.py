#!/usr/bin/env python
# -*- coding: utf-8 -*-

# PArtitioning Model for ELemental Abundances = PAMELA

#To switch to alternative Si/O parametrisations:

#Badro Si/O:

# Change partition params to
#Si 4 0.364 -16520 0 0 0 0.28 716 0 0 ba
#O -2 2.736 -11439 0 0 0 0.14 387 0 0 ba
# Change gamma calculation to include np.exp(4.29 - (16500/T))
# Change interaction params to
# g0 0 0 -6.65 0 0 -5.52 0 0 0.36 0 0 0 0 -3.3 -1.61 0 0 -3.22 0 0 2.14943391349987 0 0

# Or for Siebert:
#Si 4 0.549 -12324 0 0 0 0.303 1689 0 0 si
#O -2 0.986 -3275 0 0 0 0.272 502 0 0 si
#Set g0_Si = 0, not -6.65 as above

#General:
# in calculate_abundances in geology_info, comment / uncomment the relevant normalisation section
# Remove Si/O from reading fischer es if not using fischer Si/O
# Change ci.Element.Placeholder to ci.Element.O in this file (and calculate_abundances in geology_info)
# Change geology model config name
# Change o_cap in geology_info

import csv
import numpy as np

import chemistry_info as ci
import geology_info as gi # gi imports this file - ideally remove circularity

class PartitionModel:
    
    def __init__(self, file_params, file_interactions, file_composition, file_es=None):

        # non-bridging oxygens over tetrahedral cations
        self.nbot = 2.7
        
        # fixed activity coefficient of Fe in silicate (Rudge 2010, just above equation G.3)
        # (Actually FeO rather than Fe!)
        self.log10_gammaFe_sil = np.log10(3.0)
        
        self.partitioners = [
            ci.Element.Hf,
            ci.Element.U,
            ci.Element.Ta,
            ci.Element.Pb,
            ci.Element.Nb,
            ci.Element.Si,
            ci.Element.Mn,
            ci.Element.Zn,
            ci.Element.Ga,
            ci.Element.V,
            ci.Element.Cr,
            ci.Element.Cu,
            #ci.Element.Ti,   # Removing Ti - it's extremely erratic and ought to be essentially non-partitioning anyway. NB This should be revisited if Sulfur included, as the Ti/S interaction is strong
            ci.Element.Fe,
            ci.Element.W,
            ci.Element.P,
            ci.Element.Co,
            ci.Element.Ni,
            ci.Element.O,
            ci.Element.C,
            ci.Element.S
        ]
        
        self.non_partitioners = [ # We will assume these are all purely lithophilic
            ci.Element.N,
            ci.Element.Na,
            ci.Element.Mg,
            ci.Element.Al,
            ci.Element.Ti,
            ci.Element.Ca
        ]
        
        self.ele_set = self.partitioners + self.non_partitioners
        
        self.load_partition_parameters(file_params)
        
        self.load_activity_coefficients(file_interactions)

        self.load_default_alloy_composition(file_composition)
        
        self.load_fischer_es(file_es)

    def load_partition_parameters(self, file_params):
        with open(file_params, encoding='utf-8') as csvfile:
            read = csv.reader(csvfile, delimiter=' ')
            i = 0
            for row in read:
                if i == 0:
                    self.params = list()
                    for j in range(1, len(row)):
                        self.params.append(row[j].lower())
                        setattr(self, row[j].lower(), {})
                if i>0:
                    for j in range(1, len(row) - 1):
                        #This logic is a fancy way of setting the class member variable names to match the column headings, but a bit pointless given that in mkd we assume they're called self.a etc anyway...
                        getattr(self, self.params[j-1])[ci.Element[row[0]]] = float(row[j])
                    getattr(self, self.params[len(row)-2])[ci.Element[row[0]]] = str(row[len(row) - 1])
                i += 1

    def load_activity_coefficients(self, file_interactions):
        # get activity coefficients for elements as a function of T (assuming comp)
        with open(file_interactions, encoding='utf-8') as csvfile:
            read = csv.reader(csvfile, delimiter=' ')
        
            i = 0
            for row in read:
                if i == 0:
                    # record reference temperature
                    self.T0 = float(row[0])
                    eles = []
                    for j in range(1, len(row)):
                        eles.append(ci.Element[row[j]])
                    # It might be more readable to make self.eps a dict with (Element, Element) keys, but this is slow
                    self.eps = np.zeros((len(eles), len(eles)))
                if i > 0 and row[0] != "g0":
                    self.eps[i-1] = row[1:]
                if row[0] == "g0":
                    self.logg0 = np.array([])
                    self.logg0name = []
                    for j in range(1, len(row)):
                        self.logg0name.append(eles[j-1])
                        self.logg0 = np.append(self.logg0, float(row[j]))
                i += 1

    def load_fischer_es(self, file_es=None):
        self.fischer_es = dict()
        if file_es is None:
            return
        forbidden_e_els = [ci.Element.Si, ci.Element.O]  # Ignore es for these elements! We will stick with the default interactions
        with open(file_es, encoding='utf-8') as csvfile:
            read = csv.reader(csvfile, delimiter=' ')
            i = 0
            for row in read:
                if i == 0:
                    # record reference temperature
                    self.e_T0 = float(row[0])
                    eles = []
                    for j in range(1, len(row)):
                        eles.append(ci.Element[row[j]])
                if i > 0:
                    for j in range(1, len(row)):
                        el_i = ci.Element[row[0]]
                        el_j = eles[j-1]
                        if (el_i not in forbidden_e_els) and (el_j not in forbidden_e_els):
                            self.fischer_es[(el_i, el_j)] = float(row[j])
                i += 1
            
    def load_default_alloy_composition(self, file_composition):
        with open(file_composition, encoding='utf-8') as csvfile:
            read = csv.reader(csvfile, delimiter=' ')
    
            self.alloy = np.array([])
            self.alloy_els = list()
            for row in read:
                self.alloy_els.append(ci.Element[row[0]])
                self.alloy = np.append(self.alloy, np.array(float(row[1])))
    
    def calculate_gammas(self, T=None, composition=None):
        if T is None:
            T = self.T0
        if composition is None:
            x = self.alloy[1:]
        else:
            x = composition[1:]
        xk = np.tile(x,[len(x),1])
        xj = np.transpose(xk)
        eps = self.construct_eps_table(T, x)
        with np.errstate(divide='ignore',invalid='ignore'):
            log_1_minus_xj_over_xj = np.log(1-xj)/xj
            log_1_minus_xk_over_xk = np.log(1-xk)/xk
            one_over_1_minus_xj = 1/(1-xj)
            one_over_1_minus_xk = 1/(1-xk)
            eps_xj2_xk2 = eps * xj**2 * xk**2
            one = np.sum( np.diag(eps)*(x + np.log(1-x)) )
            
            two = - eps*xk*xj * (1 + log_1_minus_xj_over_xj + log_1_minus_xk_over_xk) + \
                  0.5 * eps_xj2_xk2 * (one_over_1_minus_xj + one_over_1_minus_xk - 1)

            two[np.isnan(two)] = 0.
            two = np.triu(two,1)
            two = np.sum(two)
            
            three = eps*xj*xk*(1+ log_1_minus_xk_over_xk -1/(1-xj)) - \
                    eps_xj2_xk2 * (one_over_1_minus_xj + one_over_1_minus_xk + xj/(2.*(1-xj)**2) - 1)
            three[np.isnan(three)] = 0.
            three = three-np.diag(np.diag(three))
            three = np.sum(three)

            loggFe = one + two + three
            one = - eps*xk*(1 + log_1_minus_xk_over_xk - 1/(1-xj)) +\
                  eps*xk**2*xj*(one_over_1_minus_xj + one_over_1_minus_xk + xj/(2*(1-xj)**2)-1)
            one[np.isnan(one)] = 0.
            one = one - np.diag(np.diag(one))
            one = np.sum(one,axis=1)
            logg = loggFe + ((self.T0/T)*self.logg0) - np.transpose(np.diag(eps))*np.log(1-x) + np.transpose(one)
            gammas = dict(zip(self.logg0name, np.exp(logg)))
            gammas[ci.Element.Fe] = np.exp(loggFe)
            #gammas[ci.Element.O] *= np.exp(4.29 - (16500/T))  # This is a hack to account for Badro's gamma parametrisation for O
            return gammas

    def construct_eps_table(self, T, alloy_composition_without_Fe):
        scaling_factor = self.T0/T  # T0 is the reference temperature
        eps = scaling_factor*self.eps
        for el_i, el_j in self.fischer_es.keys():
            i = self.logg0name.index(el_i)
            j = self.logg0name.index(el_j)
            eps[i][j] = self.calculate_eps_from_e(el_i, el_j, T)
        
        # Damp Si/O interaction above certain core concentrations:
        # This causes convergence problems! Disabling it for now
        damping = False
        if damping:
            si_index = self.alloy_els.index(ci.Element.Si) - 1  # minus 1 because alloy_composition_without_Fe is the alloy without the first element (i.e. Fe)
            o_index = self.alloy_els.index(ci.Element.O) - 1
            X_Si_core = alloy_composition_without_Fe[si_index]
            X_O_core = alloy_composition_without_Fe[o_index]
            damping_factor = self.interaction_damping_function(X_Si_core, X_O_core)
            si_index_2 = self.logg0name.index(ci.Element.Si)
            o_index_2 = self.logg0name.index(ci.Element.O)
            eps[si_index_2][o_index_2] *= damping_factor
            eps[o_index_2][si_index_2] *= damping_factor
        return eps
    
    def interaction_damping_function(self, X_Si_core, X_O_core):
        # How much to damp interactions between Si and O
        # This is arbitrary
        # Should return a value from 0 to 1
        total_SiO = X_Si_core + X_O_core
        return max(0, 1 - (4*total_SiO))
    
    def calculate_eps_from_e(self, el_i, el_j, T):
        # Fischer+ 2015 equation 4
        e = self.fischer_es[(el_i, el_j)]
        M_j = ci.get_element_mass(el_j) # Needs to be the upper index
        return ((e*M_j*self.e_T0)/(0.242*T)) - (M_j/55.85) + 1
    
    def calculate_D(self, P, T, dIW, nbot, log10_gammaFe_sil, core_composition=None, mantle_composition=None, params=None):
        first_iteration = core_composition is None # This should be a functional proxy
        
        if core_composition is None:
            alloy_composition = self.alloy
            core_composition = dict()
            for el_index, el in enumerate(self.alloy_els):
                core_composition[el] = alloy_composition[el_index]
        else:
            alloy_composition = self.convert_core_composition_to_alloy_composition(core_composition)
        
        gammas = self.calculate_gammas(T, alloy_composition)
        els_to_iterate_over = [el for el in self.partitioners if el in self.logg0name + [ci.Element.Fe]]  # Filter out elements that we don't have partitioning info for
        d_dict = dict()
        for e in els_to_iterate_over:
            p_v = params[4] if params is not None else self.v
            if e != ci.Element.Fe:
                logkd_app, sigma_logkd_app_2 = self.calculate_logKD_app(P, T, nbot, e, params) # When iterating over composition, we don't need to recalculate this
                loggamma = np.log10(gammas[e])
                loggammaFe = np.log10(gammas[ci.Element.Fe])
                if e == ci.Element.Placeholder:  # This used to be Oxygen - for Fischer O parametrisation we need to use this logic!
                    logX_met = logkd_app - loggamma - log10_gammaFe_sil + (0.5*dIW) + (2*loggammaFe)
                    #sigma_log_dfe_2 = ((self.sigv[e]/4.)**2)*dIW**2
                    #sigma_logX_met = np.sqrt(sigma_logkd_app_2 + sigma_log_dfe_2)
                    # Need to cap this somewhere below 100% of the core:
                    logX_met = min(-1, logX_met)
                    d_dict[e] = 10**logX_met
                    #mkdSd_o[e] = (10**(logX_met + sigma_logX_met), 10**(logX_met - sigma_logX_met))
                elif e == ci.Element.S:
                    # See Boujibar+ 2014 eqs 6 and 11
                    # gammas not needed - but they are used for how other elements interact with S
                    
                    if mantle_composition is None:
                        # Use a very rough approximation of Earth. Doesn't matter much - this is just an initial guess
                        X_FeO = 0.06
                        X_FeO_sil_wt_percent = 2
                        X_CaO = 0.03
                        X_MgO = 0.2
                        X_TiO2 = 0.001
                        X_Na2O = 0.002
                        X_K2O = 0.001
                    else:
                        # The logic here is that I'm assuming all these elements in the silicate are fully oxidised.
                        # (And that this accounts for all the Oxygen)
                        # So to get X for FeO, it's just the fraction of Fe out of everything that's not oxygen
                        # So X_FeO = X_Fe/(1 - X_O)
                        # This can probably afford to be v. rough: it's only to do part of the calculation for Sulfur after all
                        X_O = mantle_composition[ci.Element.O]
                        oxygen_modifier = 1/(1-X_O)
                        
                        # If these are not present in the dict at this stage, it means they are not present
                        # in the bulk composition. Hence a default value of 0
                        
                        X_FeO = oxygen_modifier*mantle_composition.get(ci.Element.Fe, 0)
                        X_CaO = oxygen_modifier*mantle_composition.get(ci.Element.Ca, 0)
                        X_MgO = oxygen_modifier*mantle_composition.get(ci.Element.Mg, 0)
                        X_TiO2 = oxygen_modifier*mantle_composition.get(ci.Element.Ti, 0)
                        # Dividing the next ones by 2 because there's half as many Na2(O) as Na(O)
                        X_Na2O = 0.5*oxygen_modifier*mantle_composition.get(ci.Element.Na, 0)
                        X_K2O = 0.5*oxygen_modifier*mantle_composition.get(ci.Element.K, 0)
                        # I'm assuming the FeO wt % is directly proportional to the molar %
                        # and taking the constant from cells K3 and G3 here: https://www.researchgate.net/publication/308404683_Sulfur_partitioning_calculator
                        # This is potentially inaccurate
                        X_FeO_sil_wt_percent = 127.5862069*X_FeO
                    
                    core_composition_by_mass = self.convert_composition_to_mass(core_composition)
                    X_Si = core_composition_by_mass[ci.Element.Si]
                    X_C = core_composition_by_mass[ci.Element.C]
                    X_Fe = core_composition_by_mass[ci.Element.Fe]
                    X_Ni = core_composition_by_mass[ci.Element.Ni]
                    X_O = core_composition_by_mass[ci.Element.O]
                    
                    
                    
                    logCs = self.calculate_logCs(X_FeO, X_CaO, X_MgO, X_TiO2, X_Na2O, X_K2O)
                    S_interactions = self.calculate_S_interactions(X_Si, X_C, X_Fe, X_Ni, X_O)
                    
                    logD = np.log10(X_FeO_sil_wt_percent) - logCs + logkd_app + S_interactions
                    d_dict[e] = 10**logD
                else:
                    additional_term = ((1 - (0.5*p_v[e]))*loggammaFe) if self.source[e] in ['f1', 'f2', 'b'] else 0  # We take care of gamma_0 elsewhere by setting it to 0 in the input file
                    logD = logkd_app - loggamma + (0.5*p_v[e]*(log10_gammaFe_sil - (0.5*dIW))) + additional_term
                    #sigma_log_dfe_2 = ((self.sigv[e]/4.)**2)*dIW**2
                    #sigma_logD = np.sqrt(sigma_logkd_app_2 + sigma_log_dfe_2)
                    d_dict[e] = 10**logD
                    #mkdSd_o[e] = (10**(logD + sigma_logD), 10**(logD - sigma_logD))
            else:
                loggammaFe = np.log10(gammas[ci.Element.Fe])
                logdfe = self.calculate_logDFe(log10_gammaFe_sil, dIW, loggammaFe) # When iterating over composition, we don't need to recalculate this
                #sigma_log_dfe_2 = ((self.sigv[e]/4.)**2)*dIW**2 #This isn't very nice
                d_dict[e] = 10**logdfe
                #mkdSd_o[e] = (10**(logdfe + sigma_log_dfe_2), 10**(logdfe - sigma_log_dfe_2))
        return d_dict
        
    # NB: Corgne+ 2008 parametrise KD_app as shown here, but Fischer+ 2015 actually parametrise KD instead, which is why the fischer_update corrections exist
    def calculate_logKD_app(self, P, T, nbot, e, params=None):
        # allow parameters to be updated to propagate error through model
        if params is None:  # Could remove this
            p_a = self.a
            p_b = self.b
            p_c = self.c
            p_d = self.d
        else:
            p_a = params[0]
            p_b = params[1]
            p_c = params[2]
            p_d = params[3]
        logkd_app = p_a[e] + p_b[e]/T + p_c[e]*P/T + p_d[e]*nbot
        sigma_logkd_app_2 = self.siga[e]**2 + nbot**2*self.sigd[e]**2 + self.sigb[e]**2*(1/T**2) + self.sigc[e]**2*(P/T)**2  # Error squared
        return logkd_app, sigma_logkd_app_2
        
    def calculate_logDFe(self, log10_gammaFe_sil, dIW, loggammaFe):
        # From Equation G.1 in Rudge 2010:
        # dIW := 2log10(gammaFe_sil/gammaFe) + 2log10(1/D_Fe) --> 0.5dIW = log10(gammaFe_sil) - loggammaFe - logD_Fe --> logD_Fe = log10(gammaFe_sil) - 0.5*dIW - loggammaFe
        logdfe = log10_gammaFe_sil - ((0.5*dIW) + loggammaFe)
        return logdfe
        
    def calculate_logCs(self, X_FeO, X_CaO, X_MgO, X_TiO2, X_Na2O, X_K2O): # Args are all number fractions in the silicate
        #Boujibar+ 2014 eq 6
        return (3.15*X_FeO) + (2.65*X_CaO) + (0.12*X_MgO) + (0.77*X_TiO2) + (0.75*(X_Na2O + X_K2O)) - 5.704
        
    def calculate_S_interactions(self, X_Si, X_C, X_Fe, X_Ni, X_O): # Args are all mass fractions in metal
        d = 32
        e = 181
        f = 305
        g = 30.2
        h = 1.13
        i = 10.7
        j = 31.4
        k = -3.72
        si_terms = d*np.log10(1 - X_Si) + e*((np.log10(1 - X_Si))**2) + f*((np.log10(1 - X_Si))**3)
        toret = si_terms + g*np.log10(1 - X_C) + h*np.log10(1 - X_Fe) + i*np.log10(1 - X_Ni) + j*np.log10(1 - X_O) + k
        return toret
        
    def peridotite_liquidus(self, P):
        # Schaefer 2016 ApJ construction
        # built on fitting Hirschmann 2000 data (high P, >3.8785673242673817 GPa, corresponds to 100km depth)
        # high P: a = 26.53 K Gpa−1, b = 1825 K
        # low P: a = 104.42 K Gpa−1, b = 1420 K
        # liquidus assumed = solidus + 600
        # Check out Righter, K. Prediction of metal–silicate partition coefficients for siderophile elements: an update and assessment of PT conditions for metal–silicate equilibrium during accretion of the Earth.
        a_hp = 26.53
        b_hp = 1825.0
        a_lp = 104.42
        b_lp = 1420.0
        liquidus_minus_solidus = 600
        
        pc = (b_hp - b_lp)/a_lp
        if P > pc:
            T = b_hp + (a_hp*(P - pc))
        else:
            T = b_lp + (a_lp*P)
        T += liquidus_minus_solidus
        return T
        
        #return 1621 + (38.415*P) + (0.00038369*P*P*P) - (0.1958*P*P)   # eqn. 7 in Badro 2015 appendix
        
        #return 2022 + (54.21*P) + (0.00090747*P*P*P) - (0.34*P*P)  # eqn 8
        
        #return 1940*((1+(P/29))**(1/1.9))  # eqn 9
        
        #a = 2022 + (54.21*P) + (0.00090747*P*P*P) - (0.34*P*P)
        #b = 1940*((1+(P/29))**(1/1.9))
        #return 0.5*(a + b)  # eqn 10
        
        #return 1973 + (28.57*P)  # Eqn H.1 in Rudge 2010
    
    def calculate_partition_coefficients(self, pressure, fO2, core_composition=None, mantle_composition=None, temp=None, nbot=None):
        # for now assume that T is fixed at the pressure of ...
        T = self.peridotite_liquidus(pressure) if temp is None else temp
        nbot_to_use = self.nbot if nbot is None else nbot
        Ds = self.calculate_D(pressure, T, fO2, nbot_to_use, self.log10_gammaFe_sil, core_composition, mantle_composition)
        return Ds

    def get_all_partition_coefficients(self, pressure, fO2, abundances=None, temp=None, nbot=None):
        core_composition = self.convert_abundance_dict_to_layer_composition(gi.Layer.core, abundances)
        mantle_composition = self.convert_abundance_dict_to_layer_composition(gi.Layer.mantle, abundances)
        Ds = self.calculate_partition_coefficients(pressure, fO2, core_composition, mantle_composition, temp, nbot)
        for element in self.non_partitioners:
            Ds[element] = 0
        return Ds
    
    def convert_abundance_dict_to_layer_composition(self, layer, abundances=None):
        if abundances is None or len(abundances) == 0:
            return None
        composition = dict()
        for element in self.ele_set:
            try:
                composition[element] = abundances.get(element, {layer: 0})[layer]
            except KeyError:
                # This means that the layer information is missing for some/all elements. i.e. there is None
                return None
        return composition
    
    def convert_core_composition_to_alloy_composition(self, core_composition):
        alloy = np.array([core_composition.get(element, 0) for element in self.alloy_els])
        return alloy
        
    def convert_composition_to_mass(self, composition):
        # Converts a composition by number to a composition by mass
        # TODO add a test for this
        total_mass = 0
        for el, number_frac in composition.items():
            total_mass += number_frac*ci.get_element_mass(el)
        toret = dict()
        for el, number_frac in composition.items():
            toret[el] = (number_frac*ci.get_element_mass(el))/total_mass
        return toret
