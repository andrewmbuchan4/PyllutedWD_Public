import numpy as np

def standalone_gamma_calc():
    composition_80 = np.array([7.73709352e-01, 8.97435091e-03, 1.14860312e-01, 4.58674373e-02, 9.01592320e-03, 4.01384975e-02])
    composition_81 = np.array([6.38889145e-01, 7.90524292e-03, 3.94100562e-02, 2.61550191e-01, 8.63055458e-03, 3.50771442e-02])
    composition = composition_81
    x = composition[1:]
    xk = np.tile(x,[len(x),1])
    xj = np.transpose(xk)
    print('sgc5')
    print(composition)
    print(x)
    print(xk)
    print(xj)
    eps = np.array(
        [[12.80560104, -20.37210851, 9.750942788, -4.856474886, 2.329824276],
        [-20.04102659, -10.48822455, -7.137370937, -11.70920226, 1.399379354],
        [9.689128839, -7.127704002, 12.41143057, -0.01126355494, 1.157649745],
        [-4.854408094, -9.829298952, 0.01891228893, 0.004678592533, -0.002652377115],
        [2.331376474, 1.372436207, 1.191112544, 0.004678592533, 0.1182124273]]
    )
    logg0 = np.array([0, 0, 0, 0, 0])
    logg0name = np.array(['C', 'O', 'Si', 'Cr', 'Ni'])
    fischer_update = True
    fischer_elements = ['Ni', 'Co', 'V', 'Cr', 'Si', 'O', 'C']
    print(eps)
    with np.errstate(divide='ignore',invalid='ignore'):
        print(np.diag(eps))
        print(np.log(1-x))
        print(np.diag(eps)*(x + np.log(1-x)))
        one = np.sum( np.diag(eps)*(x + np.log(1-x)) )
        print(one)
        
        two = - eps*xk*xj * (1+np.log(1-xj)/xj+np.log(1-xk)/xk) + \
              0.5 * eps * xj**2 * xk**2 * (1/(1-xj) + 1/(1-xk) - 1)
        print(-eps*xk*xj)
        print((1+np.log(1-xj)/xj+np.log(1-xk)/xk))
        print(0.5 * eps * xj**2 * xk**2 * (1/(1-xj) + 1/(1-xk) - 1))
        print(two)

        two[np.isnan(two)] = 0.
        two = np.triu(two,1)
        two = np.sum(two)
        
        print(two)
        three = eps*xj*xk*(1+np.log(1-xk)/xk-1/(1-xj)) - \
                eps* xj**2 * xk**2 * (1/(1-xj) + 1/(1-xk) + xj/(2.*(1-xj)**2) - 1)
        print(eps*xj*xk*(1+np.log(1-xk)/xk-1/(1-xj)))
        print(eps* xj**2 * xk**2 * (1/(1-xj) + 1/(1-xk) + xj/(2.*(1-xj)**2) - 1))
        print(three)
        three[np.isnan(three)] = 0.
        three = three-np.diag(np.diag(three))
        three = np.sum(three)
        print(three)
        
        loggFe = one + two + three
        print(loggFe)
        one = - eps*xk*(1 + np.log(1-xk)/xk - 1/(1-xj)) +\
              eps*xk**2*xj*(1/(1-xj) + 1/(1-xk) + xj/(2*(1-xj)**2)-1)
        print(eps*xk*(1 + np.log(1-xk)/xk - 1/(1-xj)))
        print(eps*xk**2*xj*(1/(1-xj) + 1/(1-xk) + xj/(2*(1-xj)**2)-1))
        print(one)
        one[np.isnan(one)] = 0.
        one = one - np.diag(np.diag(one))
        one = np.sum(one,axis=1)
        print(one)
        print('logg cont')
        print(loggFe)
        print(logg0)
        print(- np.transpose(np.diag(eps))*np.log(1-x))
        print(np.transpose(one))
        logg = loggFe + logg0 - np.transpose(np.diag(eps))*np.log(1-x) + np.transpose(one)
        print(logg0name)
        print(np.exp(logg))
        gammas = dict(zip(logg0name, np.exp(logg)))
        # If we're using either fischer update, certain elements are
        # parametrised assuming loggFe and logg0 are absorbed into a, b, c.
        # logg0 is set to 0 in the input file. loggFe is effectively set to 0 here
        # TODO: Double check this correction.
        if fischer_update:
            for fu_element in fischer_elements:
                try:
                    gammas[fu_element] /= np.exp(loggFe)
                except KeyError:
                    pass
        gammas['Fe'] = np.exp(loggFe)
        print('Final gammas')
        print(gammas)
        return gammas

def main():
    standalone_gamma_calc()

if __name__ == '__main__':
    main()
