def E_Al(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.0153-(0.0641)*(N_o)-(0)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.0641
            x_c = 0
        else:
            x_m = 0
            x_o = 0.0153/N_o
            x_c = 0
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.0641)+0.829*((0.0153-(0.0641)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Ti(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.00044-(0.0042)*(N_o)-(0)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.0042
            x_c = 0
        else:
            x_m = 0
            x_o = 0.00044/N_o
            x_c = 0
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.0042)+0.829*((0.00044-(0.0042)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Ca(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.011099-(0.04452)*(N_o)-(0)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.04452
            x_c = 0
        else:
            x_m = 0
            x_o = 0.011099/N_o
            x_c = 0
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.04452)+0.829*((0.011099-(0.04452)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Na(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >=N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.002037-(0.01771)*(N_o)-(0)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.01771
            x_c = 0
        else:
            x_m = 0
            x_o = 0.002037/N_o
            x_c = 0
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.01771)+0.829*((0.002037-(0.01771)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_O(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.482879-(0.6011)*(N_o)-(0)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.6011
            x_c = 0
        else:
            x_m = 0
            x_o = 0.482879/N_o
            x_c = 0
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.6011)+0.829*((0.482879-(0.6011)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final


def E_Mg(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.16482-(0.04167)*(N_o)-(0)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.04167
            x_c = 0
        else:
            x_m = 0
            x_o = 0.16482/N_o
            x_c = 0
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.04167)+0.829*((0.16482-(0.04167)*(0.001)-(0)*(0.17))/(0.829))+0.17*(0))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0.0001
    return E_final

def E_Fe(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.149057-(0.0314)*(N_o)-(0.7676)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.0314
            x_c = 0.7676
        else:
            x_m = 0
            x_o = (0.149057-0.7676*(N_c))/N_o
            x_c = 0.7676
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.0314)+0.829*((0.149057-(0.0314)*(0.001)-(0.7676)*(0.17))/(0.829))+0.17*(0.7676))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Ni(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.008066-(0.0000371)*(N_o)-(0.0444)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.0000371
            x_c = 0.0444
        else:
            x_m = 0
            x_o = 0.0000371
            x_c = (0.008066-(0.0000371)*N_o)/(N_c)
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.0000371)+0.829*((0.008066-(0.0000371)*(0.001)-(0.0444)*(0.17))/(0.829))+0.17*(0.0444))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Cr(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.002351-(0.000139)*(N_o)-(0.00868)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.000139
            x_c = 0.00868
        else:
            x_m = 0
            x_o = (0.002351-0.00868*(N_c))/(N_o)
            x_c = 0.00868
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.000139)+0.829*((0.002351-(0.000139)*(0.001)-(0.00868)*(0.17))/(0.829))+0.17*(0.00868))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Si(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.149118-(0.181509)*(N_o)-(0.1071)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0.181509
            x_c = 0.1071
        else:
            x_m = 0
            x_o = (0.149118-0.1071*(N_c))/(N_o)
            x_c = 0.1071
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0.181509)+0.829*((0.149118-(0.181509)*(0.001)-(0.1071)*(0.17))/(0.829))+0.17*(0.1071))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_C(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.001581-(0)*(N_o)-(0.0083482)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0
            x_c = 0.0083482
        else:
            x_m = 0
            x_o = 0
            x_c = 0.001581/N_c
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0)+0.829*((0.001581-(0)*(0.001)-(0.0083482)*(0.17))/(0.829))+0.17*(0.0083482))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final

def E_Nz(N_c,N_o,f_c,f_o):
    if f_c >= 0 and f_o >= 0 and f_c + f_o <= 1 and 0.19 >= N_c >= 0 and N_o >= 0 and N_c + N_o <= 1:
        f_m = 1 - f_c - f_o
        N_m = 1 - N_c - N_o
        x_m = (0.000046429-(0)*(N_o)-(0.0002684)*(N_c))
        if x_m >= 0:
            x_m = x_m
            x_o = 0
            x_c = 0.0002684
        else:
            x_m = 0
            x_o = 0
            x_c = 0.000046429/N_c
        E_1 = (f_o*(x_o)+f_m*(x_m/(N_m))+f_c*(x_c))
        E_2 = (0.001*(0)+0.829*((0.000046429-(0)*(0.001)-(0.0002684)*(0.17))/(0.829))+0.17*(0.0002684))
        E_final = (E_1)/(E_2)
    else:
        E_final = 0
    return E_final
