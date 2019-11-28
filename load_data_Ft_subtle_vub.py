"""
This code is fitting (1) n2 on in vivo data and (2) n1 on in vitro data
"""
#https://matplotlib.org/3.1.0/gallery/ticks_and_spines/tick-locators.html
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
# https://pythonspot.com/read-excel-with-pandas/
from scipy.optimize import curve_fit
from scipy.optimize import dual_annealing

from math import log10, floor
import matplotlib.ticker as ticker

import sympy as sym
from sympy.solvers import solve

def round_sig(x, sig=2):
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


# for latex-rendering. (I guess with a lot of redundancy.. but it works!)


import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from matplotlib import rc
from matplotlib import rc

if False:

    import latex_rendering
    rc('text', usetex=True)  #
    plt.rcParams[
        'text.usetex'] = True  # this line turns "ticks" on axis into latex  (i think is it the same function as the line above. not sure. )
    plt.rcParams['text.latex.unicode'] = True
    # for the size of the exported figures
    plt.rcParams['figure.figsize'] = (24, 18)
    rc('text.latex', preamble=r'\usepackage{cmbright}')  # this line turns text in legends for instance into latex


    # these definitions will be use to have compact format for the "ticks" on plots.
    def round_to_n(x, n):
        " Round x to n significant figures "
        return round(x, -int(math.floor(np.sign(x) * math.log10(abs(x)))) + n)


    def str_fmt(x, n=2):
        " Format x into nice Latex rounding to n"
        if x < 0.1:
            if x != .01 and x != .001 and x != .0001 and x != .00001:
                power = int(np.log10(round_to_n(x, 0))) - 1
                f_SF = round_to_n(x, n) * pow(10, -power)
                return r"${} \, 10^{{ {} }}$".format(f_SF, power)
            else:
                power = int(np.log10(round_to_n(x, 0)))
                f_SF = round_to_n(x, n) * pow(10, -power)
                return r"${} \, 10^{{ {} }}$".format(f_SF, power)
        elif x > 99:
            power = int(np.log10(round_to_n(x, 0)))
            f_SF = round_to_n(x, n) * pow(10, -power)
            return r"${} \, 10^{{ {} }}$".format(f_SF, power)
        else:
            return r"${}$".format(x)


    def str_fmt_round_numb(x, n=1):
        " Format x into nice Latex rounding to n"
        if x < 0.1:
            power = int(np.log10(round_to_n(x, 0))) - 1
            f_SF = round_to_n(x, n) * pow(10, -power)
            return r"${} \, 10^{{ {} }}$".format(f_SF, power)
        elif x > 99:
            power = int(np.log10(round_to_n(x, 0)))
            f_SF = round_to_n(x, n) * pow(10, -power)
            return r"${} \, 10^{{ {} }}$".format(f_SF, power)
        else:
            return r"${}$".format(x)

cellnames = ['a6.5', 'a6.6', 'a6.7', 'a6.8']
numcells = 4
surf_inputs = np.zeros(numcells)
cell_surfaces_g = {'a6.5':3164.5,'a6.6':1086,'a6.7':1750.5,'a6.8':275}
for kk in range(numcells):
    surf_inputs[kk] = cell_surfaces_g[cellnames[kk]]

hill_params = ['nn', 'KK', 'intercept' ,'xmax']
full_params = ['n1','K2','cst','erkmax','n2','K1']
full_params_t = ['n1','Kt','erkmax','cst']
subtle_params = ['fgf', 'c', 'Kb2', 'erkmax', 'cst']


symbolic_solve = False
if symbolic_solve == True:
    Fs = sym.Symbol('F')
    Fts = sym.Symbol('Ft')
    Rs = sym.Symbol('Rs')
    Rts = sym.Symbol('Rt')
    Kb2s = sym.Symbol('Kb2')
    cs = sym.Symbol('c')

    myequation = (Fts-Fs)*(1+2*Fs/cs/Kb2s + Fs**2/cs/Kb2s**2) - Fs*(2*Rts/cs/Kb2s + 2*Fs*Rts/cs/Kb2s**2)

    sol=solve(myequation,Fs,dict=True)

    print('my sols', sol)


def F1(Ft,Rt,Kb2,c):
    temp1 = (-27*Ft*Kb2**2*c - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt) + 2*(-Ft + 2*Kb2 + 2*Rt)**3)**2 - 4*(6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)**3
    print('temp1', temp1)
    temp2 = (-27*Ft*Kb2**2*c - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt) + 2*(-Ft + 2*Kb2 + 2*Rt)**3)**2 - 4*(6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)**3
    return Ft/3 - 2*Kb2/3 - 2*Rt/3 - (6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)/(3*(-27*Ft*Kb2**2*c/2 + np.sqrt(temp1)/2  - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt)/2 + (-Ft + 2*Kb2 + 2*Rt)**3)**(1/3)) - (-27*Ft*Kb2**2*c/2 + np.sqrt(temp1)/2 - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt)/2 + (-Ft + 2*Kb2 + 2*Rt)**3)**(1/3)/3


def F2(Ft,Rt,Kb2,c):
    return Ft/3 - 2*Kb2/3 - 2*Rt/3 - (6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)/(3*(-1/2 - np.sqrt(3)*I/2)*(-27*Ft*Kb2**2*c/2 + np.sqrt((-27*Ft*Kb2**2*c - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt) + 2*(-Ft + 2*Kb2 + 2*Rt)**3)**2 - 4*(6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)**3)/2 - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt)/2 + (-Ft + 2*Kb2 + 2*Rt)**3)**(1/3)) - (-1/2 - np.sqrt(3)*I/2)*(-27*Ft*Kb2**2*c/2 + np.sqrt((-27*Ft*Kb2**2*c - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt) + 2*(-Ft + 2*Kb2 + 2*Rt)**3)**2 - 4*(6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)**3)/2 - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt)/2 + (-Ft + 2*Kb2 + 2*Rt)**3)**(1/3)/3

def F3(Ft,Rt,Kb2,c):
    return Ft/3 - 2*Kb2/3 - 2*Rt/3 - (6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)/(3*(-1/2 + np.sqrt(3)*I/2)*(-27*Ft*Kb2**2*c/2 + np.sqrt((-27*Ft*Kb2**2*c - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt) + 2*(-Ft + 2*Kb2 + 2*Rt)**3)**2 - 4*(6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)**3)/2 - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt)/2 + (-Ft + 2*Kb2 + 2*Rt)**3)**(1/3)) - (-1/2 + np.sqrt(3)*I/2)*(-27*Ft*Kb2**2*c/2 + np.sqrt((-27*Ft*Kb2**2*c - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt) + 2*(-Ft + 2*Kb2 + 2*Rt)**3)**2 - 4*(6*Ft*Kb2 - 3*Kb2**2*c - 6*Kb2*Rt + (-Ft + 2*Kb2 + 2*Rt)**2)**3)/2 - (-9*Ft + 18*Kb2 + 18*Rt)*(-2*Ft*Kb2 + Kb2**2*c + 2*Kb2*Rt)/2 + (-Ft + 2*Kb2 + 2*Rt)**3)**(1/3)/3

def hill(x, nn, KK, intercept ,xmax):
    return xmax*x**nn/(x**nn + KK**nn) + intercept
def hill2(x, nn, KK):
    return (x)**nn/((x)**nn + (KK)**nn)

def fracact(L,K,c):
    return L**2/(c*K**2+2*K*L + L**2)

def fracact2(Ft,K,c,Rt):
    R = Rt
    bb = 1+2*R/(c*K)
    aa = 2*R/(c*K**2)
    cc = -Ft
    myF = (-bb + np.sqrt(bb**2-4*aa*cc))/(2*aa)
    return myF**2/(c*K**2+2*K*myF + myF**2)

def fracact2_vivo(Rt,Ft,K,c):
    R = Rt# /(1+2f/c/K + f**2/c/K**2)
    bb = 1+2*R/(c*K)
    aa = 2*R/(c*K**2)
    cc = -Ft
    myF = (-bb + np.sqrt(bb**2-4*aa*cc))/(2*aa)
    return myF**2/(c*K**2+2*K*myF + myF**2)

myL = np.logspace(-3,12)
#myY = np.array([fracact(myL[kk],100,.2) for kk in range(len(myL))])
#myY2 = np.array([fracact2(myL[kk],100,.2) for kk in range(len(myL))])
# def erk_act(fgf,n1,K2,cst,erkmax,n2,K1): #(fgf,n1,RT,K1,K2,cst,erkmax)
#     RT =1
#     fr = RT*fgf**n1/(K1**n1 + fgf**n1)
#     return erkmax*fr**n2/(K2**n2 + fr**n2) + cst

def erk_act(fgf,n1,Kt,erkmax,cst): #(fgf,n1,RT,K1,K2,cst,erkmax)
    n2= 2.6
    K1t =  (3.78/Kt)**(1/n2)
    return erkmax*fgf**(n1*n2)/(Kt*(K1t+fgf**n1)**n2 + fgf**(n1*n2)) + cst

def erk_act_subtle(Rt,fgf,c,Kb2,erkmax,cst,K2,n2): #(fgf,n1,RT,K1,K2,cst,erkmax)
    myFt = F1(fgf,Rt,Kb2,c)
    myFR = Rt*myFt**2/(myFt**2+2*myFt*Kb2+c*Kb2**2)
    return erkmax*myFR**n2/(myFR**n2 + K2**n2) + cst


effective_function = True

if effective_function==True:
    xi = (-1.0 + (-3.0 + 0j) ** (1 / 2.0)) / 2.0
    Rt_f = 1000000.0
    Kb2_f = 100000.0
    c_f = 0.01
    fgf_list = np.logspace((-3), 9, 100)
    outputs = np.ones(len(fgf_list))
    outputs1 = np.ones(len(fgf_list))
    outputs3 = np.ones(len(fgf_list))
    free_receptors = np.ones(len(fgf_list))
    bound_receptors = np.ones(len(fgf_list))
    doubly_bound_receptors = np.ones(len(fgf_list))
    for kk in range(len(fgf_list)):
        fgf_f = fgf_list[kk]
        AA = 1 /c_f/Kb2_f**2
        BB = (2/c_f/Kb2_f + 2*Rt_f/c_f/Kb2_f**2 - fgf_f/c_f/Kb2_f**2)
        CC = (1 + 2*Rt_f/c_f/Kb2_f -2*fgf_f/c_f/Kb2_f)
        DD = - fgf_f

        delta0 = BB**2 - 3*AA*CC
        delta1 = 2*BB**3-9*AA*BB*CC + 27*AA**2*DD

        preC = delta1**2 - 4*delta0**3 + 0j
        Croot = ((delta1 + preC**(.5))/2)**(1/3.0)



        solX1 = -1/3/AA*(BB + Croot + delta0/Croot)
        solX2 = -1 / 3 / AA * (BB + xi*Croot + delta0 / Croot/xi)
        solX3 = -1 / 3 / AA * (BB + xi**2*Croot + delta0 / Croot/xi**2)
        outputs[kk]=solX2.real
        outputs1[kk]=solX1.real
        outputs3[kk]=solX3.real
        free_receptors = Rt_f/(1 + 2*outputs/c_f/Kb2_f + outputs**2/c_f/Kb2_f**2)
        bound_receptors = 2*outputs*free_receptors/c_f/Kb2_f
        doubly_bound_receptors = outputs**2*free_receptors/c_f/Kb2_f**2
        sum_receptors = free_receptors + bound_receptors + doubly_bound_receptors
        frac_free_receptors = free_receptors/sum_receptors
        frac_bound_receptors = bound_receptors/sum_receptors
        frac_doubly_bound_receptors = doubly_bound_receptors/sum_receptors




        mybounds_s = [(0,5),(0,10**9)]
        #pann = dual_annealing(error_func_hill_s,bounds =mybounds_s)
        #print(pann.x)pann.x[1]

        model1 = frac_doubly_bound_receptors+1/2*frac_bound_receptors


        def error_func_hill_s(p, Y=model1, X=fgf_list):
            Ypredicted = [hill2(X[kk], p[0], p[1]) for kk in range(len(X))]  #
            diff = Y - Ypredicted
            return np.sum(np.square(diff))

        #pann = dual_annealing(error_func_hill_s,bounds =mybounds_s)
        #print(pann.x)#[4.45043503e-01 1.04100690e+04]
        #predicted = [hill2(fgf_list[kdk],pann.x[0],pann.x[1]) for kdk in range(len(fgf_list))]

        if solX2.imag > 10**(-8):
            print('should be small', solX2.imag)
        if False:
            print('solX', solX1)
            print('sol2 ',solX2)
            print('sol3', solX3)
            print('check')
            print('sol 1 : ', AA*solX1**3 + BB*solX1**2 + CC*solX1+DD)
            print('sol 2 : ', AA * solX2 ** 3 + BB * solX2 ** 2 + CC * solX2 + DD)
            print('sol 3 : ', AA * solX3 ** 3 + BB * solX3 ** 2 + CC * solX3 + DD)
        #print('test', F1(fgf_f, Rt_f, Kb2_f, c_f))

    if False:
        plt.plot(fgf_list,outputs/fgf_list,'o')
        plt.xlabel(r'$F_t$')
        plt.ylabel(r'$F/F_t$')
        plt.ylim([0,1])
        plt.show()
    if True:
        plt.plot(fgf_list,frac_free_receptors,label='free receptors')
        plt.plot(fgf_list,frac_bound_receptors,label='singly-bound receptors')
        plt.plot(fgf_list,frac_doubly_bound_receptors,label='doubly-bound receptors')#plt.plot(fgf_list, fgf_list, 'x')

        plt.plot(fgf_list,model1,label='activity - Hill coefficient  0.44')#plt.plot(fgf_list, fgf_list, 'x')
        #plt.plot(fgf_list,predicted,label = 'predicted')#plt.plot(fgf_list,fgf_list**2*Rt_f/(1 + 2*fgf_list/c_f/Kb2_f + fgf_list**2/c_f/Kb2_f**2)/c_f/Kb2_f**2,'x')
        #plt.plot(fgf_list, 10**(2)+fgf_list, 'x')
        plt.xlabel(r'$F_t$')
        plt.legend(title='Fraction of :')

        plt.xscale('log')
        #plt.yscale('log')
        #plt.grid()
        plt.show()

    if False:
        print('fgf list', fgf_list)


        def myFR1(fgf, Rt, Kb2, c):
            myFt = F1(fgf, Rt, Kb2, c)
            print('my myFt', myFt)
            return Rt * myFt ** 2 / (myFt ** 2 + 2 * myFt * Kb2 + c * Kb2 ** 2)


        myFR1_list = [myFR1(fgf_list[kk], Rt_f, Kb2_f, c_f) for kk in range(len(fgf_list))]
        print('outputs', myFR1_list)
        plt.plot(fgf_list,myFR1_list,'o')
        plt.xscale('log')
        plt.show()






#[n = 0.59, K= 3.8, intercept = 0.66, xmax = 3.6]
mybounds_coop = [(0,10),(10**4,10**10),(0,0.001),(0.98,1)]
mybounds_vivo = [(.5, 50), (10 ** (-6), 10 ** 3), (0, 300), (10, 10 ** 3)] #got n = 2.6 !!
mybounds_vivo_subtle = [(10**(-9),10**6),(10**(-9),10),(0,10**9),(0,10),(0,.1),(0,10**9),(0,10)]
mybounds = [(10**(-8),10),(10**(-6),10**3),(10**(-8),10**3),(1.6,10)]
#['n1','K2','cst','erkmax','n2','K1']
mybounds_invitro = [(10**(-8),10),(10*(-6),10**3),(0,300),(10**(-8),10**6),(0.1,4),(10**3,10**9)]

mybounds_invitro_t = [(0,2),(0,10**(-6)),(3,4),(0.5,1)]#(10**6,10**12),

importdata = True
if importdata:
    halfhalf_fit = False
    if halfhalf_fit == True:
        halfhalf_xlsx = pd.ExcelFile('./data/dnEph3 half-half experimentspooled.xlsx')
        halfhalf = pd.read_excel(halfhalf_xlsx, 'averages') #to get columns names: halfhalf.columns


        for kk in range(5):
            erk = halfhalf.iloc[kk,[4,5,6,7]]
            erk = np.array(erk,dtype='float64')
            #popt, pcov = curve_fit(hill, surf_inputs, erk)
            def error_func2(p, Y=erk, X=surf_inputs):
                Ypredicted = [hill(surf_inputs[kk], p[0], p[1], p[2], p[3]) for kk in range(4)]
                diff = Y - Ypredicted
                return np.sum(np.square(diff))

            p0ann = dual_annealing(error_func2, bounds=mybounds)
            myp = p0ann.x
            print('para estim ',myp)
            plotexpe = plt.plot(surf_inputs,erk,'x')
            mycolor = plotexpe[0].get_color()
            erk_th = [hill(surf_inputs[kk],myp[0],myp[1],myp[2],myp[3]) for kk in range(4)]
            myorder = np.argsort(surf_inputs)
            surf_inputs_sorted = np.zeros(len(surf_inputs), dtype=float)
            for i in range(0, len(surf_inputs)):
                surf_inputs_sorted[i] = surf_inputs[myorder[i]]
            c = np.zeros(len(surf_inputs), dtype=float)
            for i in range(0, len(c)):
                c[i] = erk_th[myorder[i]]
            new_inputs = np.linspace(0,6000)

            erk_th = [hill(new_inputs[kk],myp[0],myp[1],myp[2],myp[3]) for kk in range(len(new_inputs))]
            plt.plot(new_inputs, erk_th, '--',c=mycolor,label='Hill coef %.2f '%(myp[0]))
            plt.title('ERK activation for a65, a66, a67 and a68 -- Half half manips no ephrine')
            plt.xlabel('surface')
            plt.ylabel('ERK activation')
            #plt.plot(surf_inputs[myorder],c,c = mycolor,label = 'Hill coef '+str(myp[0]))
            plt.legend()
        plt.show()

    exp190606 = True

    if exp190606 == True:
        cellbycell_xlsx = pd.ExcelFile('./data/cellbycell.xlsx')
        cellbycell = pd.read_excel(cellbycell_xlsx,'190606 ctrl-NVP')
        cellbycell_surf_inputs = cellbycell.iloc[0:48,9]# change here +2 for the second index corresponds to actual surfaces
        cellbycell_surf_inputs = cellbycell_surf_inputs.to_numpy()
        cellbycell_erk_no_eph_expe = cellbycell.iloc[0:48,13]
        cellbycell_erk_no_eph_expe = cellbycell_erk_no_eph_expe.to_numpy()

        select_cells = np.arange(0,48) # only right cells (0,48,2) -- only left cells  (0,48,2)+1 -- all cell (0,48,)

    exp190607 = False

    if exp190607 == True:
        cellbycell_xlsx = pd.ExcelFile('./data/cellbycell.xlsx')
        cellbycell = pd.read_excel(cellbycell_xlsx, '190607 ctrl-NVP')
        cellbycell_surf_inputs = cellbycell.iloc[0:40,8]#+2]  # change here +2 for the second index corresponds to actual surfaces
        cellbycell_surf_inputs = cellbycell_surf_inputs.to_numpy()
        cellbycell_erk_no_eph_expe = cellbycell.iloc[0:40, 12]
        cellbycell_erk_no_eph_expe = cellbycell_erk_no_eph_expe.to_numpy()

        select_cells = np.arange(0, 39)

    clare = True  # of Western blot no ephrin explants
    if clare == True:
        clare_xlsx = pd.ExcelFile('./data/clare.xlsx')
        clare = pd.read_excel(clare_xlsx, 'exp1')
        fgf_doses = clare.iloc[0:9,0]
        fgf_doses = fgf_doses.to_numpy()

        clare_erk_no_eph_expe = clare.iloc[0:9, 2:13]
        clare_erk_no_eph_expe  = clare_erk_no_eph_expe .to_numpy()
        temp = np.nanmean(clare_erk_no_eph_expe, axis = 1)
        fgf_doses = np.array([10**(-8),0.03,0.1,0.3,1,3,10,30,100,300]) #,
        # Western blot no ephrin explants
        temp = np.array([0.647767072,0.846152094,0.806518884,1.083480714,1.421676865,1.851103653,2.460453121,3.120250566,3.417251476,3.878459957])#

def_error_fct = True
if def_error_fct:
    #neg cooperativity
    def error_func_hill(p, Y, X):
        Ypredicted = [hill(X[kk], p[0], p[1], p[2], p[3]) for kk in range(len(X))]  #
        diff = Y - Ypredicted
        return np.sum(np.square(diff))

    #in_vivo
    def cellcell_error_func(p, Y=cellbycell_erk_no_eph_expe[select_cells], X=cellbycell_surf_inputs[select_cells]):
        Ypredicted = [hill(X[kk], p[0], p[1], p[2], p[3]) for kk in range(len(X))]
        diff = Y - Ypredicted
        return np.sum(np.square(diff))


    def cellcell_error_func_sublte(p, Y=cellbycell_erk_no_eph_expe[select_cells], X=cellbycell_surf_inputs[select_cells]):
        Ypredicted = [erk_act_subtle(X[kk],p[0],p[1],p[2],p[3],p[4],p[5],p[6]) for kk in range(len(X))] #K,c,Rt
        diff = Y - Ypredicted
        return np.sum(np.square(diff))
    #in vitro - fit with simple hill funct
    def error_func_erk_act_hill(p, Y=temp,X=fgf_doses):
        Ypredicted = [hill(X[kk],p[0], p[1], p[2], p[3]) for kk in range(len(X))] #
        diff = Y - Ypredicted
        return np.sum(np.square(diff))
    #fit with full model
    def error_func_erk_act(p,Y=temp,X=fgf_doses):
        Ypredicted = [erk_act(X[kk],p[0],p[1],p[2],p[3]) for kk in range(len(X))]
        diff = Y - Ypredicted
        return np.sum(np.square(diff))
    def error_func_erk_act_fit_curve(p,Y=temp,X=fgf_doses):
        Ypredicted = [erk_act(X[kk],p) for kk in range(len(X))]
        diff = Y - Ypredicted
        return np.sum(np.square(diff))


neg_coop_c_varied = False

if neg_coop_c_varied == True:
    c_list = np.array([10**(-3),1,10**6]) #1,100,1000
    fig = plt.figure(figsize=(10, 3), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    for jj in range(len(c_list)):
        c = c_list[jj]

        myY = np.array([fracact(myL[kk], 100, c) for kk in range(len(myL))])
        myY2 = np.array([fracact2(myL[kk], 100, c,100) for kk in range(len(myL))])

        def error_func_hill_c(p):
            return error_func_hill(p, Y=myY, X=myL)


        def error_func_hill_c2(p):
            return error_func_hill(p, Y=myY2, X=myL)


        print('Fit with a simple Hull function ')
        pann= dual_annealing(error_func_hill_c, mybounds_coop)
        pann2 = dual_annealing(error_func_hill_c2, mybounds_coop)
        print('my params hill',pann.x)
        y_fit1 = hill(myL, pann.x[0],pann.x[1],pann.x[2],pann.x[3])
        y_fit2 = hill(myL, pann2.x[0],pann2.x[1],pann2.x[2],pann2.x[3])

        l = ax.plot(myL, myY2)
        mycolor = l[0].get_color()
        #xx = Decimal(c)integer=True
        #     plt.title(r'Temporal evolution of Ras-GTP for sos = {}'.format(str_fmt(sostotlist[i])))
        print('test c',str_fmt(c))
        ax.plot(myL,y_fit2,'o',c=mycolor,label=r'c = {} and Hill coef is {}'.format(str_fmt(c),round_sig(pann2.x[0])))
        #plt.plot(myL,myY2,'x')
        xticks = ax.xaxis.get_major_ticks()
        xticks[1].label1.set_visible(False)
        xticks[3].label1.set_visible(False)
        xticks[5].label1.set_visible(False)
        xticks[7].label1.set_visible(False)
        #xticks[9].label1.set_visible(False)

        #ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
        ax.set_xlabel(r'$F_t$',fontsize=14)
        ax.set_ylabel(r'ERK activity $Y_{FR_2}$',fontsize=14)
        ax.set_title(r'ERK activity as a function of the total ligand concentration',fontsize=14)
        majors = [0, .5, 1]
        ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        plt.legend(fontsize=14)  # using a size in points
        #ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
        #ax.yaxis.set_minor_locator(ticker.MaxNLocator(40))
        #ax.text(0.0, 0.1, "MaxNLocator(n=4)", fontsize=14, transform=ax.transAxes)
        #ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=4))
    plt.legend(fontsize=12)
    ax.set_xscale('log')
    plt.savefig('anticoop.pdf')


neg_coop_Rt_varied = False

if neg_coop_Rt_varied == True:
    myL = np.logspace(-3, 18)

    Rt_list = np.array([10,10**2,10**3]) #1,100,1000
    fig = plt.figure(figsize=(10, 3), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    for jj in range(len(Rt_list)):
        Rt = Rt_list[jj]
        myY2 = np.array([fracact2(myL[kk], 1000, 0.01,Rt) for kk in range(len(myL))])

        def error_func_hill_c(p):
            return error_func_hill(p, Y=myY, X=myL)


        def error_func_hill_c2(p):
            return error_func_hill(p, Y=myY2, X=myL)


        print('Fit with a simple Hull function ')
        #pann= dual_annealing(error_func_hill_c, mybounds_coop)
        pann2 = dual_annealing(error_func_hill_c2, mybounds_coop,maxiter=10000)
        #y_fit1 = hill(myL, pann.x[0],pann.x[1],pann.x[2],pann.x[3])
        y_fit2 = hill(myL, pann2.x[0],pann2.x[1],pann2.x[2],pann2.x[3])

        l = ax.plot(myL, myY2)
        mycolor = l[0].get_color()
        #xx = Decimal(c)integer=True
        #     plt.title(r'Temporal evolution of Ras-GTP for sos = {}'.format(str_fmt(sostotlist[i])))
        ax.plot(myL,y_fit2,'o',c=mycolor,label=r'Rt = {} and Hill coef is {}'.format(str_fmt(Rt),round_sig(pann2.x[0])))
        #plt.plot(myL,myY2,'x')
        xticks = ax.xaxis.get_major_ticks()
        xticks[1].label1.set_visible(False)
        xticks[3].label1.set_visible(False)
        xticks[5].label1.set_visible(False)
        xticks[7].label1.set_visible(False)
        #xticks[9].label1.set_visible(False)

        #ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
        ax.set_xlabel(r'$F_t$',fontsize=14)
        ax.set_ylabel(r'ERK activity $Y_{FR_2}$',fontsize=14)
        ax.set_title(r'ERK activity as a function of the total ligand concentration',fontsize=14)
        majors = [0, .5, 1]
        ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        plt.legend(fontsize=14)  # using a size in points
        #ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
        #ax.yaxis.set_minor_locator(ticker.MaxNLocator(40))
        #ax.text(0.0, 0.1, "MaxNLocator(n=4)", fontsize=14, transform=ax.transAxes)
        #ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=4))
    plt.legend(fontsize=12)
    ax.set_xscale('log')
    plt.savefig('anticoopRt.pdf')

fitinvivo = False
if fitinvivo:
    p0ann = dual_annealing(cellcell_error_func, bounds=mybounds_vivo)
    myp = p0ann.x
    for kk in range(len(hill_params)):
        print('%s is %.2f'%(hill_params[kk],myp[kk]))
    input_th_vivo = np.linspace(0.01, .6, 100)
    output_th_vivo = [hill(input_th_vivo[kk], myp[0], myp[1], myp[2], myp[3]) for kk in range(len(input_th_vivo))]

    # print('Fit with a simple Hull function ')
    # popt, pcov = curve_fit(hill, cellbycell_surf_inputs, cellbycell_erk_no_eph_expe)
    # for kk in range(len(hill_params)):
    #     print('%s is %.2f' % (hill_params[kk], popt[kk]))
    # print('quality of the fit given by', pcov)

fitinvivo_subtle = False
if fitinvivo_subtle:
    p0ann = dual_annealing(cellcell_error_func_sublte, bounds=mybounds_vivo_subtle)
    myp = p0ann.x
    for kk in range(len(subtle_params)):
        print('%s is %.2f'%(subtle_params[kk],myp[kk]))
    input_th_vivo = np.linspace(0.01, .6, 100)
    output_th_vivo = [erk_act_subtle(input_th_vivo[kk], myp[0], myp[1], myp[2],myp[3],myp[4],myp[5],myp[6]) for kk in range(len(input_th_vivo))]

    # print('Fit with a simple Hull function ')
    # popt, pcov = curve_fit(hill, cellbycell_surf_inputs, cellbycell_erk_no_eph_expe)
    # for kk in range(len(hill_params)):
    #     print('%s is %.2f' % (hill_params[kk], popt[kk]))
    # print('quality of the fit given by', pcov)

fit_invitro = False
if fit_invitro:
    print('Fit with a simple Hull function ')
    popt, pcov = curve_fit(hill, fgf_doses, temp)
    for kk in range(len(hill_params)):
        print('%s is %.2f'%(hill_params[kk],popt[kk]))
    # print('quality of the fit given by', pcov)
    #
    # print('Fit with a simple Hull function - with dual_annealing ')
    # p0ann_hill = dual_annealing(error_func_erk_act_hill,bounds = mybounds)
    # print('second hill fit',p0ann_hill.x)
    #

    print('Fit with full model ')
    p0ann_invitro = dual_annealing(error_func_erk_act,bounds = mybounds_invitro_t)
    #myp_invitro, pcovbis = curve_fit(error_func_erk_act_fit_curve, fgf_doses, temp)
    myp_invitro = p0ann_invitro.x   #[n = 0.59, K= 3.8, intercept = 0.66, xmax = 3.6]
    #K1v = 10**12
    #myp_invitro =  [0.22,3.8/K1v**(.22),0.66,3.6,2.6,K1v]
    print('numb of param is',len(full_params_t))
    for kk in range(len(full_params_t)):
         print('%s is %.2f'%(full_params_t[kk],myp_invitro[kk]))

if False:
    p0ann_invitro = dual_annealing(error_func_erk_act,bounds = mybounds_invitro)
    myp_invitro = p0ann_invitro.x
    print('params in vitro ', myp_invitro)
    input_th = np.linspace(0.01, .6, 100)
    output_th = [hill(input_th[kk], myp[0], myp[1], myp[2], myp[3]) for kk in range(len(input_th))]

do_plots_vitro = False

if do_plots_vitro:
    #print('ratio is', (myp_invitro[2]**(myp_invitro[0])*myp_invitro[3]/myp_invitro[1])**2.6)
    th_inputs = np.logspace(-8,5,100)
    th_outputs = [hill(th_inputs[kk],popt[0],popt[1],popt[2],popt[3]) for kk in range(len(th_inputs))]
    th_outputs_full = [erk_act(th_inputs[kk],myp_invitro [0],myp_invitro [1],myp_invitro[2],myp_invitro[3]) for kk in range(len(th_inputs))]
    plt.plot(fgf_doses,temp,'x')
    #plt.plot(th_inputs,th_outputs)
    plt.plot(th_inputs,th_outputs_full,'o')
    plt.xscale('log')
    plt.show()

do_plots_vivo = False
if do_plots_vivo :
    markers = ['x', 'o', '^', 'v','+','*']
    mycolors = ['r','r','b','b','g','g','m','m'] # l & r cells
    # mycolors = ['r','b','g','m']# only l or only r cells
    number_of_repeats = 6 # for 190607 it is 5 -- for 06 it is 6
    for kk in range(number_of_repeats):
        indices_for_one_exp = np.array([0,1,2,3,4,5,6,7])+kk*8 # l and r : [0,1,2,3,4,5,6,7]  -- l or r : [0,2,4,6]
        for ll in range(len(indices_for_one_exp)):
            plt.scatter(cellbycell_surf_inputs[indices_for_one_exp[ll]],cellbycell_erk_no_eph_expe[indices_for_one_exp[ll]],marker = markers[kk],c=mycolors[ll])
    plt.plot(input_th_vivo,output_th_vivo)
    plt.title('nn = 2.6 dual_ann exp190606 (left and right cells) -  erk vs normalized surface of contact with FGFcells \n (nn,KK,interc,xmax) [(.5, 50), (10 ** (-6), 10 ** 3), (0, 300), (10, 10 ** 3)] ')
    plt.show()

do_plots_vivo_subtle = False
if do_plots_vivo_subtle:
    markers = ['x', 'o', '^', 'v','+','*']
    mycolors = ['r','r','b','b','g','g','m','m'] # l & r cells
    # mycolors = ['r','b','g','m']# only l or only r cells
    number_of_repeats = 6 # for 190607 it is 5 -- for 06 it is 6
    for kk in range(number_of_repeats):
        indices_for_one_exp = np.array([0,1,2,3,4,5,6,7])+kk*8 # l and r : [0,1,2,3,4,5,6,7]  -- l or r : [0,2,4,6]
        for ll in range(len(indices_for_one_exp)):
            plt.scatter(cellbycell_surf_inputs[indices_for_one_exp[ll]],cellbycell_erk_no_eph_expe[indices_for_one_exp[ll]],marker = markers[kk],c=mycolors[ll])
    plt.plot(input_th_vivo,output_th_vivo)
    plt.title('nn = 2.6 dual_ann exp190606 (left and right cells) -  erk vs normalized surface of contact with FGFcells \n (nn,KK,interc,xmax) [(.5, 50), (10 ** (-6), 10 ** 3), (0, 300), (10, 10 ** 3)] ')
    plt.show()

# if False:
#     colors = []
#     markers = ['x', 'o', '^', 'v']
#     indices = [0, 2, 4, 6]
#     for kk in range(6):
#         print(kk)
#         subsample = np.array([0,2,4,6])+kk*8#([0,2,4,6])[1,3,5,7]
#         subsample_l = np.array([1,3,5,7])+kk*8
#         print('subsample',subsample)
#         cellbycell_surf_inputs_t = cellbycell_surf_inputs[subsample]
#         cellbycell_erk_no_eph_expe_t = cellbycell_erk_no_eph_expe[subsample]
#         cellbycell_surf_inputs_tl = cellbycell_surf_inputs[subsample_l]
#         cellbycell_erk_no_eph_expe_tl = cellbycell_erk_no_eph_expe[subsample_l]
#         print(cellbycell_surf_inputs_t)
#         print(cellbycell_erk_no_eph_expe_t)
#         l = None
#         for aa in range(4):
#             jj = indices[aa] + 8*kk
#             if l is not None:
#                 plt.plot(np.array([cellbycell_surf_inputs_t[jj]]), np.array([cellbycell_erk_no_eph_expe_t[jj]]),
#                                  marker=markers[aa], color=l.get_color())
#             else:
#                 l, = plt.plot(np.array([cellbycell_surf_inputs_t[jj]]),np.array([cellbycell_erk_no_eph_expe_t[jj]]),marker=markers[aa])
#     plt.show()
#
# if False:
#
#     def cellcell_error_func(p, Y=cellbycell_erk_no_eph_expe, X=cellbycell_surf_inputs):
#         Ypredicted = [hill(X[kk],p[0],p[1],p[2],p[3]) for kk in range(len(X))]
#         diff = Y - Ypredicted
#         return np.sum(np.square(diff))
#
#     p0ann = dual_annealing(cellcell_error_func, bounds=mybounds)
#     myp = p0ann.x
#
#     print('para estim ', myp)
#     plotexpe = plt.plot(cellbycell_surf_inputs, cellbycell_erk_no_eph_expe, 'x')
#     mycolor = plotexpe[0].get_color()
#     #erk_th = [hill(surf_inputs[kk], myp[0], myp[1], myp[2], myp[3]) for kk in range(4)]
#     # myorder = np.argsort(surf_inputs)
#     # surf_inputs_sorted = np.zeros(len(surf_inputs), dtype=float)
#     # for i in range(0, len(surf_inputs)):
#     #     surf_inputs_sorted[i] = surf_inputs[myorder[i]]
#     # c = np.zeros(len(surf_inputs), dtype=float)
#     # for i in range(0, len(c)):
#     #     c[i] = erk_th[myorder[i]]
#     new_inputs = np.linspace(0, 6000)
#     erk_th = [hill(new_inputs[kk], myp[0], myp[1], myp[2], myp[3]) for kk in range(len(new_inputs))]
#     plt.plot(new_inputs, erk_th, '--', label='Hill coef %.2f ' % (myp[0]))
#     plt.title('ERK activation for a65, a66, a67 and a68 -- Half half manips no ephrine')
#     plt.xlabel('surface')
#     plt.ylabel('ERK activation')
#     # plt.plot(surf_inputs[myorder],c,c = mycolor,label = 'Hill coef '+str(myp[0]))
#     plt.legend()
#     plt.show()