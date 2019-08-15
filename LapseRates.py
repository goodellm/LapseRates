import numpy as np
import matplotlib.pyplot as plt

##################################################

#Set to True to create the corresponding plot
plotLapseRate = True
plotAdiabats = True
plotClimateSensitivity = True

#Set to True to save plots generated (saves to current directory)
savePlots = False

#Values for Lapse Rates:
Trange = range(200,311,1) #K
prange = range(20000,100000,10) #Pa

#Values for Adiabats:
ps = 1000.e2 #Pa
p_end = 10.e2 #Pa
deltap = 100.
Tsvals = range(200,320,10) #K

#Values for Climate Sensitivity:
Tsvals = range(200,320,10) #K
prad = 670.e2 #Pa
psfc = 1000.e2 #Pa

##################################################

L = 2.5e6 #J/kg
Ra = 287 #J/kg/K
Rc= 461 #J/kg/K
cpa= 1005 #J/kg/K
cpc= 1847 #J/kg/K
g = 9.8 #m/s^2
sigma = 5.67e-8 #W/m^2K^4

def pcsat(T):
    return 611.*(T/273.15)**(L/(Rc*273.15))

def rsat(T,pa):
    return Ra*pcsat(T)/(Rc*pa)

def gamma(T,p):      #lapse rate
    xa = L/(Ra*T)
    xc = L/(Rc*T)
    r = rsat(T,p) #assuming pa = p here
    num = xa*r
    denom = 1.+((cpc/cpa)+(xc-1.)*L/(cpa*T))*r
    result = (g*num/denom+g/denom)/cpa #approx cp as cpa
    return 1000.*result #multiply by 1000 to change units to K/km instead of K/m

def print_gamma(T, p):
    g = gamma(T, p)
    print "gamma(", T, ", ", p, ") = ", g
    print ""

if plotLapseRate:
    gammavals = []

    for p in prange:
        vals = []
        for T in Trange:
            vals.append(gamma(T,p))
        gammavals.append(vals)

    plt.figure(1)
    CS = plt.contour(Trange, prange, gammavals)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Moist Adiabatic Lapse Rate (K/km)')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (Pa)')
    if savePlots:
        plt.savefig('LapseRates.png',dpi=150)

def dTdp(T,p):
    return(T/p)*(Ra/g)*gamma(T,p)/1000.

if plotAdiabats or plotClimateSensitivity:
    Tval_list = []
    pval_list = []

    for Ts in Tsvals:
        pvals = []
        Tvals = []
        
        currentp = ps
        currentT = Ts
        
        while currentp >= p_end:
            pvals.append(currentp)
            Tvals.append(currentT)
            currentT = currentT-deltap*dTdp(currentT,currentp)
            currentp = currentp-deltap
        
        Tval_list.append(Tvals)
        pval_list.append(pvals)

if plotAdiabats:
    plt.figure(2)
    for i in range(len(Tval_list)):
        plt.plot(Tval_list[i],pval_list[i])
    
    plt.xlabel('Temperature(K)')
    plt.ylabel('Pressure(Pa)')
    plt.title('Moist Adiabatic Temperature Profiles')
    plt.gca().invert_yaxis()
    if savePlots:
        plt.savefig('MoistAdiabats.png',dpi=150)

def Trad(Ts,prad):
    ind = Tsvals.index(Ts)
    Tvals = Tval_list[ind]
    pvals = pval_list[ind]
    ind2 = pvals.index(prad)
    return Tvals[ind2]

def dTsdTrad(Ts,prad):
    deltaTs = 10
    deltaTrad = Trad(Ts+deltaTs,prad)-Trad(Ts,prad)
    return deltaTs/deltaTrad

def dTsdF(Ts,prad):
    delta = dTsdTrad(Ts,prad)
    T_rad = Trad(Ts,prad)
    return (1/(4*sigma*T_rad**3))*delta

def dry_climate_sensitivity(Ts,psfc,prad):
    T_rad = Ts/((psfc/prad)**(Ra/cpa))
    return (1/(4*sigma*T_rad**3))*((psfc/prad)**(Ra/cpa))

if plotClimateSensitivity:
    mcs = []
    dcs = []

    for Ts in Tsvals[0:-1]:
        moist_climate_sensitivity = dTsdF(Ts,prad)
        mcs.append(moist_climate_sensitivity)
        dcs.append(dry_climate_sensitivity(Ts,psfc,prad))
    
    plt.figure(3)
    plt.plot(Tsvals[0:-1],mcs,'b',label='Moist Adiabatic')
    plt.plot(Tsvals[0:-1],dcs,'r',label='Dry Adiabatic')
    plt.xlabel('Surface Temperature (K)')
    plt.ylabel('Climate Sensitivity dTs/dF')
    plt.legend(loc='best')
    plt.title('Climate Sensitivity')
    if savePlots:
        plt.savefig('ClimateSensitivity.png',dpi=150)

plt.show()
