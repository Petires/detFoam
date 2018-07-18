'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
'''

from SDToolbox import *
from CJ2 import *
import datetime
from scipy.integrate import ode

class cvoutput:
    def __init__(self,j):
        self.exo_time = 0
        self.ind_time = 0
        self.ind_time_10 = 0
        self.ind_time_90 = 0
        
        self.time = []
        self.T = []
        self.P = []
        self.species = []

def cv_CJ(plt_num, P1, T1, q, mech, fname):
    """

    cv_CJ.m
    Computes the time evolution of a constant volume explosion with
    shock (at CJ speed) heated reactants as initial conditions
    
    FUNCTION
    SYNTAX
    [gas] = cv_CJ(plt_num,P1,T1,q,mech,fname)
    
    INPUT
    plt_num = unused
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    fname = output file name (0 for no output file)
    
    OUTPUT
    gas = gas object at final equilibrium state

    """
    gas1 = Solution(mech);
    gas1.TPX = T1,P1,q; 
    r = gas1.density;
    [gas, cj_speed] = CJspeed2(P1, T1, q, mech);
    gas = PostShock_fr(cj_speed, P1, T1, q, mech);
    fname = fname + '_%d' % cj_speed
    b = 10000; j = gas.n_species;
    out = cvoutput(b,j)
    out = explosion(gas,fname,out);
    return [cj_speed, gas]

def cv_shk(U1, P1, T1, q, mech, fname):
    """

    cv_shk.m
    Computes the time evolution of a constant volume explosion with
    shock heated reactants as initial conditions
    
    FUNCTION
    SYNTAX
    [gas] = cv_shk(U1,P1,T1,q,mech,fname)
    
    INPUT
    U1 = shock speed (m/s)
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    fname = output file name (0 for no output file)
    
    OUTPUT
    gas = gas object at final equilibrium state

    """
    gas1 = Solution(mech);
    gas1.TPX = T1,P1,q; 
    r = gas1.density;
    gas = PostShock_fr(U1, P1, T1, q, mech);
    fname = fname + '_%d' % U1
    b = 10000; j = gas1.n_species;
    out = cvoutput(b,j)
    out = explosion(gas,fname,out);
    return gas

def uvsys(t,y,gas,mw):

    # Set the state of the gas, based on the current solution vector.
    gas.TDY = y[0], gas.density, y[1:len(y)]
    nsp = gas.n_species;

    # energy equation
    wdot = gas.net_production_rates
    tdot = - gas.T * gas_constant * np.dot((gas.standard_enthalpies_RT- [1.0]*nsp),wdot) / (gas.density*gas.cv_mass)
    # set up column vector for dydt
    dydt = [tdot]+[0.0]*nsp

    # species equations
    rrho = 1.0/gas.density
    for i in range (0,nsp):
        dydt[i+1] = rrho*mw[i]*wdot[i]
    return dydt

def explosion(gas,fig_num):

    ##################################################################################################
    #RUN CONSTANT VOLUME EXPLOSION
    ##################################################################################################
    temp_grad = []
    mw = gas.molecular_weights
    r = gas.density
    nsp = gas.n_species
    y = [0.0]*(nsp+1)
    y[0] = gas.T
    y[1:len(y)] = gas.Y
    t0 = 0.0
    t1 = 1.0
    output=cvoutput(nsp)
    out = ode(uvsys).set_integrator('vode', method='bdf', nsteps=1, with_jacobian=True, atol=1e-15, rtol=1e-9)
    out.set_initial_value(y, t0).set_f_params(gas, mw)
    out._integrator.iwork[2] = -1
    warnings.filterwarnings("ignore", category=UserWarning)
    output.time.append(out.t)
    output.T.append(out.y[0])
    output.species.append(out.y[1:len(out.y)])
    y=out.y[1:len(out.y)]
    gas.TDY = output.T[len(output.T)-1],r,y
    P=gas.P/one_atm
    output.P.append(P)
    cv = gas.cv_mass
    wdot = gas.net_production_rates
    hs = gas.standard_enthalpies_RT
    R = gas_constant
    wt = gas.mean_molecular_weight
    sumT = 0.0
    for z in range(0,nsp):
        w = mw[z]
        e = R*output.T[len(output.T)-1]*(hs[z]/w - 1/wt)
        wd = wdot[z]
        sumT = sumT + e*wd*w;
    temp_grad.append(-sumT/(r*cv))
    while out.t < t1:
        out.integrate(t1, step=True)
        output.time.append(out.t)
        output.T.append(out.y[0])
        output.species.append(out.y[1:len(out.y)])
        y=out.y[1:len(out.y)]
        gas.TDY = output.T[len(output.T)-1],r,y
        P=gas.P/one_atm
        output.P.append(P)
        cv = gas.cv_mass
        wdot = gas.net_production_rates
        hs = gas.standard_enthalpies_RT
        R = gas_constant
        wt = gas.mean_molecular_weight
        sumT = 0.0
        for z in range(0,nsp):
            w = mw[z]
            e = R*output.T[len(output.T)-1]*(hs[z]/w - 1/wt)
            wd = wdot[z]
            sumT = sumT + e*wd*w;
        temp_grad.append(-sumT/(r*cv))
    warnings.resetwarnings()
    del out
    ##################################################################################################
    #FIND INDUCTION TIME - MAXIMUM TEMPERATURE GRADIENT
    ##################################################################################################
    MAX = max(temp_grad)
    n = temp_grad.index(MAX)
    b = len(temp_grad)

    #FIND INDUCTION TIME - MAXIMUM TEMPERATURE GRADIENT
    k = 0; MAX = max(temp_grad); d = temp_grad[0]; HMPWt = zeros(2,float)
    if d == MAX:
        print 'Initial Gradient is Maximum - post shock temperature may be too low'
        return gas
    while d < MAX:
        k = k + 1; d = temp_grad[k];
    output.ind_time = output.time[k]; k1 = k; k = 0;
    MAX10 = 0.1*MAX; d = temp_grad[0];
    while(d < MAX10 and k < b-1):
        k = k + 1; d = temp_grad[k];
    if(k == b):
        print 'MAX10 may be incorrect - reached end of array'
    output.ind_time_10 = output.time[k]; k = 0;
    MAX90 = 0.9*MAX; d = temp_grad[0];
    while(d < MAX90 and k < b-1):
        k = k + 1; d = temp_grad[k];
    if(k == b-1):
        print 'MAX90 may be incorrect - reached end of array'
    output.ind_time_90 = output.time[k];
    output.exo_time = 0 
    
    if fig_num==0:
        return output

    else:
	k = 0; MAX = max(output.T); d = output.T[0];
	while d < MAX:
            k = k + 1; d = output.T[k];
	if output.time[k] == 0:
            maxt = output.ind_time*5;
        elif output.time[k] >= output.ind_time*50:
            maxt = output.ind_time*5;
	else:
            maxt = output.time[k] + 0.1*output.time[k];
    	mint = 0;
 	maxT = max(output.T)+0.1*min(output.T); minT = min(output.T)-0.1*min(output.T); 
 	maxP = max(output.P)+0.1*min(output.P); minP = min(output.P)-0.1*min(output.P); 
        maxpw = HMPWt[1] + 0.1*HMPWt[1]; minpw = HMPWt[0] - 0.1*HMPWt[0]; 
        maxTG = max(temp_grad) + 0.1*abs(max(temp_grad));
        minTG = min(temp_grad)-0.1*abs(min(temp_grad));
	d = datetime.date.today(); P = output.P[0]/one_atm;

        fn = str(fig_num) + '_CVprofile.plt';
        outputfile = file(fn, 'w');
        outputfile.write('# CONSTANT VOLUME PROFILES\n');
        outputfile.write('# CALCULATION RUN ON %s\n\n' % d);
        outputfile.write('# Maximum time calculated = %.4e\n' % max(output.time))
        outputfile.write('# t_min = %.4f, t_max = %.4e\n' % (mint, maxt))
        outputfile.write('# T_min = %.2f, T_max = %.2f\n' % (minT, maxT))
        outputfile.write('# P_min = %.2f, P_max = %.2f\n' % (minP, maxP))
        outputfile.write('# TG_min = %.2f, TG_max = %.2f\n' % (minTG, maxTG))
        outputfile.write('# THE OUTPUT DATA COLUMNS ARE:\n');
        outputfile.write('Variables = "Time", "Temperature", "Pressure", "temp_grad"\n');
        for i in range(b):
            outputfile.write('%.4E \t %.4E \t %.4E \t %.4E\t'% (output.time[i],output.T[i],output.P[i],temp_grad[i]));
            for j in range(nsp):
                outputfile.write('%.4E \t'% (output.species[i][j]))
            outputfile.write('\n');
    return output